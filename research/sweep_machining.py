#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sweep_machining.py — 扫描总分片数 → 加工时间/误差/耗时数据，输出 xlsx，支持断点重启。

工作流（单文件叶片）：
  1) auto-identify 找叶盆/叶背面 + split-blade 提取叶盆/叶背区间
  2) 对每个总分片数 total（枚举所有 (nSplitU, nSplitV) 组合）调 simple.exe 做固定分片拟合
  3) 读 *_params.txt + meta.json，用 compute_machining 计算侧铣/点铣加工时间
  4) 记录 (分片数, 拟合误差, 过切/欠切, 扭转角, 侧铣时间, 点铣时间, 提速比, 拟合/刀轨耗时)
结果实时写入 checkpoint JSON（断点重启用）和 xlsx（含 matplotlib 图表）。

用法:
  python sweep_machining.py <file> --out <result.xlsx> \
      [--checkpoint sweep_ckpt.json] [--tolerance 0.1] \
      [--total-min 4] [--total-max 100] [--workdir tmp]
"""
import os
import sys
import json
import time
import argparse
import subprocess
import shutil
import tempfile
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, str(Path(__file__).parent))
import compute_machining as cm
import toolpath_config as tcfg  # 连续刀轨参数配置（tool_len/tool_r/o_min/o_max/eps 等）
import toolpath_continuous as tc  # 连续刀轨（迭代条带法）
import machining_model as mm     # 加工参数/时间估计（工程公式）

PROJECT_DIR = Path(__file__).resolve().parent.parent
BUILD_EXE = PROJECT_DIR / "build" / "Release" / "simple.exe"

HEADERS = [
    "total", "num_combos", "num_patches",
    "min_error_mm", "max_error_mm", "mean_error_mm", "rms_error_mm",
    "max_overcut_mm", "max_undercut_mm",
    "max_twist_deg", "flank_err_mm",
    "flank_overcut_mm", "flank_overcut_mean_mm", "flank_overcut_std_mm",
    "flank_undercut_mm", "flank_undercut_mean_mm", "flank_undercut_std_mm",
    "flank_total_s", "point_total_s", "speedup",
    "fit_time_s", "toolpath_time_s",
]


def _key(total):
    return str(total)


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def run_fitting(file1, file2, outdir, tol, nu, nv, ps_ranges=None, face_idx=0):
    cmd = [str(BUILD_EXE), file1, file2, "--mode", "ruled",
           "--outdir", outdir, "--tolerance", str(tol),
           "--nsplit-u", str(nu), "--nsplit-v", str(nv),
           "--no-refine", "--balance-edges", "--max-cells", "200000",
           "--blend", "--blend-width", "0.05", "--continuity", "3"]
    if ps_ranges:
        for k, (dir_, ranges) in enumerate(ps_ranges, start=1):
            if not ranges:
                continue
            rk = 'v' if dir_ == 'V' else 'u'
            us, ue = ranges[0]   # 取第一个叶盆/叶背区间
            cmd += [f"--face-idx{k}", str(face_idx), f"--{rk}-range{k}", f"{us},{ue}"]
    t0 = time.time()
    r = subprocess.run(cmd, cwd=str(PROJECT_DIR), capture_output=True,
                       text=True, encoding="utf-8", errors="replace", timeout=900)
    fit_sec = time.time() - t0
    if r.returncode != 0:
        err = (r.stderr or "").strip().splitlines()
        if err:
            log(f"  [fitting stderr] {'; '.join(err[-3:])}")
    return r.returncode, fit_sec


def run_exe_json(args_list):
    """运行 simple.exe 并从 stdout 解析 JSON（跳过前面的日志行）。"""
    r = subprocess.run([str(BUILD_EXE)] + args_list, cwd=str(PROJECT_DIR),
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=300)
    s = r.stdout.strip()
    i = s.find('{')
    if i >= 0:
        s = s[i:]
    try:
        return json.loads(s)
    except Exception:
        return None


def auto_identify(filepath):
    """auto-identify 找叶片面，返回 (pressureIndex, suctionIndex)。"""
    info = run_exe_json([filepath, "--mode", "auto-identify"])
    if not info or not info.get('success'):
        return 0, 0
    return info.get('pressureIndex', 0), info.get('suctionIndex', 0)


def split_blade(filepath, face_idx):
    """split-blade 一个面，返回 (dir, [regions])。"""
    info = run_exe_json([filepath, "--mode", "split-blade",
                         "--face-idx", str(face_idx)])
    if not info or not info.get('success'):
        return None, []
    return info.get('dir', 'V'), info.get('regions', [])


def identify_ps(regions):
    """从 split-blade 结果里挑出叶盆/叶背区间，返回 (pressure, suction) 两个 (u_start, u_end)。"""
    pressure = [r for r in regions if r.get('label') == 'pressure']
    suction = [r for r in regions if r.get('label') == 'suction']
    nonedge = [r for r in regions if r.get('label') != 'edge']
    if not pressure and not suction:
        # 标签不全时，按宽度取两个非 edge 区间
        nonedge.sort(key=lambda r: r.get('uEnd', 0) - r.get('uStart', 0), reverse=True)
        if len(nonedge) >= 2:
            return (nonedge[0]['uStart'], nonedge[0]['uEnd']), \
                   (nonedge[1]['uStart'], nonedge[1]['uEnd'])
        if nonedge:
            return (nonedge[0]['uStart'], nonedge[0]['uEnd']), None
        return None, None
    p = (pressure[0]['uStart'], pressure[0]['uEnd']) if pressure else None
    s = (suction[0]['uStart'], suction[0]['uEnd']) if suction else None
    return p, s


def read_error_stats(outdir):
    """从 meta.json 读取误差分布统计（max/mean/Q95/rms + 过切/欠切）。"""
    meta = os.path.join(outdir, "meta.json")
    if not os.path.exists(meta):
        return {}
    try:
        with open(meta, encoding="utf-8") as f:
            d = json.load(f)
    except Exception:
        return {}
    errs = []
    rmses = []
    overcuts = []
    undercuts = []
    for s in d.get("surfaces", []):
        for c in s.get("cells", []):
            errs.append(c.get("maxErr", 0.0))
            rmses.append(c.get("rmsErr", 0.0))
            overcuts.append(c.get("maxOvercut", 0.0))
            undercuts.append(c.get("maxUndercut", 0.0))
    if not errs:
        return {}
    errs.sort()
    n = len(errs)
    return {
        "max_error_mm": round(errs[-1], 6),
        "mean_error_mm": round(sum(errs) / n, 6),
        "q95_error_mm": round(errs[min(n - 1, int(n * 0.95))], 6),
        "rms_error_mm": round((sum(e * e for e in rmses) / n) ** 0.5, 6),
        "max_overcut_mm": round(max(overcuts, default=0.0), 6),
        "max_undercut_mm": round(max(undercuts, default=0.0), 6),
    }


def collect(outdir, args, fit_sec=0.0):
    """与 UI 一致：用连续刀轨（toolpath_continuous）计算侧铣/点铣加工时间。"""
    t0 = time.time()
    models = sorted(fn for fn in os.listdir(outdir) if fn.endswith("_surface_model.json"))
    if not models:
        return None
    flank_time = 0.0
    point_time = 0.0
    num_strips = 0
    err_mean = 0.0
    err_overcut = 0.0
    err_undercut = 0.0
    err_overcut_mean = 0.0
    err_overcut_std = 0.0
    err_undercut_mean = 0.0
    err_undercut_std = 0.0
    for mp in models:
        try:
            res = tc.plan_toolpaths(os.path.join(outdir, mp), eps=args.eps)
        except Exception:
            res = None
        if res is None or not res.get("machining"):
            continue
        flank_time += res["machining"]["flank_time_s"]
        point_time += res["machining"]["point_time_s"]
        num_strips += res["N"]
        err_mean = max(err_mean, res["err_mean"])
        err_overcut = max(err_overcut, res["machining"].get("err_overcut_mm", 0.0))
        err_undercut = max(err_undercut, res["machining"].get("err_undercut_mm", 0.0))
        err_overcut_mean = max(err_overcut_mean, res["machining"].get("err_overcut_mean_mm", 0.0))
        err_overcut_std = max(err_overcut_std, res["machining"].get("err_overcut_std_mm", 0.0))
        err_undercut_mean = max(err_undercut_mean, res["machining"].get("err_undercut_mean_mm", 0.0))
        err_undercut_std = max(err_undercut_std, res["machining"].get("err_undercut_std_mm", 0.0))
    toolpath_sec = time.time() - t0
    speedup = point_time / flank_time if flank_time > 0 else 0.0
    # 拟合格数 + 最大扭转角（meta.json）
    num_patches, max_twist = _meta_stats(outdir)
    rec = {
        "num_patches": num_patches if num_patches else num_strips,
        "max_twist_deg": round(max_twist, 3),
        "flank_cut_s": round(flank_time, 2),
        "flank_total_s": round(flank_time, 2),
        "point_cut_s": round(point_time, 2),
        "point_total_s": round(point_time, 2),
        "speedup": round(speedup, 3),
        "flank_err_mm": round(err_mean, 4),
        "flank_overcut_mm": round(err_overcut, 4),
        "flank_undercut_mm": round(err_undercut, 4),
        "flank_overcut_mean_mm": round(err_overcut_mean, 4),
        "flank_overcut_std_mm": round(err_overcut_std, 4),
        "flank_undercut_mean_mm": round(err_undercut_mean, 4),
        "flank_undercut_std_mm": round(err_undercut_std, 4),
        "fit_time_s": round(fit_sec, 4),
        "toolpath_time_s": round(toolpath_sec, 4),
    }
    rec.update(read_error_stats(outdir))
    return rec


def _meta_stats(outdir):
    """从 meta.json 统计拟合直纹面格数（nRows×nCols 求和）与最大扭转角。"""
    meta_path = os.path.join(outdir, "meta.json")
    if not os.path.exists(meta_path):
        return 0, 0.0
    try:
        with open(meta_path, encoding="utf-8") as f:
            meta = json.load(f)
    except Exception:
        return 0, 0.0
    total = 0
    max_twist = 0.0
    for s in meta.get("surfaces", []):
        total += int(s.get("nRows", 0)) * int(s.get("nCols", 0))
        max_twist = max(max_twist, float(s.get("maxTwist", 0.0)))
    return total, max_twist


def _run_combo(task):
    """并行 worker：对单个 (nu,nv) 组合跑拟合 + 刀轨计算。返回 (nu, nv, rec 或 None)。"""
    file, workdir_base, total, nu, nv, tol, ps_ranges, face_idx, mach_args = task
    outdir = os.path.join(workdir_base, f"t{total}_u{nu}_v{nv}")
    if os.path.isdir(outdir):
        shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir, exist_ok=True)
    try:
        rc, fit_sec = run_fitting(file, file, outdir, tol, nu, nv, ps_ranges, face_idx)
    except Exception:
        return (nu, nv, None)
    if rc != 0:
        return (nu, nv, None)
    rec = collect(outdir, mach_args, fit_sec)
    return (nu, nv, rec)


def write_xlsx(rows, out_path, eps=None):
    import openpyxl
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "machining_sweep"
    ws.append(HEADERS)
    for r in rows:
        ws.append([r.get(h) for h in HEADERS])
    if eps is not None:
        meta = wb.create_sheet("meta")
        meta["A1"] = "eps"
        meta["B1"] = float(eps)
    wb.save(out_path)


def add_charts(out_path, eps=None):
    """用 matplotlib 生成 PNG 图表并嵌入 xlsx（样式完全可控，无图例重叠问题）。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for _f in ["Microsoft YaHei", "SimHei", "SimSun"]:
        if any(f.name == _f for f in font_manager.fontManager.ttflist):
            plt.rcParams["font.sans-serif"] = [_f]
            break
    plt.rcParams["axes.unicode_minus"] = False
    import openpyxl
    from openpyxl.drawing.image import Image as XLImage

    try:
        wb = openpyxl.load_workbook(out_path)
    except Exception:
        return
    # 用数据表 machining_sweep（不能用 wb.active，因为 charts/meta 表可能被设为 active）
    ws = wb["machining_sweep"]
    if ws.max_row < 2:
        return

    # 点铣基线误差：优先用传入的 eps，否则从 xlsx 的 meta 表读，再回退到 config.scallop
    if eps is None:
        try:
            eps = float(wb["meta"]["B1"].value)
        except Exception:
            mcfg = cm.load_config()
            eps = mcfg.get("scallop", 0.1)
    eps = float(eps)

    headers = [c.value for c in ws[1]]
    rows = [list(r) for r in ws.iter_rows(min_row=2, values_only=True)]
    try:
        x_idx = headers.index("total")
        speedup_idx = headers.index("speedup")
        err_cols = [headers.index(x) for x in
                    ("min_error_mm", "mean_error_mm", "rms_error_mm")]
    except ValueError:
        return

    data = sorted((r for r in rows if r[x_idx] is not None),
                  key=lambda r: r[x_idx])
    totals = [r[x_idx] for r in data]
    speedup = [r[speedup_idx] for r in data]

    tmp = os.path.join(tempfile.gettempdir(), "sweep_charts")
    os.makedirs(tmp, exist_ok=True)

    # 图1：提速 vs 总分片数（单色折线，无图例）
    fig, ax = plt.subplots(figsize=(8, 4), dpi=120)
    ax.plot(totals, speedup, color="#1f77b4", linewidth=1.6)
    ax.set_xlabel("Total Patches")
    ax.set_ylabel("Speedup")
    ax.set_title("Speedup vs Total Patches")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    p1 = os.path.join(tmp, "speedup.png")
    fig.savefig(p1)
    plt.close(fig)

    # 图2：误差 vs 总分片数（统一取最优组合的误差：最大/平均/RMS）
    fig, ax = plt.subplots(figsize=(8, 4), dpi=120)
    err_labels = ["max error (best)", "mean error (best)", "rms error (best)"]
    for ci, lab in zip(err_cols, err_labels):
        ys = [r[ci] for r in data]
        ax.plot(totals, ys, linewidth=1.3, label=lab)
    ax.set_xlabel("Total Patches")
    ax.set_ylabel("Error (mm)")
    ax.set_title("Error vs Total Patches (best split per total)")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    p2 = os.path.join(tmp, "error.png")
    fig.savefig(p2)
    plt.close(fig)

    # 图3：提速 vs 平均误差（侧铣 vs 点铣基线）
    fig, ax = plt.subplots(figsize=(6, 4), dpi=120)
    mean_err = [r[headers.index("mean_error_mm")] for r in data]
    ax.scatter(mean_err, speedup, color="#1f77b4", s=18, label="侧铣(直纹面拟合)")
    ax.scatter([eps], [1.0], color="red", marker="*", s=200, zorder=5,
               label=f"点铣基线(残留{eps}mm, 提速1)")
    ax.axhline(1.0, color="red", linestyle="--", linewidth=1, alpha=0.6)
    ax.axvline(eps, color="green", linestyle=":", linewidth=1, alpha=0.6)
    ax.set_xlabel("Mean Error (mm)")
    ax.set_ylabel("Speedup")
    ax.set_title("Speedup vs Mean Error (flank vs point)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    p3 = os.path.join(tmp, "speedup_vs_error.png")
    fig.savefig(p3)
    plt.close(fig)

    # 图4：过切/欠切 vs 总分片数（符号误差对比维度）
    try:
        oc_idx = headers.index("max_overcut_mm")
        uc_idx = headers.index("max_undercut_mm")
        fig, ax = plt.subplots(figsize=(8, 4), dpi=120)
        ax.plot(totals, [r[oc_idx] for r in data], linewidth=1.4,
                color="#d62728", label="Max Over-cut (过切)")
        ax.plot(totals, [r[uc_idx] for r in data], linewidth=1.4,
                color="#2ca02c", label="Max Under-cut (欠切)")
        ax.set_xlabel("Total Patches")
        ax.set_ylabel("Signed error (mm)")
        ax.set_title("Over-cut / Under-cut vs Total Patches")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        p4 = os.path.join(tmp, "overcut_undercut.png")
        fig.savefig(p4)
        plt.close(fig)
    except (ValueError, IndexError):
        p4 = None

    # 图5：拟合 / 刀轨计算耗时 vs 总分片数
    try:
        ft_idx = headers.index("fit_time_s")
        tt_idx = headers.index("toolpath_time_s")
        fig, ax = plt.subplots(figsize=(8, 4), dpi=120)
        ax.plot(totals, [r[ft_idx] for r in data], linewidth=1.4,
                color="#1f77b4", label="Fitting time (分区拟合)")
        ax.plot(totals, [r[tt_idx] for r in data], linewidth=1.4,
                color="#ff7f0e", label="Toolpath time (刀轨计算)")
        ax.set_xlabel("Total Patches")
        ax.set_ylabel("Time (s)")
        ax.set_title("Runtime vs Total Patches")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        p5 = os.path.join(tmp, "runtime.png")
        fig.savefig(p5)
        plt.close(fig)
    except (ValueError, IndexError):
        p5 = None

    # 图6：加工过切/欠切 vs 总分片数（刀轨符号误差，随 eps 变化）
    try:
        foc_idx = headers.index("flank_overcut_mm")
        fuc_idx = headers.index("flank_undercut_mm")
        fig, ax = plt.subplots(figsize=(8, 4), dpi=120)
        ax.plot(totals, [r[foc_idx] for r in data], linewidth=1.4,
                color="#d62728", label="Machining Over-cut (加工过切)")
        ax.plot(totals, [r[fuc_idx] for r in data], linewidth=1.4,
                color="#2ca02c", label="Machining Under-cut (加工欠切)")
        ax.set_xlabel("Total Patches")
        ax.set_ylabel("Signed toolpath error (mm)")
        ax.set_title("Machining Over-cut / Under-cut vs Total Patches")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        p6 = os.path.join(tmp, "machining_overcut_undercut.png")
        fig.savefig(p6)
        plt.close(fig)
    except (ValueError, IndexError):
        p6 = None

    if "charts" in wb.sheetnames:
        del wb["charts"]
    cs = wb.create_sheet("charts")
    cs.add_image(XLImage(p1), "A1")
    cs.add_image(XLImage(p2), "A22")
    cs.add_image(XLImage(p3), "H1")
    if p4:
        cs.add_image(XLImage(p4), "A44")
    if p5:
        cs.add_image(XLImage(p5), "H22")
    if p6:
        cs.add_image(XLImage(p6), "A66")
    # 把 eps 写回 meta 表，保证下次 --charts-only 不带 --eps 也能读对
    if eps is not None:
        if "meta" not in wb.sheetnames:
            wb.create_sheet("meta")
        wb["meta"]["A1"] = "eps"
        wb["meta"]["B1"] = float(eps)
    wb.save(out_path)


def load_xlsx_rows(path):
    """读取已有 xlsx 的行（用于断点追加，不丢失历史数据）。"""
    import openpyxl
    if not os.path.exists(path):
        return []
    try:
        wb = openpyxl.load_workbook(path)
        ws = wb.active
        headers = None
        rows = []
        for row in ws.iter_rows(values_only=True):
            if headers is None:
                headers = list(row)
                continue
            d = dict(zip(headers, row))
            if d.get("total") is not None:
                rows.append(d)
        return rows
    except Exception:
        return []


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        try:
            sys.stdout.reconfigure(encoding='utf-8', errors='replace')
            sys.stderr.reconfigure(encoding='utf-8', errors='replace')
        except Exception:
            pass
    ap = argparse.ArgumentParser(description="单文件叶片：auto-identify→split-blade→拟合扫描(总分片数范围)，输出 xlsx")
    ap.add_argument("file", nargs="?", help="叶片文件（单文件，如 blade.igs；--charts-only 时可不填）")
    ap.add_argument("--out", required=True, help="输出 xlsx 路径")
    ap.add_argument("--checkpoint", default=None, help="断点检查点 JSON 路径")
    ap.add_argument("--tolerance", type=float, default=0.1, help="拟合容差 mm")
    ap.add_argument("--total-min", type=int, default=4, help="总分片数下限（含）")
    ap.add_argument("--total-max", type=int, default=100, help="总分片数上限（含）")
    ap.add_argument("--workdir", default=None, help="临时工作目录（默认系统临时目录）")
    ap.add_argument("--workers", type=int, default=None, help="并行进程数（默认 CPU 核数）")
    ap.add_argument("--eps", type=float, default=None, help="侧铣可接受误差（默认取 toolpath_config 的 eps）")
    ap.add_argument("--charts-only", action="store_true", help="仅从 --out 的 xlsx 重新生成图表（不重跑拟合）")
    args = ap.parse_args()

    # 仅重新生成图表
    if args.charts_only:
        add_charts(args.out, args.eps)
        log(f"已从 {args.out} 重新生成图表")
        sys.exit(0)

    if not args.file:
        log("[ERROR] 需要叶片文件（或使用 --charts-only 仅重生成图表）")
        sys.exit(1)

    # 关键：把叶片文件转成绝对路径。子进程 simple.exe 以 PROJECT_DIR 为 cwd 运行，
    # 相对路径若按 PROJECT_DIR 解析会失效（例如在 research 下用 "..\Blade.igs"）。
    args.file = os.path.abspath(args.file)

    if not BUILD_EXE.exists():
        log(f"[ERROR] simple.exe not found: {BUILD_EXE}")
        sys.exit(1)

    # 模拟 UI 手动流程：auto-identify → split-blade → identify pressure/suction，缓存结果
    log("auto-identify...")
    pidx, sidx = auto_identify(args.file)
    log(f"  pressureIndex={pidx} suctionIndex={sidx}")

    split_face_idx = pidx if pidx >= 0 else 0
    log(f"split-blade face {split_face_idx}...")
    split_dir, regions = split_blade(args.file, split_face_idx)
    log(f"  dir={split_dir} regions={len(regions)}")
    for reg in regions:
        log(f"    {reg.get('label')} V[{reg.get('uStart')},{reg.get('uEnd')}]")

    pressure, suction = identify_ps(regions)
    log(f"  pressure={pressure} suction={suction}")

    ps_ranges = None
    if pressure or suction:
        ps_ranges = [
            (split_dir, [pressure] if pressure else []),
            (split_dir, [suction] if suction else []),
        ]
    else:
        log("  未识别到叶盆/叶背区间，将拟合整面")

    tcfg_ = tcfg.load_config()
    mcfg = cm.load_config()
    eps = args.eps if args.eps is not None else tcfg_.get("eps", 0.1)
    mach_args = argparse.Namespace(
        feed=mcfg.get("feed", 500.0), tool_r=mcfg.get("tool_r", 5.0),
        ball_r=mcfg.get("ball_r", 5.0), scallop=mcfg.get("scallop", 0.1),
        twist_limit=mcfg.get("twist_limit", 2.0),
        taper_angle=mcfg.get("taper_angle", 3.0),
        overhead=mcfg.get("overhead", 4.0),
        point_overhead=mcfg.get("point_overhead", 10.0),
        eps=eps)

    # 断点重启：合并「已有 xlsx」+「checkpoint」，按 total 去重
    def key_of(r):
        if r.get("key"):
            return r["key"]
        return _key(r.get("total"))

    merged = {}
    combo_ckpt = {}  # {(total, nu, nv): rec}，用于组合级断点恢复
    for r in load_xlsx_rows(args.out):
        r.setdefault("key", key_of(r))
        merged[r["key"]] = r
    if args.checkpoint and os.path.exists(args.checkpoint):
        try:
            with open(args.checkpoint, encoding="utf-8") as f:
                ck = json.load(f)
            for r in ck.get("results", []):
                r.setdefault("key", key_of(r))
                merged[r["key"]] = r
            for k, v in ck.get("combos", {}).items():
                combo_ckpt[k] = v
        except Exception:
            pass
    results = list(merged.values())
    done_keys = {r["key"] for r in results if r.get("key")}

    workdir_base = args.workdir or tempfile.mkdtemp(prefix="sweep_")
    os.makedirs(workdir_base, exist_ok=True)

    totals = list(range(args.total_min, args.total_max + 1))
    log(f"扫描总分片数 {args.total_min}~{args.total_max}（{len(totals)} 个），"
        f"每个 total 枚举所有乘法组合并取平均")
    log(f"simple.exe: {BUILD_EXE}")

    def enumerate_splits(total):
        combos = []
        for nu in range(1, total + 1):
            if total % nu == 0:
                combos.append((nu, total // nu))
        return combos

    def save_checkpoint():
        if not args.checkpoint:
            return
        try:
            with open(args.checkpoint, "w", encoding="utf-8") as f:
                json.dump({"results": results, "combos": combo_ckpt}, f,
                          ensure_ascii=False, indent=1)
        except Exception as e:
            log(f"  [warn] 写 checkpoint 失败: {e}")

    def save_xlsx():
        try:
            write_xlsx(results, args.out, eps)
        except Exception as e:
            log(f"  [warn] 写 xlsx 失败: {e}")

    def aggregate_rows(combo_rows):
        if not combo_rows:
            return None
        # 统一取「最优组合」（误差最小的切分）：误差/过切/欠切/扭转/侧铣误差/时间/提速
        # 全部来自同一个组合，与提速比口径一致，避免把最差组合的误差和最优组合的提速混在一起。
        best = min(combo_rows, key=lambda r: r["max_error_mm"])
        return {
            "num_patches": best["num_patches"],
            "min_error_mm": round(best["max_error_mm"], 4),
            "max_error_mm": round(best["max_error_mm"], 4),
            "mean_error_mm": round(best["mean_error_mm"], 4),
            "rms_error_mm": round(best["rms_error_mm"], 4),
            "max_overcut_mm": round(best["max_overcut_mm"], 4),
            "max_undercut_mm": round(best["max_undercut_mm"], 4),
            "max_twist_deg": round(best["max_twist_deg"], 3),
            "flank_err_mm": round(best["flank_err_mm"], 4),
            "flank_overcut_mm": round(best.get("flank_overcut_mm", 0.0), 4),
            "flank_overcut_mean_mm": round(best.get("flank_overcut_mean_mm", 0.0), 4),
            "flank_overcut_std_mm": round(best.get("flank_overcut_std_mm", 0.0), 4),
            "flank_undercut_mm": round(best.get("flank_undercut_mm", 0.0), 4),
            "flank_undercut_mean_mm": round(best.get("flank_undercut_mean_mm", 0.0), 4),
            "flank_undercut_std_mm": round(best.get("flank_undercut_std_mm", 0.0), 4),
            "flank_total_s": best["flank_total_s"],
            "point_total_s": best["point_total_s"],
            "speedup": best["speedup"],
            "fit_time_s": round(best["fit_time_s"], 4),
            "toolpath_time_s": round(best["toolpath_time_s"], 4),
        }

    from concurrent.futures import as_completed

    def _combo_key(total, nu, nv):
        return f"{total}_{nu}_{nv}"

    i = 0
    for total in totals:
        key = _key(total)
        i += 1
        if key in done_keys:
            log(f"[{i}/{len(totals)}] skip (已完成): total={total}")
            continue
        combos = enumerate_splits(total)
        # 断点恢复：已缓存的组合直接复用，其余进队列重跑
        combo_rows = []
        pending = []
        for (nu, nv) in combos:
            ck = _combo_key(total, nu, nv)
            if ck in combo_ckpt:
                combo_rows.append(combo_ckpt[ck])
            else:
                pending.append((nu, nv))
        log(f"[{i}/{len(totals)}] total={total} → {len(combos)} 组合 "
            f"(已缓存 {len(combo_rows)}，待跑 {len(pending)}) {combos}")
        if pending:
            tasks = [(args.file, workdir_base, total, nu, nv, args.tolerance,
                      ps_ranges, split_face_idx, mach_args) for (nu, nv) in pending]
            with ProcessPoolExecutor(max_workers=args.workers) as ex:
                futures = [ex.submit(_run_combo, t) for t in tasks]
                for fut in as_completed(futures):
                    nu, nv, rec = fut.result()
                    if rec is None:
                        log(f"    [{nu}x{nv}] 拟合失败/未解析")
                        continue
                    combo_ckpt[_combo_key(total, nu, nv)] = rec
                    combo_rows.append(rec)
                    save_checkpoint()  # 每个组合完成即落盘，中断后可续跑
        if not combo_rows:
            log(f"  total={total} 所有组合失败，跳过")
            continue
        agg = aggregate_rows(combo_rows)
        agg["key"] = key
        agg["total"] = total
        agg["num_combos"] = len(combo_rows)
        results.append(agg)
        log(f"  → minErr={agg['min_error_mm']}mm maxErr={agg['max_error_mm']}mm "
            f"meanErr={agg['mean_error_mm']}mm rmsErr={agg['rms_error_mm']}mm "
            f"最优提速={agg['speedup']}x")
        save_xlsx()
        save_checkpoint()

    add_charts(args.out, eps)
    log(f"完成。共 {len(results)} 条结果 → {args.out}")
    if args.checkpoint:
        log(f"断点文件 → {args.checkpoint}")


if __name__ == "__main__":
    main()
