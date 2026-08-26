#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
nx_export.py — 把直纹面拟合结果导出为 Siemens NX 原生 .prt。

读取 simple.exe 生成的 `<blade>_directrices.json`（每个直纹格的两条优化准线
采样点），通过 NX Open 逐格重建「直纹面」特征（Insert→Mesh Surface→Ruled，
即 RuledBuilder），保存为 .prt。直纹语义被保留，NX CAM 可用侧刃铣
（Multi-Axis Variable Contour + Tanto Fan / Flank Milling）。

★★ 运行方式（脚本必须在 NX 进程内运行，不是外部进程）★★
  1. 打开 NX，随便新建或打开一个 Part。
  2. File → Preferences → User Interface → 把 Journal Language 改成 Python。
  3. 顶部 Ribbon 右键 → 勾选 Developer（开发）选项卡。
  4. Developer → Journal → 选择本脚本运行（或 File→Execute→NX Open 选择 .py）。
  5. 脚本读完 JSON 后自动建样条 + 直纹面，最后 File→Save As 保存 .prt。

依赖：本机已安装 Siemens NX（脚本只在 NX 内运行，无外部 pip 依赖）。

说明：
  - --resample 把每条准线均匀下采样到 N 个点，减少特征数提速。
  - 默认把叶片平移到原点（NX CAM 默认 MCS 在原点，叶片太远刀轨全 0）。
  - 默认按准线长度降序（长=易失败）优先建，减少个别失败。
"""

import os
import sys
import json
import math

# 在 NX 内运行时可导入 NXOpen；若在外部 python 运行（--dry-run）则跳过。
HAS_NXOPEN = False
try:
    import NXOpen
    import NXOpen.Features
    HAS_NXOPEN = True
except Exception:
    pass


# --------------------------------------------------------------------------
# 纯数据处理（与 CATIA 版一致，无 NX 依赖）
# --------------------------------------------------------------------------
def resample(pts, target):
    n = len(pts)
    if n <= 2 or target >= n:
        return pts
    out = []
    for i in range(target):
        idx = int(round(i * (n - 1) / float(target - 1)))
        out.append(pts[min(idx, n - 1)])
    dedup = []
    for p in out:
        if not dedup or p != dedup[-1]:
            dedup.append(p)
    return dedup if len(dedup) >= 2 else pts


def directrix_length(c):
    c0 = c["c0"]
    total = 0.0
    for i in range(1, len(c0)):
        a, b = c0[i - 1], c0[i]
        dx = a[0] - b[0]; dy = a[1] - b[1]; dz = a[2] - b[2]
        total += (dx * dx + dy * dy + dz * dz) ** 0.5
    return total


def bbox_center(cells):
    mn = [1e30, 1e30, 1e30]
    mx = [-1e30, -1e30, -1e30]
    for c in cells:
        for p in c["c0"] + c["c1"]:
            for k in range(3):
                if p[k] < mn[k]: mn[k] = p[k]
                if p[k] > mx[k]: mx[k] = p[k]
    return [(mn[0] + mx[0]) / 2, (mn[1] + mx[1]) / 2, (mn[2] + mx[2]) / 2]


def translate_cells(cells, center):
    out = []
    for c in cells:
        c0 = [[p[0] - center[0], p[1] - center[1], p[2] - center[2]] for p in c["c0"]]
        c1 = [[p[0] - center[0], p[1] - center[1], p[2] - center[2]] for p in c["c1"]]
        out.append(dict(c, c0=c0, c1=c1))
    return out


# --------------------------------------------------------------------------
# NX Open：建样条 + 直纹面
# --------------------------------------------------------------------------
def create_spline_through_points(work_part, points):
    """用 Studio Spline（Through Points）在 NX 里建一条通过点列的样条。

    返回 NXOpen.Spline（非关联曲线对象，供直纹面使用）。
    """
    builder = work_part.Features.CreateStudioSplineBuilderEx(None)
    try:
        builder.Type = NXOpen.Features.StudioSplineBuilderEx.Types.ThroughPoints
        builder.Degree = 3
        builder.IsAssociative = False   # 只需要曲线，不需要历史特征
        for x, y, z in points:
            pt = work_part.Points.CreatePoint(NXOpen.Point3d(x, y, z))
            gc = builder.ConstraintManager.CreateGeometricConstraintData()
            gc.Point = pt
            builder.ConstraintManager.Append(gc)
        builder.Evaluate()
        spline = builder.Curve
        return spline
    finally:
        try:
            builder.Destroy()
        except Exception:
            pass


def add_curve_to_section(section, curve):
    """把一条曲线加入 NXOpen.Section（不同 NX 版本签名略有差异，做了兜底）。"""
    errors = []
    # 方式 1：最简单（部分新版本支持）
    try:
        section.AddToSection(curve)
        return
    except Exception as e:
        errors.append(str(e))
    # 方式 2：带规则 + 起终点（journal 常见形式）
    try:
        start = curve.StartPoint if hasattr(curve, "StartPoint") else None
        end = curve.EndPoint if hasattr(curve, "EndPoint") else None
        section.AddToSection(None, curve, start, end, start, end,
                             False, False, False)
        return
    except Exception as e:
        errors.append(str(e))
    # 方式 3：只带规则
    try:
        section.AddToSection(None, curve)
        return
    except Exception as e:
        errors.append(str(e))
    raise RuntimeError("AddToSection 失败: " + " | ".join(errors))


def create_ruled_surface(work_part, spline0, spline1):
    """由两条准线样条建「直纹面」特征（RuledBuilder）。"""
    builder = work_part.Features.CreateRuledBuilder(None)
    try:
        s1 = work_part.Sections.CreateSection(0.00095, 0.00095, 0.05)
        s2 = work_part.Sections.CreateSection(0.00095, 0.00095, 0.05)
        add_curve_to_section(s1, spline0)
        add_curve_to_section(s2, spline1)
        builder.FirstSection = s1
        builder.SecondSection = s2
        # 参数对齐（Parameter）：对应 CATIA 的"直纹"参数对应
        try:
            builder.AlignmentMethod.Method = \
                NXOpen.GeometricUtilities.AlignmentMethodBuilder.MethodType.Parameter
        except Exception:
            pass
        feature = builder.CommitFeature()
        return feature
    finally:
        try:
            builder.Destroy()
        except Exception:
            pass


def run_nx_export(cells, out_prt):
    session = NXOpen.Session.GetSession()
    work_part = session.Parts.Work
    if work_part is None:
        raise RuntimeError("没有活动的 Work Part，请先新建或打开一个 Part")

    ok = 0
    failed = []
    for idx, c in enumerate(cells):
        try:
            s0 = create_spline_through_points(work_part, c["c0"])
            s1 = create_spline_through_points(work_part, c["c1"])
            create_ruled_surface(work_part, s0, s1)
            ok += 1
        except Exception as e:
            failed.append((c["index"], str(e)))
        if (idx + 1) % 50 == 0 or (idx + 1) == len(cells):
            print(f"[nx] progress {idx + 1}/{len(cells)} (ok={ok}, failed={len(failed)})")

    # 保存
    if out_prt:
        base = os.path.splitext(os.path.abspath(out_prt))[0]
        save_path = out_prt
        if os.path.exists(save_path):
            save_path = f"{base}_1.prt"
        work_part.SaveAs(os.path.abspath(save_path))
        print(f"[nx] saved {save_path} ({ok}/{len(cells)} ruled surfaces, "
              f"{len(failed)} failed)")

    if failed:
        print("[nx] failed cells:")
        for i, e in failed[:50]:
            print(f"  cell[{i}]: {e}")
    return ok, failed


def main():
    import argparse
    ap = argparse.ArgumentParser(description="Export ruled fitting to NX .prt")
    ap.add_argument("json_path", help="path to <blade>_directrices.json")
    ap.add_argument("--out", default=None, help="output .prt path")
    ap.add_argument("--dry-run", action="store_true",
                    help="only validate JSON and print plan (no NX)")
    ap.add_argument("--resample", type=int, default=25,
                    help="downsample each directrix to N points (default 25, 0=off)")
    ap.add_argument("--no-prioritize", action="store_true",
                    help="keep original cell order")
    ap.add_argument("--no-translate", dest="translate", action="store_false",
                    help="keep original coordinates (default: translate to origin)")
    args = ap.parse_args()

    with open(args.json_path) as f:
        data = json.load(f)
    cells = data["cells"]
    name = data.get("name", "Blade")
    print(f"[nx] {name}: {len(cells)} ruled cells, resample={args.resample}")

    if args.resample > 0:
        cells = [dict(c, c0=resample(c["c0"], args.resample),
                      c1=resample(c["c1"], args.resample)) for c in cells]

    if args.translate:
        center = bbox_center(cells)
        cells = translate_cells(cells, center)
        print(f"[nx] translated to origin by -({center[0]:.2f}, "
              f"{center[1]:.2f}, {center[2]:.2f})")

    if not args.no_prioritize:
        cells = sorted(cells, key=directrix_length, reverse=True)

    if args.dry_run or not HAS_NXOPEN:
        for c in cells[:10]:
            print(f"  cell[{c['index']}] fitDir={c['fitDir']} "
                  f"c0={len(c['c0'])}pts c1={len(c['c1'])}pts")
        print("[nx] dry-run done (no NX connection)")
        return 0

    out = args.out or (os.path.splitext(args.json_path)[0] + ".prt")
    run_nx_export(cells, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
