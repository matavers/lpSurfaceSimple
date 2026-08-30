#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compute_machining.py — 直纹面分片拟合 vs 传统点铣 加工时间仿真（自算刀轨）

读取 API 输出的直纹面准线参数文件 (*_params.txt)，
1. 计算每块的可展性指标（母线扭转角）
2. 生成侧铣刀轨（刀具轴线 ∥ 母线，沿准线进给）并计算切削/非切削时间
3. 生成点铣刀轨（球头刀沿准线方向行切）并计算切削时间
4. 输出对比表 + 可展性分类 + JSON 结果 + 可选 VTK 刀轨

用法:
  python compute_machining.py <output_dir> [--feed 500] [--rapid 5000] \
      [--ball-r 5] [--scallop 0.01] [--overhead 2] [--twist-limit 1.0] \
      [--json out.json] [--vtk vtk_dir]
"""
import os
import sys
import math
import json
import argparse
from pathlib import Path
from dataclasses import dataclass, asdict

# ===================== 几何 =====================

def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

def _norm(a):
    return math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])

def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])

def _lerp(a, b, t):
    return tuple(a[k] * (1.0 - t) + b[k] * t for k in range(3))

def load_directrix(path):
    """读取 _params.txt，返回 (C0, C1)，各为 [(x,y,z),...]"""
    C0, C1, cur = [], [], None
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('['):
                cur = line.strip('[]')
            elif cur in ('C0', 'C1') and line and not line.startswith(('n', 'd', 'm')):
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        p = tuple(float(x) for x in parts[:3])
                        (C0 if cur == 'C0' else C1).append(p)
                    except ValueError:
                        pass
    return C0, C1

def curve_length(C):
    return sum(_norm(_sub(C[i + 1], C[i])) for i in range(len(C) - 1))

def mean_ruling_length(C0, C1):
    n = min(len(C0), len(C1))
    if n == 0:
        return 0.0
    return sum(_norm(_sub(C1[i], C0[i])) for i in range(n)) / n

def twist_angle_deg(C0, C1):
    """最大母线扭转角（度）：母线方向相对平均方向的偏差，衡量可展性。"""
    n = min(len(C0), len(C1))
    avg = [0.0, 0.0, 0.0]
    cnt = 0
    for i in range(n):
        r = _sub(C1[i], C0[i])
        m = _norm(r)
        if m < 1e-12:
            continue
        avg = [avg[k] + r[k] / m for k in range(3)]
        cnt += 1
    if cnt == 0:
        return 0.0
    mavg = _norm(avg)
    if mavg < 1e-12:
        return 0.0
    avg = [avg[k] / mavg for k in range(3)]
    mx = 0.0
    for i in range(n):
        r = _sub(C1[i], C0[i])
        m = _norm(r)
        if m < 1e-12:
            continue
        d = max(-1.0, min(1.0, _dot(avg, [r[k] / m for k in range(3)])))
        mx = max(mx, math.degrees(math.acos(d)))
    return mx

def surface_area(C0, C1):
    """直纹面面积 = Σ |r × dC|（数值积分）。"""
    n = min(len(C0), len(C1))
    area = 0.0
    for i in range(n - 1):
        r = _sub(C1[i], C0[i])
        dC = _sub(C0[i + 1], C0[i])
        area += _norm(_cross(r, dC))
    return area

# ===================== 刀轨 + 时间 =====================

def flank_milling(C0, C1, feed, twist_limit):
    """
    侧铣：刀具轴线 ∥ 母线，沿准线进给。
    可展（扭转角 <= 阈值）：单刀成形，刀轨长 = 准线长。
    不可展：按扭转角/阈值近似放大刀数（锯齿多刀）。
    返回 (切削时间 s, 刀轨长 mm, 刀数, 是否可展, 扭转角 deg)
    """
    L = curve_length(C0)
    twist = twist_angle_deg(C0, C1)
    developable = twist <= twist_limit
    passes = 1.0 if developable else max(1.0, twist / twist_limit)
    path_len = L * passes
    return path_len / feed * 60.0, path_len, passes, developable, twist

def point_milling(C0, C1, feed, ball_r, scallop):
    """
    点铣：球头刀沿准线方向行切，行距由残留高度反算。
    返回 (切削时间 s, 行距 mm, 总刀轨长 mm, 面积 mm²)
    """
    stepover = 2.0 * math.sqrt(max(0.0, 2.0 * ball_r * scallop - scallop * scallop))
    A = surface_area(C0, C1)
    total_len = A / stepover if stepover > 0 else 0.0
    return total_len / feed * 60.0, stepover, total_len, A

def flank_toolpath_lines(C0, C1, n_along=50):
    """侧铣刀轨：沿准线的刀具中心路径 + 每条母线的轴线段（可视化）。"""
    n = len(C0)
    if n < 2:
        return []
    idx = ([int(round(i * (n - 1) / (n_along - 1))) for i in range(n_along)]
           if n > n_along else list(range(n)))
    lines = [[C0[i] for i in idx]]  # 进给路径（下准线）
    for i in idx:
        lines.append([C0[i], C1[i]])  # 母线（刀具轴线方向）
    return lines

def point_toolpath_lines(C0, C1, stepover, n_along=50):
    """点铣刀轨：沿准线方向行切，行距 stepover（沿母线方向）。"""
    n = len(C0)
    if n < 2:
        return []
    idx = ([int(round(i * (n - 1) / (n_along - 1))) for i in range(n_along)]
           if n > n_along else list(range(n)))
    mean_r = mean_ruling_length(C0, C1)
    n_across = max(1, int(round(mean_r / stepover))) if stepover > 0 else 1
    lines = []
    for j in range(n_across + 1):
        t = j / n_across if n_across > 0 else 0.0
        lines.append([_lerp(C0[i], C1[i], t) for i in idx])
    return lines

def write_vtk_polylines(path, lines):
    pts, segs = [], []
    for line in lines:
        base = len(pts)
        for p in line:
            pts.append(p)
        for k in range(len(line) - 1):
            segs.append((base + k, base + k + 1))
    with open(path, 'w', encoding='utf-8') as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("toolpath\nASCII\nDATASET POLYDATA\n")
        f.write(f"POINTS {len(pts)} float\n")
        for p in pts:
            f.write(f"{p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
        f.write(f"LINES {len(segs)} {len(segs) * 3}\n")
        for a, b in segs:
            f.write(f"2 {a} {b}\n")

# ===================== 数据模型 =====================

@dataclass
class Patch:
    name: str
    C0: list = None
    C1: list = None
    directrix_len: float = 0.0
    mean_ruling: float = 0.0
    area: float = 0.0
    twist: float = 0.0
    developable: bool = True
    flank_time: float = 0.0
    flank_len: float = 0.0
    flank_passes: float = 0.0
    point_time: float = 0.0

# ===================== 标定参数 =====================

def load_config(path=None):
    """读取标定参数 JSON 配置，返回 dict（不存在则空）。"""
    if path is None:
        path = Path(__file__).parent / "machining_config.json"
    try:
        with open(path, encoding='utf-8') as f:
            return json.load(f)
    except Exception:
        return {}

# ===================== 主流程 =====================

def compute(input_dir, args):
    """扫描 _params.txt，计算每块指标，返回 patches 列表。"""
    files = sorted(Path(input_dir).glob("*_params.txt"))
    if not files:
        print(f"[Error] 未找到 *_params.txt 文件: {input_dir}")
        return []
    patches = []
    for fp in files:
        C0, C1 = load_directrix(fp)
        if len(C0) < 2 or len(C1) < 2:
            print(f"  [skip] {fp.name}: 准线点数不足")
            continue
        p = Patch(name=fp.name.replace('_params.txt', ''), C0=C0, C1=C1)
        p.directrix_len = curve_length(C0)
        p.mean_ruling = mean_ruling_length(C0, C1)
        p.area = surface_area(C0, C1)
        p.flank_time, p.flank_len, p.flank_passes, p.developable, p.twist = \
            flank_milling(C0, C1, args.feed, args.twist_limit)
        p.point_time, _, _, _ = point_milling(C0, C1, args.feed, args.ball_r, args.scallop)
        patches.append(p)
    return patches

def summarize(patches, args):
    """汇总并返回结果 dict。"""
    n_dev = sum(1 for p in patches if p.developable)
    n_nondev = len(patches) - n_dev
    flank_cut = sum(p.flank_time for p in patches)
    flank_overhead = len(patches) * args.overhead
    flank_total = flank_cut + flank_overhead
    point_cut = sum(p.point_time for p in patches)
    point_overhead = getattr(args, 'point_overhead', 10.0)
    point_total = point_cut + point_overhead
    total_area = sum(p.area for p in patches)
    return {
        "num_patches": len(patches), "developable": n_dev, "non_developable": n_nondev,
        "total_area": total_area,
        "flank": {"cut": flank_cut, "overhead": flank_overhead, "total": flank_total},
        "point": {"cut": point_cut, "overhead": point_overhead, "total": point_total},
    }

def print_report(patches, args, summary):
    print("=" * 92)
    print(f"直纹面拟合 vs 传统点铣 加工时间仿真  (feed={args.feed} mm/min, 球刀R={args.ball_r}, 残留={args.scallop})")
    print(f"面片数: {summary['num_patches']}")
    print("=" * 92)
    print(f"{'面片':<30}{'可展':<5}{'扭转°':>8}{'准线mm':>9}{'母线mm':>9}{'侧铣s':>9}{'点铣s':>9}")
    for p in patches:
        print(f"{p.name:<30}{'是' if p.developable else '否':<5}{p.twist:>8.2f}"
              f"{p.directrix_len:>9.1f}{p.mean_ruling:>9.1f}"
              f"{p.flank_time:>9.2f}{p.point_time:>9.2f}")
    print("-" * 92)
    print(f"可展面片: {summary['developable']}  不可展面片: {summary['non_developable']}  总曲面面积: {summary['total_area']:.1f} mm²")
    print()
    f, p = summary["flank"], summary["point"]
    print(f"侧铣（直纹面拟合后）: 切削 {f['cut']:8.1f}s + 非切削 {f['overhead']:8.1f}s = {f['total']:8.1f}s")
    print(f"点铣（传统整体）:     切削 {p['cut']:8.1f}s + 非切削 {p['overhead']:8.1f}s = {p['total']:8.1f}s")
    if f['total'] > 0:
        print(f"侧铣相对点铣提速: {p['total'] / f['total']:.1f}x")

def main():
    cfg = load_config()
    ap = argparse.ArgumentParser(description="直纹面拟合加工时间仿真（自算刀轨）")
    ap.add_argument("input", help="输出目录（含 *_params.txt）")
    ap.add_argument("--feed", type=float, default=cfg.get("feed", 500.0), help="切削进给 mm/min")
    ap.add_argument("--rapid", type=float, default=cfg.get("rapid", 5000.0), help="快进 mm/min")
    ap.add_argument("--ball-r", type=float, default=cfg.get("ball_r", 5.0), help="球头刀半径 mm")
    ap.add_argument("--scallop", type=float, default=cfg.get("scallop", 0.01), help="点铣残留高度 mm")
    ap.add_argument("--overhead", type=float, default=cfg.get("overhead", 2.0), help="每块侧铣进退刀+衔接时间 s")
    ap.add_argument("--point-overhead", type=float, default=cfg.get("point_overhead", 10.0), help="点铣整体进退刀时间 s")
    ap.add_argument("--twist-limit", type=float, default=cfg.get("twist_limit", 1.0), help="可展判定阈值 度")
    ap.add_argument("--json", help="可选：输出 JSON 结果文件路径")
    ap.add_argument("--vtk", help="可选：输出刀轨 VTK 到指定目录")
    args = ap.parse_args()

    patches = compute(args.input, args)
    if not patches:
        sys.exit(1)
    summary = summarize(patches, args)
    print_report(patches, args, summary)

    if args.vtk:
        os.makedirs(args.vtk, exist_ok=True)
        for p in patches:
            base = os.path.join(args.vtk, p.name)
            write_vtk_polylines(base + "_flank.vtk", flank_toolpath_lines(p.C0, p.C1))
            stepover = 2.0 * math.sqrt(max(0.0, 2.0 * args.ball_r * args.scallop - args.scallop * args.scallop))
            write_vtk_polylines(base + "_point.vtk", point_toolpath_lines(p.C0, p.C1, stepover))
        print(f"刀轨 VTK 已写入: {args.vtk}")

    if args.json:
        result = {
            "feed": args.feed, "ball_r": args.ball_r, "scallop": args.scallop,
            "twist_limit": args.twist_limit, "overhead": args.overhead,
            **summary,
            "patches": [{k: v for k, v in asdict(p).items() if k not in ('C0', 'C1')} for p in patches],
        }
        with open(args.json, 'w', encoding='utf-8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
        print(f"结果已写入: {args.json}")

if __name__ == "__main__":
    main()
