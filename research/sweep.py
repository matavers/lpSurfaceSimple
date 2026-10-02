#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sweep.py — 参数扫描：容差 / 点铣工艺参数 对加工时间的影响

用法:
  python sweep.py <blade_file> [--exe path/to/simple.exe] \
      [--tolerances 0.05,0.1,0.2,0.5,1.0] [--ball-r 5] [--scallop 0.01] \
      [--twist-limit 1.0] [--feed 500] [--json sweep.json]
"""
import os
import sys
import subprocess
import shutil
import json
import argparse
import tempfile
from pathlib import Path

try:
    from research import compute_machining as cm
except ImportError:
    import compute_machining as cm


def run_fitting(exe, blade, out_dir, tolerance):
    """调用 simple.exe --mode ruled（井字形网格细分）生成输出。返回错误码。"""
    os.makedirs(out_dir, exist_ok=True)
    cmd = [exe, blade, blade, '--mode', 'ruled',
           '--tolerance', str(tolerance), '--outdir', out_dir]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True)
        return r.returncode
    except Exception:
        return -1


def main():
    cfg = cm.load_config()
    ap = argparse.ArgumentParser(description="参数扫描：容差 vs 加工时间")
    ap.add_argument("blade", help="叶片 STEP/IGES 文件路径")
    ap.add_argument("--exe", default=r"D:\Projects\lpSurface\Simple-gr\build\Release\simple.exe")
    ap.add_argument("--tolerances", default=",".join(str(x) for x in cfg.get("tolerances", [0.05, 0.1, 0.2, 0.5, 1.0])),
                    help="容差列表（逗号分隔，mm）")
    ap.add_argument("--feed", type=float, default=cfg.get("feed", 500.0))
    ap.add_argument("--ball-r", type=float, default=cfg.get("ball_r", 5.0))
    ap.add_argument("--scallop", type=float, default=cfg.get("scallop", 0.01))
    ap.add_argument("--twist-limit", type=float, default=cfg.get("twist_limit", 1.0))
    ap.add_argument("--overhead", type=float, default=cfg.get("overhead", 2.0))
    ap.add_argument("--point-overhead", type=float, default=cfg.get("point_overhead", 10.0))
    ap.add_argument("--json", default="sweep.json")
    ap.add_argument("--workdir", default=None, help="临时输出根目录")
    args = ap.parse_args()

    tolerances = [float(x) for x in args.tolerances.split(',')]
    workdir = args.workdir or tempfile.mkdtemp(prefix="sweep_")

    rows = []
    for tol in tolerances:
        out_dir = os.path.join(workdir, f"tol_{tol:.3f}")
        if os.path.isdir(out_dir):
            shutil.rmtree(out_dir)
        code = run_fitting(args.exe, args.blade, out_dir, tol)
        if code != 0:
            print(f"  tolerance={tol}: 拟合失败 (code={code})")
            continue

        # 用与主脚本相同的参数
        a = argparse.Namespace(feed=args.feed, ball_r=args.ball_r, scallop=args.scallop,
                               twist_limit=args.twist_limit, overhead=args.overhead,
                               point_overhead=args.point_overhead)
        patches = cm.compute(out_dir, a)
        if not patches:
            continue
        s = cm.summarize(patches, a)
        rows.append({
            "tolerance": tol,
            "num_patches": s["num_patches"],
            "developable": s["developable"],
            "non_developable": s["non_developable"],
            "max_twist": max((p.twist for p in patches), default=0.0),
            "flank_total": s["flank"]["total"],
            "point_total": s["point"]["total"],
            "speedup": s["point"]["total"] / s["flank"]["total"] if s["flank"]["total"] > 0 else 0.0,
        })

    print("=" * 80)
    print("参数扫描结果（容差 vs 加工时间）")
    print("=" * 80)
    print(f"{'容差mm':>8}{'面片数':>8}{'可展':>6}{'不可展':>7}{'最大扭转°':>10}{'侧铣s':>9}{'点铣s':>9}{'提速':>7}")
    for r in rows:
        print(f"{r['tolerance']:>8.3f}{r['num_patches']:>8}{r['developable']:>6}{r['non_developable']:>7}"
              f"{r['max_twist']:>10.2f}{r['flank_total']:>9.1f}{r['point_total']:>9.1f}{r['speedup']:>7.1f}x")

    if args.json:
        with open(args.json, 'w', encoding='utf-8') as f:
            json.dump({"blade": args.blade, "params": {
                "feed": args.feed, "ball_r": args.ball_r, "scallop": args.scallop,
                "twist_limit": args.twist_limit, "overhead": args.overhead,
            }, "rows": rows}, f, ensure_ascii=False, indent=2)
        print(f"\n结果已写入: {args.json}")

if __name__ == "__main__":
    main()
