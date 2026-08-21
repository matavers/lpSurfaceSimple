#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
catia_export.py — 把直纹面拟合结果导出为 CATIA V5 原生 CATPart。

读取 simple.exe 生成的 `<blade>_directrices.json`（每个直纹格的两条优化准线
采样点），通过 CATIA COM 自动化逐格重建「直纹面（HybridShapeRuledSurface）」，
保存为原生 CATPart。直纹语义被完整保留，CATIA 加工模块可直接用侧刃铣
（Multi-Axis Flank Contouring / Sweeping）策略。

用法:
    python catia_export.py <blade1_directrices.json> [--out out.CATPart] [--dry-run]

依赖:
    - 本机已安装 CATIA V5（脚本会自动启动/连接）
    - pip install pywin32

说明:
    - 直纹面 = 两条准线 + 母线。脚本为每条准线建三次样条（HybridShapeSpline），
      再用 HybridShapeRuledSurface 由两条样条生成直纹面。
    - 若你的 CATIA 版本对 AddPointWithConstraintExplicit 签名不同，可调整
      make_spline 里的调用参数。
"""

import json
import os
import sys
import argparse


def make_spline(hsf, gs, part, pts):
    """由点列在 CATIA 里建一条三次样条曲线（作为直纹面准线）。"""
    spline = hsf.AddNewSpline()
    spline.SetSplineType(0)   # 0 = 三次样条
    for x, y, z in pts:
        pt = hsf.AddNewPointCoord(x, y, z)
        gs.AppendHybridShape(pt)
        part.Update()
        ref = part.CreateReferenceFromObject(pt)
        try:
            spline.AddPointWithConstraintExplicit(ref, None, -1.0, 1, 0.0, 0.0, 0.0)
        except Exception:
            # 兼容更简签名
            spline.AddPointWithConstraintExplicit(ref)
        part.Update()
    gs.AppendHybridShape(spline)
    part.Update()
    return spline


def build_ruled(part, hsf, gs, c0_pts, c1_pts, name):
    s0 = make_spline(hsf, gs, part, c0_pts)
    s1 = make_spline(hsf, gs, part, c1_pts)
    ruled = hsf.AddNewRuledSurface(
        part.CreateReferenceFromObject(s0),
        part.CreateReferenceFromObject(s1))
    ruled.set_Name(name)
    gs.AppendHybridShape(ruled)
    part.Update()
    return ruled


def main():
    ap = argparse.ArgumentParser(description="Export ruled fitting to CATIA CATPart")
    ap.add_argument("json_path", help="path to <blade>_directrices.json")
    ap.add_argument("--out", default=None, help="output .CATPart path")
    ap.add_argument("--dry-run", action="store_true",
                    help="only validate JSON and print plan (no CATIA)")
    args = ap.parse_args()

    with open(args.json_path) as f:
        data = json.load(f)
    cells = data["cells"]
    name = data.get("name", "Blade")
    print(f"[catia] {name}: {len(cells)} ruled cells")

    if args.dry_run:
        for c in cells:
            print(f"  cell[{c['index']}] fitDir={c['fitDir']} "
                  f"c0={len(c['c0'])}pts c1={len(c['c1'])}pts")
        print("[catia] dry-run done (no CATIA connection)")
        return 0

    import win32com.client
    catia = win32com.client.Dispatch("CATIA.Application")
    catia.Visible = True

    doc = catia.Documents.Add("Part")
    part = doc.Part
    hsf = part.HybridShapeFactory

    gs = part.HybridBodies.Add()
    gs.set_Name(f"{name}_ruled")

    ok = 0
    for c in cells:
        try:
            build_ruled(part, hsf, gs, c["c0"], c["c1"], f"ruled_{c['index']}")
            ok += 1
        except Exception as e:
            print(f"[catia] cell[{c['index']}] failed: {e}")

    part.Update()
    out = args.out or (os.path.splitext(args.json_path)[0] + ".CATPart")
    doc.SaveAs(out)
    print(f"[catia] saved {out} ({ok}/{len(cells)} ruled surfaces)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
