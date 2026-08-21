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
"""

import json
import os
import sys
import argparse


def flag_methods(obj, *names):
    """晚绑定下把名字强制标记为方法（修复 Add 等被误判成属性的问题）。"""
    try:
        obj._FlagAsMethod(*names)
    except Exception:
        pass


def set_name(obj, name):
    try:
        obj.Name = name
    except Exception:
        try:
            obj.set_Name(name)
        except Exception:
            pass


def connect_catia():
    # 用 dynamic.Dispatch 强制晚绑定：避免 gencache 早绑定把 Document 类型化后
    # 取不到 Part 属性；同时 dynamic 会带上类型信息做惰性方法解析，Add() 也能正确识别。
    from win32com.client import dynamic
    catia = dynamic.Dispatch("CATIA.Application")
    catia.Visible = True
    try:
        catia.Interactive = False      # 抑制 CATIA 错误弹窗，失败直接抛异常
    except Exception:
        pass
    try:
        catia.RefreshDisplay = False   # 关显示刷新，大幅提速
    except Exception:
        pass
    return catia


def make_spline(hsf, gs, part, pts):
    """由点列在 CATIA 里建一条三次样条曲线（作为直纹面准线）。不在每个点后
    Update，最后统一 Update，避免每个点都重算整个 part（提速关键）。"""
    spline = hsf.AddNewSpline()
    try:
        spline.SetSplineType(0)   # 0 = 三次样条
    except Exception:
        spline.SplineType = 0
    for x, y, z in pts:
        pt = hsf.AddNewPointCoord(x, y, z)
        gs.AppendHybridShape(pt)
        spline.AddPoint(part.CreateReferenceFromObject(pt))
    gs.AppendHybridShape(spline)
    return spline


def build_ruled(part, hsf, gs, c0_pts, c1_pts, name):
    s0 = make_spline(hsf, gs, part, c0_pts)
    s1 = make_spline(hsf, gs, part, c1_pts)
    # 直纹面：CATIA V5-6R2020 无 AddNewRuledSurface，用 Line Sweep（直线扫掠，两极限）
    sw = hsf.AddNewSweepLine(part.CreateReferenceFromObject(s0))
    sw.SecondGuideCrv = part.CreateReferenceFromObject(s1)
    sw.Mode = 1   # 1 = Two limits（两极限，即直纹面）
    set_name(sw, name)
    gs.AppendHybridShape(sw)
    return sw


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

    catia = connect_catia()
    doc = catia.Documents.Add("Part")
    part = doc.Part
    hsf = part.HybridShapeFactory

    # 晚绑定下把易误判的方法名显式标记为方法
    flag_methods(part.HybridBodies, "Add")
    flag_methods(hsf, "AddNewSpline", "AddNewPointCoord",
                 "AddNewSweepLine", "AddNewLoft")
    flag_methods(part, "CreateReferenceFromObject")

    gs = part.HybridBodies.Add()
    flag_methods(gs, "AppendHybridShape")
    set_name(gs, f"{name}_ruled")

    ok = 0
    failed = []
    sel = doc.Selection
    for i, c in enumerate(cells):
        sw = None
        try:
            sw = build_ruled(part, hsf, gs, c["c0"], c["c1"], f"ruled_{c['index']}")
            part.Update()   # 逐格 Update，能及时发现无效扫掠并回滚
            ok += 1
        except Exception as e:
            failed.append(c["index"])
            if sw is not None:
                try:   # 回滚：删掉无效扫掠，避免污染后续 Update
                    sel.Clear()
                    sel.Add(sw)
                    sel.Delete()
                except Exception:
                    pass
            print(f"[catia] cell[{c['index']}] failed: {e}")
        if (i + 1) % 20 == 0:
            print(f"[catia] progress {i + 1}/{len(cells)} (ok={ok}, failed={len(failed)})")

    try:
        catia.RefreshDisplay = True
    except Exception:
        pass
    out = args.out or (os.path.splitext(args.json_path)[0] + ".CATPart")
    doc.SaveAs(out)
    print(f"[catia] saved {out} ({ok}/{len(cells)} ruled surfaces, "
          f"{len(failed)} failed: {failed})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
