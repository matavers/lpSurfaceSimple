#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dump_split.py — 导出并 3D 可视化 split-blade 结果（pyvista 小窗口），
每个 region 用不同颜色渲染，便于检查叶盆/叶背/叶缘切分、以及卷曲部分是否被切出。

用法:
  python dump_split.py <file> [out.json]
"""
import os
import sys
import json
import tempfile
import subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import sweep_machining as sm

COLORS = {"pressure": "#1f77b4", "suction": "#2ca02c", "edge": "#d62728"}


def export_region_mesh(file, face_idx, dir_, us, ue, out_path):
    rk = "v" if dir_ == "V" else "u"
    subprocess.run(
        [str(sm.BUILD_EXE), file, "--mode", "face-obj-range",
         "--face-idx", str(face_idx), f"--{rk}-range1", f"{us},{ue}",
         "--face-out", out_path],
        cwd=str(sm.PROJECT_DIR), capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=60)
    return os.path.exists(out_path) and os.path.getsize(out_path) > 100


def main():
    file = sys.argv[1]
    out_json = sys.argv[2] if len(sys.argv) > 2 else "blade_split.json"

    pidx, sidx = sm.auto_identify(file)
    face_idx = pidx if pidx >= 0 else 0
    dir_, regions = sm.split_blade(file, face_idx)
    print(f"[dump] face_idx={face_idx} dir={dir_}")

    result = {"file": file, "face_idx": face_idx, "dir": dir_, "regions": regions}
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=1)
    print(f"[dump] regions -> {out_json}")
    for reg in regions:
        print(f"  {reg.get('label'):10s} V[{reg.get('uStart'):.4f}, {reg.get('uEnd'):.4f}] "
              f"avgCurv={reg.get('avgCurv', 0):.3f}")

    try:
        import pyvista as pv
    except Exception:
        print("[warn] pyvista 不可用，跳过 3D 可视化")
        return

    plotter = pv.Plotter(title="split-blade regions")
    tmp = tempfile.gettempdir()
    for k, reg in enumerate(regions):
        us, ue = reg["uStart"], reg["uEnd"]
        label = reg.get("label", "?")
        sub = os.path.join(tmp, f"dump_split_{face_idx}_{k}.obj")
        if export_region_mesh(file, face_idx, dir_, us, ue, sub):
            m = pv.read(sub)
            plotter.add_mesh(m, color=COLORS.get(label, "gray"),
                             name=f"{label}_{k}", show_edges=True, edge_color="black")
    plotter.add_legend([("pressure", COLORS["pressure"]),
                        ("suction", COLORS["suction"]),
                        ("edge", COLORS["edge"])])
    print("[dump] 打开 3D 窗口（可旋转视角检查卷曲部分是否被切出）")
    plotter.show()


if __name__ == "__main__":
    main()
