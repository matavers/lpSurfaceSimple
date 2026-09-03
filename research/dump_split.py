#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dump_split.py — 导出并可视化 split-blade 结果，便于检查叶盆/叶背/叶缘切分是否正确。

用法:
  python dump_split.py <file> [out.json] [out.png]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import sweep_machining as sm


def main():
    file = sys.argv[1]
    out_json = sys.argv[2] if len(sys.argv) > 2 else "blade_split.json"
    out_png = sys.argv[3] if len(sys.argv) > 3 else "blade_split.png"

    pidx, sidx = sm.auto_identify(file)
    face_idx = pidx if pidx >= 0 else 0
    dir_, regions = sm.split_blade(file, face_idx)

    result = {"file": file, "face_idx": face_idx, "dir": dir_, "regions": regions}
    with open(out_json, "w", encoding="utf-8") as f:
        json = __import__("json")
        json.dump(result, f, ensure_ascii=False, indent=1)
    print(f"[dump] regions -> {out_json}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = {"pressure": "#1f77b4", "suction": "#2ca02c", "edge": "#d62728"}
    fig, ax = plt.subplots(figsize=(10, 2.5), dpi=120)
    for reg in regions:
        us = reg["uStart"]
        ue = reg["uEnd"]
        if ue < us:
            ue += 1.0
        label = reg.get("label", "?")
        ax.axvspan(us, ue, color=colors.get(label, "gray"), alpha=0.6)
        ax.text((us + ue) / 2, 0.5, label, ha="center", va="center", fontsize=10)
    ax.set_xlim(0, 1.0)
    ax.set_yticks([])
    ax.set_xlabel("V")
    ax.set_title(f"split-blade face {face_idx} ({dir_})")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)
    print(f"[dump] visualization -> {out_png}")


if __name__ == "__main__":
    main()
