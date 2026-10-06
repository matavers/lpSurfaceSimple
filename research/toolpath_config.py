#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""连续刀轨规划参数配置（toolpath_config.json）读写。

toolpath_continuous.py、sweep_machining.py 与 GUI 共用本模块，
统一从 toolpath_config.json 读取/写回刀轨规划参数。
"""
import json
from pathlib import Path

CONFIG_PATH = Path(__file__).resolve().parent / "toolpath_config.json"

DEFAULTS = {
    "tool_len": 25.0,        # 刀具有效切削刃长 L_tool (mm)，只是最大切削宽度上限
    "tool_r": 5.0,           # 刀具半径 R (mm)
    "o_min": 2.0,            # 最小重叠量 (mm)
    "o_max": 5.0,            # 最大重叠量 (mm)
    "eps": 0.1,              # 可接受误差（mean|e|，mm），决定切削长度
    "n_rule": 10,            # 母线方向采样数
    "target_spacing": 1.5,   # 目标进给采样间距 (mm)
    "n_jobs": -1,            # 刀轨并行进程数（-1=全部 CPU，1=串行）
}


def load_config(path=None):
    """读取配置，返回 dict（缺失键用 DEFAULTS 补齐，文件不存在则返回 DEFAULTS）。"""
    if path is None:
        path = CONFIG_PATH
    cfg = {}
    try:
        with open(path, encoding="utf-8") as f:
            cfg = json.load(f)
    except Exception:
        cfg = {}
    merged = dict(DEFAULTS)
    merged.update({k: v for k, v in cfg.items() if k in DEFAULTS})
    return merged


def save_config(cfg, path=None):
    """写回配置（仅保留 DEFAULTS 中的键，避免写入未知字段）。"""
    if path is None:
        path = CONFIG_PATH
    merged = {k: cfg.get(k, DEFAULTS[k]) for k in DEFAULTS}
    with open(path, "w", encoding="utf-8") as f:
        json.dump(merged, f, ensure_ascii=False, indent=4)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--init":
        save_config({}, CONFIG_PATH)
        print(f"已初始化 {CONFIG_PATH}")
    else:
        print(json.dumps(load_config(), ensure_ascii=False, indent=2))
