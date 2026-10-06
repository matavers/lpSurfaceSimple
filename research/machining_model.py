# -*- coding: utf-8 -*-
"""加工参数与时间估计（工程常用公式）。

切削速度 Vc (m/min) → 主轴转速 n (rpm)   = 1000·Vc / (π·D)
每齿进给 fz (mm/齿) → 进给率 Vf (mm/min) = fz·z·n
材料去除率 MRR (cm³/min)                 = a_p·a_e·Vf / 1000
主切削力   Fc (N)                         = kc·a_p·a_e
切削功率   Pc (kW)                        = Fc·Vc / 60000
主轴扭矩   M  (N·m)                       = Fc·D / 2000
切削时间   t  (min)                       = 刀轨长(mm) / Vf(mm/min)
"""
import json
import math
from pathlib import Path

CONFIG_PATH = Path(__file__).resolve().parent / "machining_config.json"


def load_config(path=None):
    """读取加工参数配置 machining_config.json，返回 dict（不存在则空）。"""
    if path is None:
        path = CONFIG_PATH
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


def spindle_speed(Vc, D):
    """主轴转速 n (rpm)。Vc: 切削速度 m/min；D: 刀具直径 mm。"""
    if D <= 0:
        return 0.0
    return 1000.0 * Vc / (math.pi * D)


def feed_rate(fz, z, n):
    """进给率 Vf (mm/min) = 每齿进给 fz × 齿数 z × 转速 n。"""
    return fz * z * n


def material_removal_rate(a_p, a_e, Vf):
    """材料去除率 MRR (cm³/min)。a_p: 轴向切深 mm；a_e: 径向切深 mm；Vf: mm/min。"""
    return a_p * a_e * Vf / 1000.0


def cutting_force(kc, a_p, a_e):
    """主切削力 Fc (N) = 比切削力 kc (N/mm²) × 轴向切深 a_p × 径向切深 a_e。"""
    return kc * a_p * a_e


def cutting_power(Fc, Vc):
    """切削功率 Pc (kW)。Fc: 切削力 N；Vc: 切削速度 m/min。"""
    return Fc * Vc / 60000.0


def torque(Fc, D):
    """主轴扭矩 M (N·m)。Fc: 切削力 N；D: 刀具直径 mm。"""
    return Fc * D / 2000.0


def estimate_point(surface_area_mm2, cfg, ball_r):
    """点铣（球头刀）时间估计。

    surface_area_mm2: 曲面面积 mm²。
    cfg: machining_config.json 内容。
    ball_r: 球头刀半径 mm。

    行距由残留高度反算 stepover = 2√(2Rh−h²)，刀轨总长 = 面积/行距，时间 = 总长/进给率。
    返回 dict。
    """
    D = 2.0 * ball_r
    Vc = float(cfg.get("cutting_speed", 120.0))
    fz = float(cfg.get("feed_per_tooth", 0.05))
    z = int(cfg.get("num_teeth", 4))
    scallop = float(cfg.get("scallop", 0.1))

    n = spindle_speed(Vc, D)
    Vf = feed_rate(fz, z, n)
    stepover = 2.0 * math.sqrt(max(0.0, 2.0 * ball_r * scallop - scallop * scallop))
    total_len = surface_area_mm2 / stepover if stepover > 0 else 0.0
    t_min = total_len / Vf if Vf > 0 else 0.0

    return {
        "diameter_mm": D,
        "cutting_speed_m_min": Vc,
        "spindle_rpm": n,
        "feed_rate_mm_min": Vf,
        "stepover_mm": stepover,
        "scallop_mm": scallop,
        "total_path_mm": total_len,
        "surface_area_mm2": surface_area_mm2,
        "cut_time_min": t_min,
        "cut_time_s": t_min * 60.0,
    }


def estimate_flank(path_length_mm, cutting_lengths_mm, cfg, tool_r):
    """侧铣加工时间/载荷估计。

    path_length_mm: 刀轨总长 (mm)。
    cutting_lengths_mm: 各条带切削长度列表 (mm)，即轴向切深 a_p。
    cfg: machining_config.json 内容（含 cutting_speed/feed_per_tooth/num_teeth/
         specific_cutting_force/radial_depth）。
    tool_r: 刀具半径 (mm)。

    返回 dict：直径/转速/进给率/切深/时间/MRR/力/功率/扭矩。
    """
    D = 2.0 * tool_r
    Vc = float(cfg.get("cutting_speed", 120.0))
    fz = float(cfg.get("feed_per_tooth", 0.05))
    z = int(cfg.get("num_teeth", 4))
    kc = float(cfg.get("specific_cutting_force", 2500.0))
    a_e = float(cfg.get("radial_depth", 0.5))

    n = spindle_speed(Vc, D)
    Vf = feed_rate(fz, z, n)
    # 轴向切深取平均切削长度（各条带沿母线方向的切削宽度）
    a_p = float(sum(cutting_lengths_mm) / len(cutting_lengths_mm)) if cutting_lengths_mm else 0.0

    t_min = path_length_mm / Vf if Vf > 0 else 0.0
    mrr = material_removal_rate(a_p, a_e, Vf)
    Fc = cutting_force(kc, a_p, a_e)
    Pc = cutting_power(Fc, Vc)
    M = torque(Fc, D)

    return {
        "diameter_mm": D,
        "cutting_speed_m_min": Vc,
        "spindle_rpm": n,
        "feed_per_tooth_mm": fz,
        "num_teeth": z,
        "feed_rate_mm_min": Vf,
        "axial_depth_mm": a_p,
        "radial_depth_mm": a_e,
        "cut_time_min": t_min,
        "cut_time_s": t_min * 60.0,
        "mrr_cm3_min": mrr,
        "cutting_force_N": Fc,
        "cutting_power_kW": Pc,
        "torque_Nm": M,
    }
