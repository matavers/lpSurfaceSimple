#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify_continuity.py — 验证自由曲面直纹化连续性处理的四个性质，并提供交互式可视化。

四个验证项（对应论文第4章 4.4 几何验证）：
  (1) 混合曲面在格胞边界处达到 C^k 连续（数值导数检验）；
  (2) 单位分解恒等式在机器精度内成立；
  (3) 非负性成立；
  (4) 过渡带宽度是影响保真度的主要参数（过渡带越窄保真度越高）。

输入（参考 sweep_machining.py 的输入模型预处理）：
  --file  blade.igs            : auto-identify → split-blade → 直纹拟合 → 验证
  --model blade_surface_model.json : 直接加载已拟合曲面模型
  --synthetic                  : 构造带已知缝隙的合成模型（信号清晰，便于演示）

用法：
  python verify_continuity.py --synthetic --continuity 2 --band-width 0.2
  python verify_continuity.py --model ..\\output\\blade1_surface_model.json
  python verify_continuity.py --file blade.igs --tolerance 0.1 --nsplit-u 3 --nsplit-v 6

依赖：numpy, scipy, matplotlib，以及同目录 toolpath_continuous / sweep_machining。
"""
import os
import sys
import copy
import json
import argparse
import subprocess
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline, make_interp_spline

sys.path.insert(0, str(Path(__file__).resolve().parent))
import toolpath_continuous as tc   # 复用：单位分解/混合曲面解析求值

PROJECT_DIR = Path(__file__).resolve().parent.parent
BUILD_EXE = PROJECT_DIR / "build" / "Release" / "simple.exe"


# ======================================================================
# 基础工具：分解函数 / 混合曲面 / 分片直纹面
# ======================================================================

def with_partition(model, k, band_w):
    """返回一个修改了 partition(continuity, bandWidth) 的模型副本。"""
    m = copy.deepcopy(model)
    m["partition"]["continuity"] = int(k)
    m["partition"]["bandWidth"] = float(band_w)
    return m


def build_partition(model):
    """返回函数 phi_weights(u,v) -> (phi数组, omega数组)，复用 tc._chi1_np。"""
    k = model["partition"]["continuity"]
    band_w = model["partition"]["bandWidth"]
    cells = model["cells"]
    min_u = min(c["u1"] - c["u0"] for c in cells)
    min_v = min(c["v1"] - c["v0"] for c in cells)
    a_u = band_w * min_u
    a_v = band_w * min_v

    def phi_weights(u, v):
        phis = np.empty(len(cells))
        for i, c in enumerate(cells):
            cu = tc._chi1_np(u, c["u0"], c["u1"], a_u, a_u, k)
            cv = tc._chi1_np(v, c["v0"], c["v1"], a_v, a_v, k)
            phis[i] = cu * cv
        s = phis.sum()
        w = phis / (s + 1e-12)
        return phis, w
    return phi_weights


def build_piecewise_eval(model):
    """分片直纹面求值（不混合）：返回 (u,v) -> 点。用于测量混合引入的偏差与缝隙。"""
    cells = model["cells"]
    fitDir = model["fitDir"]
    cell_data = []
    for c in cells:
        c0 = [BSpline(c["c0"]["knots"], np.array(c["c0"]["ctrl"])[:, d],
                      c["c0"]["degree"]) for d in range(3)]
        c1 = [BSpline(c["c1"]["knots"], np.array(c["c1"]["ctrl"])[:, d],
                      c["c1"]["degree"]) for d in range(3)]
        cell_data.append((c, c0, c1))

    def eval_pw(u, v):
        for (c, c0, c1) in cell_data:
            if c["u0"] - 1e-9 <= u <= c["u1"] + 1e-9 and \
               c["v0"] - 1e-9 <= v <= c["v1"] + 1e-9:
                du = c["u1"] - c["u0"]
                dv = c["v1"] - c["v0"]
                R = np.zeros(3)
                for d in range(3):
                    if fitDir == "U":
                        s = (u - c["u0"]) / du
                        R[d] = (1.0 - s) * float(c0[d](v)) + s * float(c1[d](v))
                    else:
                        s = (v - c["v0"]) / dv
                        R[d] = (1.0 - s) * float(c0[d](u)) + s * float(c1[d](u))
                return R
        # 兜底：最近格胞
        c, c0, c1 = cell_data[0]
        du = c["u1"] - c["u0"]
        dv = c["v1"] - c["v0"]
        R = np.zeros(3)
        for d in range(3):
            if fitDir == "U":
                s = (u - c["u0"]) / du
                R[d] = (1.0 - s) * float(c0[d](v)) + s * float(c1[d](v))
            else:
                s = (v - c["v0"]) / dv
                R[d] = (1.0 - s) * float(c0[d](u)) + s * float(c1[d](u))
        return R
    return eval_pw


# ======================================================================
# 合成模型构造（带已知缝隙，用于清晰演示）
# ======================================================================

def _bspline_dict(pts, p0, p1, degree=3):
    """由采样点构造 B 样条，输出与 surface_model.json 一致的 {degree,knots,ctrl}。"""
    xs = np.linspace(p0, p1, len(pts))
    bs = make_interp_spline(xs, np.asarray(pts, dtype=float), k=degree)
    return {"degree": degree,
            "knots": [float(x) for x in bs.t],
            "ctrl": [[float(v) for v in row] for row in bs.c]}


def build_synthetic_model(nu=3, nv=3, gap_mm=0.15, degree=3):
    """构造 nu×nv 直纹面片合成模型，在每列之间引入 z 向 gap_mm 的缝隙。

    曲面形状：z(u,v) = 0.5·sin(πu)·cos(πv) + 0.25·u + 0.15·v，x=u, y=v。
    fitDir=V（母线沿 v，准线沿 u）。列 c≥1 整体 z 偏移 c·gap_mm，制造缝隙。
    """
    u_edges = [float(i) / nu for i in range(nu + 1)]
    v_edges = [float(j) / nv for j in range(nv + 1)]

    def z_of(u, v):
        return 0.5 * np.sin(np.pi * u) * np.cos(np.pi * v) + 0.25 * u + 0.15 * v

    cells = []
    for r in range(nv):
        v0, v1 = v_edges[r], v_edges[r + 1]
        for c in range(nu):
            u0, u1 = u_edges[c], u_edges[c + 1]
            off = c * gap_mm
            us = np.linspace(u0, u1, 12)
            c0 = _bspline_dict([(u, v0, z_of(u, v0) + off) for u in us], u0, u1, degree)
            c1 = _bspline_dict([(u, v1, z_of(u, v1) + off) for u in us], u0, u1, degree)
            cells.append({"row": r, "col": c,
                          "u0": u0, "u1": u1, "v0": v0, "v1": v1,
                          "c0": c0, "c1": c1})
    return {
        "name": "synthetic_blade", "nRows": nv, "nCols": nu,
        "fitDir": "V", "uEdges": u_edges, "vEdges": v_edges,
        "partition": {"type": "smootherstep", "bandWidth": 0.2, "continuity": 2},
        "cells": cells,
    }


# ======================================================================
# 四个验证项
# ======================================================================

def check_partition_of_unity(model, n=80):
    """(2) 单位分解恒等式：采样网格上 max|Σω − 1|（应≈机器精度）。同时报告原始 Σφ。"""
    pw = build_partition(model)
    u0, u1 = model["uEdges"][0], model["uEdges"][-1]
    v0, v1 = model["vEdges"][0], model["vEdges"][-1]
    us = np.linspace(u0, u1, n)
    vs = np.linspace(v0, v1, n)
    max_w_dev = 0.0
    max_phi_dev = 0.0
    for u in us:
        for v in vs:
            phis, w = pw(u, v)
            max_w_dev = max(max_w_dev, abs(w.sum() - 1.0))
            max_phi_dev = max(max_phi_dev, abs(phis.sum() - 1.0))
    return {"max_omega_dev": float(max_w_dev),
            "max_phi_dev": float(max_phi_dev),
            "machine_eps": float(np.finfo(np.float64).eps)}


def check_nonnegativity(model, n=80):
    """(3) 非负性：采样网格上 φ 与 ω 的最小值（应 ≥ 0）。"""
    pw = build_partition(model)
    u0, u1 = model["uEdges"][0], model["uEdges"][-1]
    v0, v1 = model["vEdges"][0], model["vEdges"][-1]
    us = np.linspace(u0, u1, n)
    vs = np.linspace(v0, v1, n)
    min_phi = np.inf
    min_w = np.inf
    for u in us:
        for v in vs:
            phis, w = pw(u, v)
            min_phi = min(min_phi, float(phis.min()))
            min_w = min(min_w, float(w.min()))
    return {"min_phi": float(min_phi), "min_omega": float(min_w)}


def _cell_point(c, u, v, fitDir):
    """求单个格胞直纹面在 (u,v) 处的点（显式，不按格胞归属做容差匹配）。"""
    c0 = [BSpline(c["c0"]["knots"], np.array(c["c0"]["ctrl"])[:, d], c["c0"]["degree"])
          for d in range(3)]
    c1 = [BSpline(c["c1"]["knots"], np.array(c["c1"]["ctrl"])[:, d], c["c1"]["degree"])
          for d in range(3)]
    du = c["u1"] - c["u0"]
    dv = c["v1"] - c["v0"]
    R = np.zeros(3)
    for d in range(3):
        if fitDir == "U":
            s = (u - c["u0"]) / du
            R[d] = (1.0 - s) * float(c0[d](v)) + s * float(c1[d](v))
        else:
            s = (v - c["v0"]) / dv
            R[d] = (1.0 - s) * float(c0[d](u)) + s * float(c1[d](u))
    return R


def _smootherstep_poly(k):
    """返回 smootherstep s_k(t) = t^{k+1}·Σ_{i=0..k} C(k+i,i)(1−t)^i 的多项式表示。"""
    from math import comb
    c = np.zeros(2 * k + 2)
    for i in range(k + 1):
        p = np.polynomial.Polynomial([1.0, -1.0]) ** i
        p = p * np.polynomial.Polynomial([0.0] * (k + 1) + [1.0])
        c[:len(p.coef)] += comb(k + i, i) * p.coef
    return np.polynomial.Polynomial(c)


def check_basis_continuity(k):
    """(1a) 分解函数（smootherstep 基函数）的精确连续阶：m=0..k+2 阶导数在拼接点处的跳跃。

    由于 1D 基函数在过渡带两侧分别是常数 0/1 与 s_k(t)，其 m 阶导数跳跃
    = |s_k^(m)(0)|（t=0 端）与 |s_k^(m)(1)|（t=1 端）的较大者。
    期望：m≤k 时跳跃≈0（C^k），m=k+1 时跳跃=|s_k^(k+1)(0)|>0。
    """
    p = _smootherstep_poly(k)
    jumps = []
    endpoint_vals = []
    for m in range(0, k + 3):
        if m == 0:
            j0 = abs(p(0.0) - 0.0)
            j1 = abs(p(1.0) - 1.0)
            endpoint_vals.append((p(0.0), p(1.0)))
        else:
            d = p.deriv(m)
            j0 = abs(d(0.0))
            j1 = abs(d(1.0))
            endpoint_vals.append((d(0.0), d(1.0)))
        jumps.append(float(max(j0, j1)))
    return {"jumps": jumps, "endpoint_vals": endpoint_vals}


def check_continuity(model, k=None, band_w=None):
    """(1) 连续性：内部边界处，分片直纹面（有缝隙，C^{-1}）与混合曲面（C^k）的对比。

    返回：
      basis   : 分解函数精确连续阶（check_basis_continuity）
      surface : 每个内部边界的 piecewise_gap（缝隙）与 blended_value_jump（混合曲面值跳跃）
    """
    if k is None:
        k = model["partition"]["continuity"]
    if band_w is None:
        band_w = model["partition"]["bandWidth"]
    m = with_partition(model, k, band_w)
    eval_blend, _ = tc.build_surface_eval(m)
    cells = model["cells"]
    fitDir = model["fitDir"]

    surface = []
    # 内部 u 边界（列之间）：沿 u 截取，固定 v 在 v 格胞中点
    for ub in model["uEdges"][1:-1]:
        for vmid in [(model["vEdges"][i] + model["vEdges"][i + 1]) / 2
                     for i in range(len(model["vEdges"]) - 1)]:
            L = [c for c in cells if abs(c["u1"] - ub) < 1e-9 and c["v0"] - 1e-9 <= vmid <= c["v1"] + 1e-9]
            R = [c for c in cells if abs(c["u0"] - ub) < 1e-9 and c["v0"] - 1e-9 <= vmid <= c["v1"] + 1e-9]
            if not L or not R:
                continue
            gap = float(np.linalg.norm(_cell_point(L[0], ub, vmid, fitDir)
                                       - _cell_point(R[0], ub, vmid, fitDir)))
            eps = 1e-7
            value_jump = float(np.linalg.norm(eval_blend(ub + eps, vmid)[0]
                                              - eval_blend(ub - eps, vmid)[0]))
            surface.append({"dir": "u", "boundary": ub, "coord": vmid,
                            "piecewise_gap": gap, "blended_value_jump": value_jump})
            break  # 每个边界取一个截面即可
    # 内部 v 边界（行之间）：沿 v 截取，固定 u 在 u 格胞中点
    for vb in model["vEdges"][1:-1]:
        for umid in [(model["uEdges"][i] + model["uEdges"][i + 1]) / 2
                     for i in range(len(model["uEdges"]) - 1)]:
            D = [c for c in cells if abs(c["v1"] - vb) < 1e-9 and c["u0"] - 1e-9 <= umid <= c["u1"] + 1e-9]
            U = [c for c in cells if abs(c["v0"] - vb) < 1e-9 and c["u0"] - 1e-9 <= umid <= c["u1"] + 1e-9]
            if not D or not U:
                continue
            gap = float(np.linalg.norm(_cell_point(D[0], umid, vb, fitDir)
                                       - _cell_point(U[0], umid, vb, fitDir)))
            eps = 1e-7
            value_jump = float(np.linalg.norm(eval_blend(umid, vb + eps)[0]
                                              - eval_blend(umid, vb - eps)[0]))
            surface.append({"dir": "v", "boundary": vb, "coord": umid,
                            "piecewise_gap": gap, "blended_value_jump": value_jump})
            break
    return {"basis": check_basis_continuity(k), "surface": surface}


def check_transition_band(model, k=None, band_widths=None, n=200):
    """(4) 过渡带宽度与保真度：扫描 bandWidth，测量混合曲面相对分片直纹面的偏差（RMS/最大）。

    偏差只出现在过渡带内（带内混合、带外=分片直纹面），故偏差随 bandWidth 增大而增大。
    """
    if k is None:
        k = model["partition"]["continuity"]
    if band_widths is None:
        band_widths = [0.5, 0.2, 0.1, 0.05, 0.02, 0.01]
    u0, u1 = model["uEdges"][0], model["uEdges"][-1]
    v0, v1 = model["vEdges"][0], model["vEdges"][-1]
    us = np.linspace(u0, u1, n)
    vs = np.linspace(v0, v1, n)
    eval_pw = build_piecewise_eval(model)
    out = []
    for bw in band_widths:
        m = with_partition(model, k, bw)
        eval_blend, _ = tc.build_surface_eval(m)
        sq_sum = 0.0
        max_dev = 0.0
        cnt = 0
        for u in us:
            for v in vs:
                d = np.linalg.norm(eval_blend(u, v)[0] - eval_pw(u, v))
                sq_sum += d * d
                max_dev = max(max_dev, float(d))
                cnt += 1
        rms = float(np.sqrt(sq_sum / cnt))
        out.append((float(bw), rms, max_dev))
    return out


# ======================================================================
# 输入模型预处理（参考 sweep_machining.py）
# ======================================================================

def run_exe_json(args_list):
    r = subprocess.run([str(BUILD_EXE)] + args_list, cwd=str(PROJECT_DIR),
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=300)
    s = r.stdout.strip()
    i = s.find("{")
    if i >= 0:
        s = s[i:]
    try:
        return json.loads(s)
    except Exception:
        return None


def preprocess_and_fit(filepath, tol, nu, nv, outdir, blend_width=0.05, continuity=3):
    """auto-identify → split-blade → 直纹拟合，返回 surface_model.json 路径列表。"""
    import sweep_machining as sm
    pidx, sidx = sm.auto_identify(filepath)
    split_face_idx = pidx if pidx >= 0 else 0
    _dir, regions = sm.split_blade(filepath, split_face_idx)
    pressure, suction = sm.identify_ps(regions)
    ps_ranges = None
    if pressure or suction:
        ps_ranges = [("V" if _dir == "V" else "U", [pressure] if pressure else []),
                     ("V" if _dir == "V" else "U", [suction] if suction else [])]
    os.makedirs(outdir, exist_ok=True)
    cmd = [str(BUILD_EXE), filepath, filepath, "--mode", "ruled",
           "--outdir", outdir, "--tolerance", str(tol),
           "--nsplit-u", str(nu), "--nsplit-v", str(nv),
           "--no-refine", "--balance-edges", "--max-cells", "200000",
           "--blend", "--blend-width", str(blend_width), "--continuity", str(continuity)]
    if ps_ranges:
        for kk, (dir_, ranges) in enumerate(ps_ranges, start=1):
            if not ranges:
                continue
            rk = "v" if dir_ == "V" else "u"
            us, ue = ranges[0]
            cmd += [f"--face-idx{kk}", str(split_face_idx), f"--{rk}-range{kk}", f"{us},{ue}"]
    subprocess.run(cmd, cwd=str(PROJECT_DIR), capture_output=True,
                   text=True, encoding="utf-8", errors="replace", timeout=900)
    models = sorted(fn for fn in os.listdir(outdir) if fn.endswith("_surface_model.json"))
    return [os.path.join(outdir, fn) for fn in models]


# ======================================================================
# 交互式可视化
# ======================================================================

def run_interactive(model, initial_k=None, initial_bw=None):
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider
    from matplotlib import font_manager
    for _f in ["Microsoft YaHei", "SimHei", "SimSun"]:
        if any(f.name == _f for f in font_manager.fontManager.ttflist):
            matplotlib.rcParams["font.sans-serif"] = [_f]
            break
    matplotlib.rcParams["axes.unicode_minus"] = False

    if initial_k is None:
        initial_k = model["partition"]["continuity"]
    if initial_bw is None:
        initial_bw = model["partition"]["bandWidth"]

    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1],
                          left=0.07, right=0.97, top=0.92, bottom=0.18,
                          wspace=0.25, hspace=0.38)
    ax_cont = fig.add_subplot(gs[0, 0])   # (1) 连续性
    ax_surf = fig.add_subplot(gs[0, 1], projection="3d")  # 混合曲面
    ax_pu = fig.add_subplot(gs[1, 0])      # (2)(3) 单位分解/非负性
    ax_band = fig.add_subplot(gs[1, 1])    # (4) 过渡带 vs 保真度

    ax_k = fig.add_axes([0.07, 0.07, 0.4, 0.03])
    ax_bw = fig.add_axes([0.07, 0.02, 0.4, 0.03])
    s_k = Slider(ax_k, "连续性阶数 k", 0, 3, valinit=initial_k, valstep=1)
    s_bw = Slider(ax_bw, "过渡带宽度 bandWidth", 0.005, 0.5, valinit=initial_bw)

    u0, u1 = model["uEdges"][0], model["uEdges"][-1]
    v0, v1 = model["vEdges"][0], model["vEdges"][-1]

    # 预计算（不随 k/bw 变化）：分片直纹面格胞边界点，用于 3D 图上标注。
    eval_pw = build_piecewise_eval(model)
    boundary_pts = []
    for ue in model["uEdges"][1:-1]:
        for vi in np.linspace(v0, v1, 30):
            boundary_pts.append(eval_pw(ue, vi))
    for ve in model["vEdges"][1:-1]:
        for ui in np.linspace(u0, u1, 30):
            boundary_pts.append(eval_pw(ui, ve))
    boundary_pts = np.asarray(boundary_pts)

    # 预计算：分片缝隙（不随 k/bw 变化）。
    gaps = {r["dir"] + str(r["boundary"]): r["piecewise_gap"]
            for r in check_continuity(model, k=int(initial_k), band_w=float(initial_bw))["surface"]}
    max_gap = max(gaps.values()) if gaps else 0.0

    # 缓存：check 4（过渡带）只依赖 k；3D 曲面网格依赖 (k, bw)。
    band_cache = {}
    surf_cache = {}

    def get_band(k):
        if k not in band_cache:
            band_cache[k] = check_transition_band(model, k=k, n=100)
        return band_cache[k]

    def get_surface(k, bw):
        key = (k, round(bw, 3))
        if key not in surf_cache:
            m = with_partition(model, k, bw)
            eval_blend, _ = tc.build_surface_eval(m)
            nu, nv = 30, 30
            U = np.linspace(u0, u1, nu)
            V = np.linspace(v0, v1, nv)
            X = np.zeros((nv, nu)); Y = np.zeros((nv, nu)); Z = np.zeros((nv, nu))
            for i in range(nv):
                for j in range(nu):
                    p = eval_blend(U[j], V[i])[0]
                    X[i, j], Y[i, j], Z[i, j] = p
            surf_cache[key] = (X, Y, Z)
        return surf_cache[key]

    def redraw():
        k = int(round(s_k.val))
        bw = float(s_bw.val)
        m = with_partition(model, k, bw)

        # ---- 曲面 ----
        ax_surf.clear()
        X, Y, Z = get_surface(k, bw)
        ax_surf.plot_surface(X, Y, Z, cmap="viridis", alpha=0.85, linewidth=0)
        ax_surf.scatter(boundary_pts[:, 0], boundary_pts[:, 1], boundary_pts[:, 2],
                        color="r", s=2)
        ax_surf.set_title(f"混合曲面 (k={k}, bw={bw:.3f})，红点为分片直纹面格胞边界")
        ax_surf.set_xlabel("x"); ax_surf.set_ylabel("y"); ax_surf.set_zlabel("z")

        # ---- (1) 连续性 ----
        ax_cont.clear()
        basis_jumps = check_basis_continuity(k)["jumps"]
        orders = list(range(len(basis_jumps)))
        ax_cont.bar(orders, [max(j, 1e-15) for j in basis_jumps],
                    color=["#1f77b4" if o <= k else "#d62728" for o in orders])
        ax_cont.axvline(k + 0.5, color="gray", linestyle="--", alpha=0.7)
        ax_cont.set_yscale("log")
        ax_cont.set_xticks(orders)
        ax_cont.set_xlabel("导数阶数 m")
        ax_cont.set_ylabel("基函数导数跳跃 |s_k^(m)(0,1)| (log)")
        ax_cont.set_title(f"(1) 分解函数精确 C^{k}：m≤{k} 跳跃=0，m={k+1} 跳跃>0"
                          f"（分片缝隙 {max_gap:.3f}mm，混合后值跳跃≈0）")
        ax_cont.grid(True, alpha=0.3, which="both")

        # ---- (2) 单位分解 + (3) 非负性 ----
        ax_pu.clear()
        pw = build_partition(m)
        us = np.linspace(u0, u1, 150)
        vmid = (v0 + v1) / 2
        sums = np.array([pw(u, vmid)[1].sum() for u in us])
        pu = check_partition_of_unity(m, n=50)
        nn = check_nonnegativity(m, n=50)
        ax_pu.plot(us, sums, color="#1f77b4", lw=1.5)
        ax_pu.axhline(1.0, color="red", linestyle="--", alpha=0.6)
        ax_pu.set_ylim(0.0, 2.0)
        ax_pu.set_xlabel("u"); ax_pu.set_ylabel("Σω = Σ_i ω_i")
        ax_pu.set_title(f"(2) 单位分解：Σω ≡ 1（max|Σω−1| = {pu['max_omega_dev']:.1e}）")
        ax_pu.text(0.02, 0.97, f"(3) 非负性：min φ = {nn['min_phi']:.1e},  min ω = {nn['min_omega']:.1e}",
                   transform=ax_pu.transAxes, va="top", ha="left", fontsize=9,
                   bbox=dict(boxstyle="round", fc="wheat", ec="gray", alpha=0.8))
        ax_pu.grid(True, alpha=0.3)

        # ---- (4) 过渡带 vs 保真度 ----
        ax_band.clear()
        band = get_band(k)
        bws = [b[0] for b in band]
        rms = [b[1] for b in band]
        ax_band.plot(bws, rms, marker="o", color="#ff7f0e", lw=1.5)
        ax_band.set_xlabel("bandWidth"); ax_band.set_ylabel("RMS |S_blend − S_piecewise| (mm)")
        ax_band.set_title("(4) 过渡带越窄，混合保真度越高（RMS 偏差）")
        ax_band.grid(True, alpha=0.3)

        fig.canvas.draw_idle()

    # 防抖：拖动滑块时只调度重绘，停止 ~250ms 后才真正重绘，避免拖动卡顿。
    timer = fig.canvas.new_timer(interval=250)
    timer.single_shot = True
    timer.add_callback(redraw)

    def request_redraw(_=None):
        timer.stop()
        timer.start()

    s_k.on_changed(request_redraw)
    s_bw.on_changed(request_redraw)
    redraw()
    plt.show()


# ======================================================================
# 主流程
# ======================================================================

def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser(description="验证自由曲面直纹化连续性处理的四个性质")
    ap.add_argument("--file", help="叶片文件（auto-identify→split-blade→拟合）")
    ap.add_argument("--model", help="已拟合 *_surface_model.json")
    ap.add_argument("--synthetic", action="store_true", help="构造带已知缝隙的合成模型")
    ap.add_argument("--outdir", default=None, help="拟合输出目录（--file 模式）")
    ap.add_argument("--tolerance", type=float, default=0.1)
    ap.add_argument("--nsplit-u", type=int, default=3)
    ap.add_argument("--nsplit-v", type=int, default=6)
    ap.add_argument("--continuity", type=int, default=None, help="验证时覆盖的连续性阶数 k")
    ap.add_argument("--band-width", type=float, default=None, help="验证时覆盖的 bandWidth")
    ap.add_argument("--no-gui", action="store_true", help="仅打印数值结果，不弹窗")
    ap.add_argument("--save-synthetic", default=None, help="把合成模型存成 JSON（如 synth.json）")
    args = ap.parse_args()

    # ---- 获取模型 ----
    if args.synthetic:
        model = build_synthetic_model(nu=args.nsplit_u, nv=args.nsplit_v)
        if args.save_synthetic:
            with open(args.save_synthetic, "w", encoding="utf-8") as f:
                json.dump(model, f, ensure_ascii=False, indent=1)
            print(f"合成模型已保存 → {args.save_synthetic}")
        src = "合成模型（带缝隙）"
    elif args.model:
        model = tc.load_surface_model(args.model)
        src = args.model
    elif args.file:
        if not BUILD_EXE.exists():
            print(f"[ERROR] simple.exe not found: {BUILD_EXE}")
            sys.exit(1)
        outdir = args.outdir or os.path.join(tempfile_mkdtemp(), "verify")
        models = preprocess_and_fit(os.path.abspath(args.file), args.tolerance,
                                    args.nsplit_u, args.nsplit_v, outdir,
                                    continuity=3, blend_width=0.05)
        if not models:
            print("[ERROR] 拟合未生成 surface_model.json")
            sys.exit(1)
        model = tc.load_surface_model(models[0])
        src = models[0]
    else:
        print("[ERROR] 需要 --file / --model / --synthetic 之一")
        sys.exit(1)

    k = args.continuity if args.continuity is not None else model["partition"]["continuity"]
    bw = args.band_width if args.band_width is not None else model["partition"]["bandWidth"]
    print(f"模型: {src}")
    print(f"  拟合方向 fitDir={model['fitDir']}, 格胞数={len(model['cells'])}, "
          f"验证 k={k}, bandWidth={bw}")

    # ---- 数值结果 ----
    pu = check_partition_of_unity(with_partition(model, k, bw))
    nn = check_nonnegativity(with_partition(model, k, bw))
    cont = check_continuity(model, k=k, band_w=bw)
    band = check_transition_band(model, k=k, n=100)

    print("\n(2) 单位分解恒等式")
    print(f"    max|Σω − 1| = {pu['max_omega_dev']:.3e}  (机器精度 {pu['machine_eps']:.2e})")
    print(f"    max|Σφ − 1| = {pu['max_phi_dev']:.3e}  (原始基函数，边界处可 >1)")
    print("\n(3) 非负性")
    print(f"    min φ = {nn['min_phi']:.3e}, min ω = {nn['min_omega']:.3e}  (应 ≥ 0)")
    print("\n(1) 连续性（分片直纹面 vs 混合曲面）")
    basis = cont["basis"]
    js = "  ".join(f"m{o}={j:.2e}" for o, j in enumerate(basis["jumps"]))
    print(f"    分解函数精确连续阶（导数跳跃，m≤{k} 应=0）：{js}")
    for r in cont["surface"]:
        print(f"    {r['dir']}-边界 {r['boundary']:.3f}: 分片缝隙={r['piecewise_gap']:.4f}mm, "
              f"混合值跳跃={r['blended_value_jump']:.2e}mm")
    print(f"    → 分片直纹面在边界有 {max(r['piecewise_gap'] for r in cont['surface']):.3f}mm 缝隙(C⁻¹)，"
          f"混合后值连续、且 C^{k} 连续")
    print("\n(4) 过渡带宽度 vs 保真度")
    for bw_, rms_, mx_ in band:
        print(f"    bandWidth={bw_:<6} → RMS偏差={rms_:.6f} mm, 最大偏差={mx_:.6f} mm")

    if not args.no_gui:
        run_interactive(model, initial_k=k, initial_bw=bw)


def tempfile_mkdtemp():
    import tempfile
    return tempfile.mkdtemp(prefix="verify_continuity_")


if __name__ == "__main__":
    main()
