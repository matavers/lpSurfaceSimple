# -*- coding: utf-8 -*-
"""连续刀轨规划算法（迭代条带法）——严格按方案执行，支持 fitDir=V/U。

曲面 S 按解析的 C^k 单位分解公式数值精确重建，法矢 nS、偏导 Su/Sv 解析求值。

当前版本简化（方案第七节 + 切削长度优化）：
- 刀轴方向直接取母线方向 T = d（d：fitDir=V 取 Sv，fitDir=U 取 Su），不优化；
- 刀心 A 真正优化：对固定 T，逐进给位置最小二乘（min Σ(ρ−R)² + 正则）；
- 切削长度 w 用梯度下降优化：最小化 (mean|e|−eps)²，求「误差≤eps 下的最大切削长度」，
  w 范围由重叠量 o_min 与刃长 L_tool 确定：[o_min, L_tool]；
- 各条带沿母线方向顺序铺满整面，条带间重叠 o_min。

参数由 research/toolpath_config.json 提供（tool_len/tool_r/o_min/o_max/eps 等）。
"""
import json
import math
import os
import sys
from pathlib import Path
import numpy as np
from scipy.interpolate import BSpline, make_lsq_spline
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parent))
import toolpath_config as tcfg
import machining_model as mm

import multiprocessing as _mp

# 并行 worker 的全局缓存（由 _pool_init 在每个子进程初始化）
_G_SURF = None
_G_RD = None
_G_FITDIR = None


def _pool_init(model):
    global _G_SURF, _G_RD, _G_FITDIR
    _G_SURF, _G_RD = build_surface_eval(model)
    _G_FITDIR = model["fitDir"]


def _scan_worker(payload):
    """并行 worker（模块级可 pickle）：payload=(b_prev, w, feed_grid, R, n_rule, n_v)。"""
    global _G_SURF, _G_RD, _G_FITDIR
    b_prev, w, feed_grid, R, n_rule, n_v = payload
    surf = _G_SURF
    rd = _G_RD
    if _G_FITDIR == "U":
        def surf_r(feed, rule):
            S, nS, _Su, _Sv = surf(rule, feed)
            return S, nS, rd(rule, feed)
        def surf_Sn(feed, rule):
            S, nS, _Su, _Sv = surf(rule, feed)
            return S, nS
    else:
        def surf_r(feed, rule):
            S, nS, _Su, _Sv = surf(feed, rule)
            return S, nS, rd(feed, rule)
        def surf_Sn(feed, rule):
            S, nS, _Su, _Sv = surf(feed, rule)
            return S, nS
    return _strip_error(surf_r, surf_Sn, b_prev, w, R, feed_grid, n_rule, n_v)


def _area_worker(args):
    """并行 worker：计算一行（固定 u，遍历 v）的面积微元。args=(v_min, n, u, du, dv)。"""
    global _G_SURF
    v_min, n, u, du, dv = args
    vs = v_min + dv * (np.arange(n) + 0.5)
    _S, _nS, Su, Sv = _G_SURF(np.full(n, u), vs)
    return float(np.sum(np.linalg.norm(np.cross(Su, Sv), axis=1)) * du * dv)


def _pose_worker(args):
    """并行 worker：单个进给位置的刀心/刀轴。args=(f, b_prev, b_next, R, n_rule)。"""
    global _G_SURF, _G_RD, _G_FITDIR
    f, b_prev, b_next, R, n_rule = args
    surf = _G_SURF
    rd = _G_RD
    if _G_FITDIR == "U":
        def surf_r(feed, rule):
            S, nS, _Su, _Sv = surf(rule, feed)
            return S, nS, rd(rule, feed)
    else:
        def surf_r(feed, rule):
            S, nS, _Su, _Sv = surf(feed, rule)
            return S, nS, rd(feed, rule)
    T = _ruling_direction(surf_r, f, b_prev, b_next)
    A = optimize_center_at_u(surf_r, f, b_prev, b_next, T, R, n_v=n_rule)
    L = np.linalg.norm(surf_r(f, b_next)[0] - surf_r(f, b_prev)[0])
    return A, T, L


def load_surface_model(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _smootherstep(t, k, deriv):
    if k <= 0:
        return t if not deriv else 1.0
    if k == 1:
        return 6.0 * t * (1.0 - t) if deriv else t * t * (3.0 - 2.0 * t)
    if k == 2:
        return 30.0 * t * t * (1.0 - t) ** 2 if deriv else t * t * t * (t * (6.0 * t - 15.0) + 10.0)
    if k == 3:
        return 140.0 * t ** 3 * (1.0 - t) ** 3 if deriv else t ** 4 * (35.0 - 84.0 * t + 70.0 * t * t - 20.0 * t ** 3)
    p = np.polynomial.Polynomial([0.0])
    for i in range(k + 1):
        p = p + math.comb(k + i, i) * np.polynomial.Polynomial([1.0, -1.0]) ** i
    p = p * np.polynomial.Polynomial([0.0] * (k + 1) + [1.0])
    return np.polynomial.Polynomial(p.coef).deriv()(t) if deriv else p(t)


def _chi1_np(u, uL, uR, a, b, k, deriv=False):
    if u <= uL - a:
        return (0.0, 0.0) if deriv else 0.0
    if u < uL:
        t = (u - (uL - a)) / a
        return (_smootherstep(t, k, False), _smootherstep(t, k, True) / a) if deriv else _smootherstep(t, k, False)
    if u <= uR:
        return (1.0, 0.0) if deriv else 1.0
    if u < uR + b:
        t = (u - uR) / b
        return (1.0 - _smootherstep(t, k, False), -_smootherstep(t, k, True) / b) if deriv else 1.0 - _smootherstep(t, k, False)
    return (0.0, 0.0) if deriv else 0.0


def _chi1_np_vec(u, uL, uR, a, b, k):
    """向量化版 _chi1_np（数组输入），返回 (chi, chi') 数组。"""
    u = np.asarray(u, dtype=float)
    chi = np.zeros_like(u)
    chi_p = np.zeros_like(u)
    # 上升过渡带 [uL-a, uL)
    m = (u > uL - a) & (u < uL)
    if np.any(m):
        t = (u[m] - (uL - a)) / a
        chi[m] = _smootherstep(t, k, False)
        chi_p[m] = (_smootherstep(t, k, True) if k > 0 else np.ones_like(t)) / a
    # 平台 [uL, uR]
    m = (u >= uL) & (u <= uR)
    chi[m] = 1.0
    # 下降过渡带 (uR, uR+b)
    m = (u > uR) & (u < uR + b)
    if np.any(m):
        t = (u[m] - uR) / b
        chi[m] = 1.0 - _smootherstep(t, k, False)
        chi_p[m] = -(_smootherstep(t, k, True) if k > 0 else np.ones_like(t)) / b
    return chi, chi_p


def build_surface_eval(model):
    k = model["partition"]["continuity"]
    band_w = model["partition"]["bandWidth"]
    cells = model["cells"]
    fitDir = model["fitDir"]
    min_u = min(c["u1"] - c["u0"] for c in cells)
    min_v = min(c["v1"] - c["v0"] for c in cells)
    a_u = band_w * min_u
    a_v = band_w * min_v
    cell_data = []
    for c in cells:
        c0 = [BSpline(c["c0"]["knots"], np.array(c["c0"]["ctrl"])[:, d], c["c0"]["degree"]) for d in range(3)]
        c1 = [BSpline(c["c1"]["knots"], np.array(c["c1"]["ctrl"])[:, d], c["c1"]["degree"]) for d in range(3)]
        c0p = [b.derivative() for b in c0]
        c1p = [b.derivative() for b in c1]
        cell_data.append((c, c0, c1, c0p, c1p))

    def eval_pt(u, v):
        """混合曲面解析求值。u/v 可为标量（返回 (S,nS,Su,Sv) 各为 3 向量）
        或等长数组（返回各为 (N,3) 数组），后者用于批量加速。"""
        u = np.atleast_1d(np.asarray(u, dtype=float)).ravel()
        v = np.atleast_1d(np.asarray(v, dtype=float)).ravel()
        N = u.size
        acc = np.zeros((N, 3)); acc_u = np.zeros((N, 3)); acc_v = np.zeros((N, 3))
        sum_phi = np.zeros(N); sum_phi_u = np.zeros(N); sum_phi_v = np.zeros(N)
        for (c, c0, c1, c0p, c1p) in cell_data:
            cu, cu_p = _chi1_np_vec(u, c["u0"], c["u1"], a_u, a_u, k)
            cv, cv_p = _chi1_np_vec(v, c["v0"], c["v1"], a_v, a_v, k)
            phi = cu * cv
            phi_u = cu_p * cv
            phi_v = cu * cv_p
            R = np.zeros((N, 3)); R_u = np.zeros((N, 3)); R_v = np.zeros((N, 3))
            du = c["u1"] - c["u0"]; dv = c["v1"] - c["v0"]
            for d in range(3):
                if fitDir == "U":
                    # 母线沿 u：R = (1-s)·c0(v) + s·c1(v)
                    c0v = c0[d](v); c1v = c1[d](v)
                    c0p_v = c0p[d](v); c1p_v = c1p[d](v)
                    s = (u - c["u0"]) / du
                    R[:, d] = (1.0 - s) * c0v + s * c1v
                    R_u[:, d] = (c1v - c0v) / du
                    R_v[:, d] = (1.0 - s) * c0p_v + s * c1p_v
                else:
                    # 母线沿 v：R = (1-s)·c0(u) + s·c1(u)
                    c0v = c0[d](u); c1v = c1[d](u)
                    c0p_v = c0p[d](u); c1p_v = c1p[d](u)
                    s = (v - c["v0"]) / dv
                    R[:, d] = (1.0 - s) * c0v + s * c1v
                    R_u[:, d] = (1.0 - s) * c0p_v + s * c1p_v
                    R_v[:, d] = (c1v - c0v) / dv
            acc += phi[:, None] * R
            acc_u += phi_u[:, None] * R + phi[:, None] * R_u
            acc_v += phi_v[:, None] * R + phi[:, None] * R_v
            sum_phi += phi
            sum_phi_u += phi_u
            sum_phi_v += phi_v
        denom = sum_phi + 1e-12
        S = acc / denom[:, None]
        Su = acc_u / denom[:, None] - acc * sum_phi_u[:, None] / (denom[:, None] ** 2)
        Sv = acc_v / denom[:, None] - acc * sum_phi_v[:, None] / (denom[:, None] ** 2)
        nS = np.cross(Su, Sv)
        nS = nS / (np.linalg.norm(nS, axis=1, keepdims=True) + 1e-12)
        if N == 1:
            return S[0], nS[0], Su[0], Sv[0]
        return S, nS, Su, Sv

    def ruling_dir(u, v):
        """母线方向 = 该格胞直纹面的 ruling（c1 − c0），即格内直纹面的直线方向。
        与分解函数求导无关，跨格不连续，但在格内为常方向（方案第七节「直接取母线方向」）。"""
        for (c, c0, c1, _c0p, _c1p) in cell_data:
            if c["u0"] - 1e-9 <= u <= c["u1"] + 1e-9 and c["v0"] - 1e-9 <= v <= c["v1"] + 1e-9:
                if fitDir == "U":
                    d = np.array([float(c1[d](v)) - float(c0[d](v)) for d in range(3)])
                else:
                    d = np.array([float(c1[d](u)) - float(c0[d](u)) for d in range(3)])
                n = np.linalg.norm(d)
                return d / n if n > 1e-12 else np.array([0.0, 0.0, 1.0])
        # 兜底：取最后一个格
        c, c0, c1, _c0p, _c1p = cell_data[-1]
        if fitDir == "U":
            d = np.array([float(c1[d](v)) - float(c0[d](v)) for d in range(3)])
        else:
            d = np.array([float(c1[d](u)) - float(c0[d](u)) for d in range(3)])
        n = np.linalg.norm(d)
        return d / n if n > 1e-12 else np.array([0.0, 0.0, 1.0])

    return eval_pt, ruling_dir


def _point_axis_dist(P, A, T):
    d = P - A
    return np.linalg.norm(d - np.dot(d, T) * T)


def _ruling_direction(surf_r, feed, rule_prev, rule_next):
    """刀轴方向直接取母线方向（方案第七节简化：T = d，d 为母线方向，不优化）。"""
    d = surf_r(feed, 0.5 * (rule_prev + rule_next))[2]
    return d / (np.linalg.norm(d) + 1e-12)


def _center_analytic(surf_r, feed, rule_prev, rule_next, R):
    """解析刀心：条带中点的曲面点沿法矢偏移 R（无需最小二乘，稳定快速，用于扫描）。"""
    S_mid = surf_r(feed, 0.5 * (rule_prev + rule_next))[0]
    n_mid = surf_r(feed, 0.5 * (rule_prev + rule_next))[1]
    return S_mid + R * n_mid


def optimize_center_at_u(surf_r, feed, rule_prev, rule_next, T, R, n_v=10, w_reg=0.2):
    vs = np.linspace(rule_prev, rule_next, n_v)
    P = np.array([surf_r(feed, v)[0] for v in vs])
    nS_mid = surf_r(feed, 0.5 * (rule_prev + rule_next))[1]
    A0 = P[n_v // 2] + R * nS_mid

    def residual(A):
        r = np.zeros(n_v + 3)
        for i in range(n_v):
            d = P[i] - A
            r[i] = np.linalg.norm(d - np.dot(d, T) * T) - R
        r[n_v:] = w_reg * (A - A0)
        return r

    return least_squares(residual, A0, max_nfev=30).x


def fit_bspline_3d(us, vals, n_ctrl=10, k=3):
    t = np.linspace(0.0, 1.0, len(us))
    t_inner = np.linspace(0.0, 1.0, n_ctrl - k + 1)[1:-1]
    knots = np.concatenate([np.zeros(k + 1), t_inner, np.ones(k + 1)])
    ctrl = np.zeros((n_ctrl, 3))
    for d in range(3):
        ctrl[:, d] = make_lsq_spline(t, vals[:, d], knots, k).c
    return BSpline(knots, ctrl, k)


def _optimize_strip_pose(surf_r, b_prev, b_next, R, feed_grid, n_rule, pool=None):
    """优化一条带 [b_prev, b_next] 的刀心 A（刀轴取母线方向），
    沿进给方向用 B 样条光顺，返回 (A_smooth, T_smooth, L_strip)。逐进给位置可并行。"""
    n_feed = len(feed_grid)
    if pool is not None:
        results = pool.map(_pose_worker,
                           [(f, b_prev, b_next, R, n_rule) for f in feed_grid])
        A_samples = np.array([r[0] for r in results])
        T_samples = np.array([r[1] for r in results])
        L_strip = np.array([r[2] for r in results])
    else:
        A_samples = np.zeros((n_feed, 3))
        T_samples = np.zeros((n_feed, 3))
        L_strip = np.zeros(n_feed)
        for i, f in enumerate(feed_grid):
            T = _ruling_direction(surf_r, f, b_prev, b_next)
            A = optimize_center_at_u(surf_r, f, b_prev, b_next, T, R, n_v=n_rule)
            T_samples[i] = T
            A_samples[i] = A
            L_strip[i] = np.linalg.norm(surf_r(f, b_next)[0] - surf_r(f, b_prev)[0])
    t = np.linspace(0.0, 1.0, n_feed)
    try:
        T_smooth = fit_bspline_3d(t, T_samples)(t)
        A_smooth = fit_bspline_3d(t, A_samples)(t)
    except Exception:
        T_smooth, A_smooth = T_samples, A_samples
    T_smooth = T_smooth / (np.linalg.norm(T_smooth, axis=1, keepdims=True) + 1e-12)
    return A_smooth, T_smooth, L_strip


def _strip_error(surf_r, surf_Sn, b_prev, w, R, feed_grid, n_rule, n_v=8):
    """给定切削长度 w（参数域），用解析刀心 A（稳定快速），返回 (mean|e|, max|e|)。
    用全部进给位置 + 母线子采样，保证最大误差（过切/欠切极值）不被采样遗漏。
    surf_Sn(feed, rule) 为批量求点接口（数组输入，返回 (S,nS) 各 (N,3)），用于向量化加速。"""
    b_next = b_prev + w
    mid = 0.5 * (b_prev + b_next)
    vs = np.linspace(b_prev, b_next, n_v)
    FF, VV = np.meshgrid(feed_grid, vs, indexing="ij")
    S_all, _ = surf_Sn(FF.ravel(), VV.ravel())            # (N,3) 全部采样点
    S_mid, n_mid = surf_Sn(feed_grid, np.full_like(feed_grid, mid))  # (n_feed,3) 中点+法矢
    es = np.empty(FF.size)
    idx = 0
    for i, f in enumerate(feed_grid):
        T = _ruling_direction(surf_r, f, b_prev, b_next)
        A = S_mid[i] + R * n_mid[i]
        for _v in vs:
            es[idx] = abs(_point_axis_dist(S_all[idx], A, T) - R)
            idx += 1
    return float(es.mean()), float(es.max())


def _optimize_strip_width(surf_r, surf_Sn, b_prev, R, eps, w_min, w_max, rule_hi, feed_grid,
                          n_rule, n_v=8, n_scan=40, pool=None, w_hint=None):
    """求「最大误差≤eps 下的最大切削长度 w」（方案 B，目标用 max|e| 保证不过切）。
    误差对 w 非单调（曲面存在特征导致误差跳变），故用线性扫描从 w_max 向下找第一个满足 max|e|≤eps 的 w，
    比梯度下降更稳健。扫描点可并行（pool 进程池）。返回 (w, mean|e|, max|e|)。

    加速：传 w_hint（上一条刀轨的切削长度）时，先在其附近（0~2×w_hint）扫描，
    扫描点数按范围等比例缩减，避免紧容差下从满刃长 w_max 一路扫到很小的 w 浪费大量求值。"""
    w_max = min(w_max, rule_hi - b_prev)
    if w_max <= w_min:
        w = w_max
        mean_e, max_e = _strip_error(surf_r, surf_Sn, b_prev, w, R, feed_grid, n_rule, n_v)
        return w, mean_e, max_e

    def _scan(ws):
        if pool is not None:
            batch = max(4, int(getattr(pool, "_processes", 4)))
            for i in range(0, len(ws), batch):
                batch_ws = ws[i:i + batch]
                payloads = [(b_prev, w, feed_grid, R, n_rule, n_v) for w in batch_ws]
                results = pool.map(_scan_worker, payloads)
                for w, (mean_e, max_e) in zip(batch_ws, results):
                    if max_e <= eps:
                        return (w, mean_e, max_e)
        else:
            for w in ws:
                mean_e, max_e = _strip_error(surf_r, surf_Sn, b_prev, w, R, feed_grid, n_rule, n_v)
                if max_e <= eps:
                    return (w, mean_e, max_e)
        return None

    # 全扫描使用固定网格（保证结果与无热启动完全一致）。
    ws_full = np.linspace(w_max, w_min, n_scan)

    # 热启动：从上一条刀轨宽度的 2 倍处开始，跳过更宽的候选（网格不变，仅截断起点）。
    start = 0
    if w_hint is not None:
        scan_hi = min(w_max, max(w_min, w_hint * 2.0))
        start = int(np.argmax(ws_full <= scan_hi))  # ws_full 降序，取第一个 ≤ scan_hi 的下标

    found = _scan(ws_full[start:])
    if found is not None:
        w_found = found[0]
        # 命中的是热启动起点（上边界）=> 真实最优可能在更宽处，补扫被跳过的上段
        if start > 0 and w_found >= ws_full[start] - 1e-12:
            top = _scan(ws_full[:start])
            if top is not None:
                return top
        return found

    # 全部都不满足（eps 太紧），取最小宽度
    w = w_min
    mean_e, max_e = _strip_error(surf_r, surf_Sn, b_prev, w, R, feed_grid, n_rule, n_v)
    return w, mean_e, max_e


def _surface_area_mm2(surf, u_min, u_max, v_min, v_max, n=120, pool=None):
    """数值积分曲面面积：∫∫ |Su×Sv| du dv（mm²）。可并行（pool 按行切分）。"""
    du = (u_max - u_min) / n
    dv = (v_max - v_min) / n
    if pool is not None:
        tasks = [(v_min, n, u_min + du * (i + 0.5), du, dv) for i in range(n)]
        area = float(sum(pool.map(_area_worker, tasks)))
    else:
        us = u_min + du * (np.arange(n) + 0.5)
        vs = v_min + dv * (np.arange(n) + 0.5)
        UU, VV = np.meshgrid(us, vs, indexing="ij")
        _S, _nS, Su, Sv = surf(UU.ravel(), VV.ravel())
        area = float(np.sum(np.linalg.norm(np.cross(Su, Sv), axis=1)) * du * dv)
    return area


def plan_toolpaths(model_path, L_tool=25.0, R=5.0, o_min=2.0, o_max=5.0,
                   eps=0.1, n_feed=None, n_rule=10, n_jobs=1):
    """迭代条带法（方案第七节 + 切削长度优化）：
    - 刀轴方向直接取母线方向 T = d（不优化）；
    - 刀心 A 逐进给位置最小二乘优化；
    - 切削长度 w 用梯度下降优化：最小化 (mean|e|−eps)²，求「误差≤eps 下的最大切削长度」；
    - w 范围由重叠量 o_min 与刃长 L_tool 确定：[o_min, L_tool]；
    - 各条带沿母线顺序铺满整面，重叠 o_min；
    - 切削长度线性扫描可并行（n_jobs 进程数，joblib）。"""
    model = load_surface_model(model_path)
    surf, ruling_dir = build_surface_eval(model)
    fitDir = model["fitDir"]
    u_min, u_max = model["uEdges"][0], model["uEdges"][-1]
    v_min, v_max = model["vEdges"][0], model["vEdges"][-1]

    # 抽象 feed(进给/准线) 与 rule(母线/刀轴) 方向
    if fitDir == "U":
        feed_lo, feed_hi = v_min, v_max
        rule_lo, rule_hi = u_min, u_max

        def surf_r(feed, rule):
            S, nS, _Su, _Sv = surf(rule, feed)
            return S, nS, ruling_dir(rule, feed)   # 母线方向 = 格内直纹面 ruling

        def surf_Sn(feed, rule):
            S, nS, _Su, _Sv = surf(rule, feed)
            return S, nS
    else:  # V
        feed_lo, feed_hi = u_min, u_max
        rule_lo, rule_hi = v_min, v_max

        def surf_r(feed, rule):
            S, nS, _Su, _Sv = surf(feed, rule)
            return S, nS, ruling_dir(feed, rule)   # 母线方向 = 格内直纹面 ruling

        def surf_Sn(feed, rule):
            S, nS, _Su, _Sv = surf(feed, rule)
            return S, nS

    # 按物理进给长度自适应采样数，使各面刀轴间距一致（避免短弦向进给时刀轴过密）。
    _fs = np.linspace(feed_lo, feed_hi, 40)
    _prev = surf_r(_fs[0], 0.5 * (rule_lo + rule_hi))[0]
    _L_feed = 0.0
    for _f in _fs[1:]:
        _p = surf_r(_f, 0.5 * (rule_lo + rule_hi))[0]
        _L_feed += float(np.linalg.norm(_p - _prev))
        _prev = _p
    if n_feed is None:
        n_feed = max(20, min(240, int(round(_L_feed / 1.5)) + 1))
    feed_grid = np.linspace(feed_lo, feed_hi, n_feed)

    # 母线长（整面最大 W）、最大切削宽度/重叠量（参数域）
    L_ruling = np.array([np.linalg.norm(surf_r(f, rule_hi)[0] - surf_r(f, rule_lo)[0])
                         for f in np.linspace(feed_lo, feed_hi, 40)])
    L_avg = float(L_ruling.mean())
    W = float(L_ruling.max())
    rule_range = rule_hi - rule_lo
    # 切削长度 w 的范围由重叠量与刃长确定：[o_min, L_tool]（参数域）
    w_min = o_min / L_avg * rule_range
    w_max = L_tool / L_avg * rule_range
    o_param = o_min / L_avg * rule_range  # 重叠量（参数域），取 o_min
    print(f"切削长度范围 w∈[{w_min:.4f}, {w_max:.4f}] (param)，W={W:.2f}mm，eps={eps}mm")

    boundaries = [rule_lo]
    strips = []
    feed_lines = []
    axis_segs = []
    intersection_lines = []
    strip_data = []
    b_prev = rule_lo
    max_strips = 200  # 安全上限

    # 并行进程池（一次创建、所有条带复用，避免逐条带重复 spawn 开销）
    pool = None
    if n_jobs is None or n_jobs > 1 or n_jobs < 0:
        if n_jobs and n_jobs > 1:
            n_proc = min(n_jobs, _mp.cpu_count())
        else:
            n_proc = min(8, _mp.cpu_count())
        try:
            pool = _mp.Pool(processes=n_proc, initializer=_pool_init, initargs=(model,))
            print(f"并行计算已启用：{n_proc} 进程", flush=True)
        except Exception:
            pool = None

    w_prev = None  # 上一条刀轨的切削长度，用于宽度扫描的热启动
    while b_prev < rule_hi - 1e-9 and len(strips) < max_strips:
        w, mean_e, max_e = _optimize_strip_width(
            surf_r, surf_Sn, b_prev, R, eps, w_min, w_max, rule_hi, feed_grid, n_rule,
            pool=pool, w_hint=w_prev)
        if w <= 1e-9:
            print(f"梯度下降收敛到 w≈0（eps={eps} 下该处无法侧铣），停止于 {len(strips)} 条带", flush=True)
            break
        w_prev = w
        b_next = b_prev + w
        intersection_lines.append(np.array([surf_r(f, b_next)[0] for f in feed_grid]))

        A_smooth, T_smooth, L_strip = _optimize_strip_pose(
            surf_r, b_prev, b_next, R, feed_grid, n_rule, pool=pool)
        feed_lines.append(A_smooth)
        half = 0.5 * L_tool + R  # 刀轴线段 = 刀具刃长 + 2R（方案 2.3）
        axis_segs.append(np.stack([A_smooth - half * T_smooth,
                                   A_smooth + half * T_smooth], axis=1))
        strips.append((b_prev, b_next))
        boundaries.append(b_next)
        strip_data.append((b_prev, b_next, A_smooth, T_smooth))

        w_mm = w * L_avg / rule_range
        print(f"[条带 {len(strips)}] 切削长度 w={w_mm:.2f}mm, mean|e|={mean_e:.4f}mm, max|e|={max_e:.4f}mm "
              f"(范围[{w_min * L_avg / rule_range:.2f}, {w_max * L_avg / rule_range:.2f}]mm, "
              f"边界 v={b_prev:.4f}→{b_next:.4f})", flush=True)

        if b_next >= rule_hi - 1e-9:
            break
        # 下一条带确定边界 = 切削终点 − 重叠量；若切削长度已到最小(≈o_min)，
        # 无法再重叠，则零重叠推进以继续覆盖、避免推进≈0 的死循环
        if w > o_param:
            b_prev = b_next - o_param
        else:
            if max_e > eps:
                # 只有当真未达标（max|e|>eps）才警告 eps 太紧；若已达准则静默推进
                print(f"  警告：eps={eps} 太紧，最小切削长度 {w_min * L_avg / rule_range:.2f}mm 下 "
                      f"max|e|={max_e:.4f} 仍>eps，零重叠推进", flush=True)
            b_prev = b_next

    all_e = []
    for (b0, b1, A, T) in strip_data:
        for i, f in enumerate(feed_grid):
            for v in np.linspace(b0, b1, n_rule):
                all_e.append(_point_axis_dist(surf_r(f, v)[0], A[i], T[i]) - R)
    all_e = np.array(all_e)
    err_mean = float(np.mean(np.abs(all_e)))
    err_std = float(np.std(all_e))
    N = len(strips)
    if strips:
        avg_width_mm = float(np.mean([b1 - b0 for b0, b1 in strips]) * L_avg / rule_range)
    else:
        avg_width_mm = 0.0
    print(f"总条带数 N={N}, 平均切削宽度={avg_width_mm:.2f}mm, "
          f"mean|e|={err_mean:.4f} mm, std(e)={err_std:.4f} mm")

    # 加工时间/载荷估计（工程公式）：刀轨总长 + 各条带切削长度 + 曲面面积
    mcfg = mm.load_config()
    path_len_mm = float(sum(np.linalg.norm(np.diff(fl, axis=0), axis=1).sum()
                            for fl in feed_lines if len(fl) >= 2))
    cut_lens_mm = [(b1 - b0) * L_avg / rule_range for (b0, b1) in strips]
    flank = mm.estimate_flank(path_len_mm, cut_lens_mm, mcfg, R)

    # 曲面面积（数值积分 |Su×Sv|）→ 点铣基线
    area_mm2 = float(_surface_area_mm2(surf, u_min, u_max, v_min, v_max, pool=pool))
    # 点铣精度与侧铣对应：残留高度 h = 侧铣可接受误差 eps（同精度下对比才公平）
    point_cfg = dict(mcfg)
    point_cfg["scallop"] = eps
    point = mm.estimate_point(area_mm2, point_cfg, mcfg.get("ball_r", 5.0))

    if pool is not None:
        pool.close()
        pool.join()

    speedup = point["cut_time_s"] / flank["cut_time_s"] if flank["cut_time_s"] > 0 else 0.0
    mach = {
        "flank": flank, "point": point, "speedup": speedup,
        "flank_time_s": flank["cut_time_s"], "point_time_s": point["cut_time_s"],
        "num_strips": N, "err_mean": err_mean, "err_std": err_std,
        "path_length_mm": path_len_mm, "surface_area_mm2": area_mm2,
        "avg_cutting_length_mm": avg_width_mm,
        "config": {k: mcfg.get(k) for k in ("cutting_speed", "feed_per_tooth",
                                            "num_teeth", "specific_cutting_force",
                                            "radial_depth", "tool_r", "ball_r", "scallop")},
    }

    print(f"加工估计: 侧铣 转速={flank['spindle_rpm']:.0f}rpm 进给率={flank['feed_rate_mm_min']:.0f}mm/min "
          f"时间={flank['cut_time_s']:.1f}s MRR={flank['mrr_cm3_min']:.2f}cm3/min "
          f"力={flank['cutting_force_N']:.0f}N 功率={flank['cutting_power_kW']:.2f}kW; "
          f"点铣(残留={eps}mm) 时间={point['cut_time_s']:.1f}s 行距={point['stepover_mm']:.2f}mm; "
          f"提速={speedup:.1f}x", flush=True)

    return _build_result(model, R, N, boundaries, strips, feed_lines, axis_segs,
                         intersection_lines, err_mean, err_std, mach)


def _build_result(model, R, N, boundaries, strips, feed_lines, axis_segs,
                  intersection_lines, err_mean, err_std, mach=None):
    return {
        "model": model, "R": R, "N": N, "boundaries": boundaries, "strips": strips,
        "feed_lines": feed_lines, "axis_segs": axis_segs,
        "intersection_lines": intersection_lines,
        "err_mean": err_mean, "err_std": err_std,
        "machining": mach,
        "meta": {"fitDir": model["fitDir"],
                 "continuity": model["partition"]["continuity"],
                 "n_cells": len(model["cells"])},
    }


def write_continuous_toolpath_vtk(res, out_path):
    # 每条带进给线独立 polyline（不串联），刀轴线段独立 2 点线，消除条带间斜向跳变
    feed_lines = [np.asarray(fl) for fl in res["feed_lines"]]
    axis_lines = [np.asarray([s[0], s[1]]) for segs in res["axis_segs"] for s in segs]
    lines = feed_lines + axis_lines
    n_pts = sum(len(l) for l in lines)
    n_segs = sum(len(l) - 1 for l in lines if len(l) >= 2)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("# vtk DataFile Version 3.0\ncontinuous toolpath\nASCII\nDATASET POLYDATA\n")
        f.write(f"POINTS {n_pts} float\n")
        for l in lines:
            for p in l:
                f.write(f"{p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
        f.write(f"LINES {n_segs} {n_segs * 3}\n")
        base = 0
        for l in lines:
            for i in range(len(l) - 1):
                f.write(f"2 {base + i} {base + i + 1}\n")
            base += len(l)
    for suffix, seg_list in (("_feed", feed_lines), ("_axis", axis_lines)):
        with open(out_path.replace(".vtk", suffix + ".vtk"), "w", encoding="utf-8") as f:
            n_pts = sum(len(l) for l in seg_list)
            n_segs = sum(len(l) - 1 for l in seg_list if len(l) >= 2)
            f.write("# vtk DataFile Version 3.0\ntoolpath\nASCII\nDATASET POLYDATA\n")
            f.write(f"POINTS {n_pts} float\n")
            for l in seg_list:
                for p in l:
                    f.write(f"{p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
            f.write(f"LINES {n_segs} {n_segs * 3}\n")
            base = 0
            for l in seg_list:
                for i in range(len(l) - 1):
                    f.write(f"2 {base + i} {base + i + 1}\n")
                base += len(l)
    inter_lines = [np.asarray(l) for l in res.get("intersection_lines", [])
                   if len(np.asarray(l)) >= 2]
    if inter_lines:
        with open(out_path.replace(".vtk", "_intersections.vtk"), "w", encoding="utf-8") as f:
            n_pts = sum(len(l) for l in inter_lines)
            n_segs = sum(len(l) - 1 for l in inter_lines)
            f.write("# vtk DataFile Version 3.0\nintersections\nASCII\nDATASET POLYDATA\n")
            f.write(f"POINTS {n_pts} float\n")
            for l in inter_lines:
                for p in l:
                    f.write(f"{p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
            f.write(f"LINES {n_segs} {n_segs * 3}\n")
            base = 0
            for l in inter_lines:
                for i in range(len(l) - 1):
                    f.write(f"2 {base + i} {base + i + 1}\n")
                base += len(l)


if __name__ == "__main__":
    import sys
    import argparse

    # 统一 stdout/stderr 为 UTF-8，避免 Windows 控制台 GBK 与父进程(UTF-8)读取不一致导致中文乱码
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    if hasattr(sys.stderr, "reconfigure"):
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")

    cfg = tcfg.load_config()
    ap = argparse.ArgumentParser(description="连续刀轨规划（迭代条带法）")
    ap.add_argument("model", help="*_surface_model.json 路径")
    ap.add_argument("--out", help="输出目录")
    ap.add_argument("--tool-len", type=float, default=cfg.get("tool_len"), help="刀具有效切削刃长 L_tool (mm)")
    ap.add_argument("--tool-r", type=float, default=cfg.get("tool_r"), help="刀具半径 R (mm)")
    ap.add_argument("--o-min", type=float, default=cfg.get("o_min"), help="最小重叠量 o_min (mm)")
    ap.add_argument("--o-max", type=float, default=cfg.get("o_max"), help="最大重叠量 o_max (mm)")
    ap.add_argument("--eps", type=float, default=cfg.get("eps"), help="可接受误差 mean|e| (mm)")
    ap.add_argument("--n-rule", type=int, default=cfg.get("n_rule"), help="母线方向采样数")
    ap.add_argument("--n-jobs", type=int, default=cfg.get("n_jobs", 1), help="并行进程数（-1=全部CPU，1=串行）")
    args = ap.parse_args()

    stem = os.path.basename(args.model).replace("_surface_model.json", "")
    res = plan_toolpaths(args.model, L_tool=args.tool_len, R=args.tool_r,
                         o_min=args.o_min, o_max=args.o_max,
                         eps=args.eps, n_rule=args.n_rule, n_jobs=args.n_jobs)
    print(f"总刀轨条数 N={res['N']}, 平均误差={res['err_mean']:.4f} mm, 标准差={res['err_std']:.4f} mm")
    if args.out:
        os.makedirs(args.out, exist_ok=True)
        write_continuous_toolpath_vtk(res, os.path.join(args.out, f"{stem}_toolpath_continuous.vtk"))
        if res.get("machining"):
            mp = os.path.join(args.out, f"{stem}_machining_summary.json")
            with open(mp, "w", encoding="utf-8") as f:
                json.dump(res["machining"], f, ensure_ascii=False, indent=2)
            print(f"已导出加工参数估计 {mp}")
        print("已导出连续刀轨")
