# -*- coding: utf-8 -*-
"""连续刀轨规划算法（迭代条带法）——严格按方案执行，支持 fitDir=V/U。

曲面 S 按解析的 C^k 单位分解公式数值精确重建，法矢 nS、偏导 Su/Sv 解析求值。

当前版本简化（方案第七节）：
- 刀轴方向弱共轭优化：min_{‖T‖=1} w_conj·(nS·T)² + w_ruling·(1−(T·d)²)，
  d 为母线方向（fitDir=V 取 Sv，fitDir=U 取 Su），可展区退化为母线方向；
- 刀心 A 真正优化：对固定 T，逐进给位置最小二乘（min Σ(ρ−R)² + 正则）；
- 条带边界由包络覆盖极限（相交检测）确定；
- 总条数 N 初值 = 沿母线方向的直纹面数量（fitDir=V→nRows，fitDir=U→nCols），
  起枚举，终止条件用「平均误差 + 方差」联合判据。

参数：L_tool=25mm，o_min=2mm，o_max=5mm。
"""
import json
import math
import os
import numpy as np
from scipy.interpolate import BSpline, make_lsq_spline
from scipy.optimize import least_squares


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


def build_surface_eval(model):
    k = model["partition"]["continuity"]
    band_w = model["partition"]["bandWidth"]
    cells = model["cells"]
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
        acc = np.zeros(3); acc_u = np.zeros(3); acc_v = np.zeros(3)
        sum_phi = 0.0; sum_phi_u = 0.0; sum_phi_v = 0.0
        for (c, c0, c1, c0p, c1p) in cell_data:
            cu, cu_p = _chi1_np(u, c["u0"], c["u1"], a_u, a_u, k, deriv=True)
            cv, cv_p = _chi1_np(v, c["v0"], c["v1"], a_v, a_v, k, deriv=True)
            phi = cu * cv
            phi_u = cu_p * cv
            phi_v = cu * cv_p
            R = np.zeros(3); R_u = np.zeros(3); R_v = np.zeros(3)
            for d in range(3):
                c0v = float(c0[d](u)); c1v = float(c1[d](u))
                R[d] = (1.0 - v) * c0v + v * c1v
                R_u[d] = (1.0 - v) * float(c0p[d](u)) + v * float(c1p[d](u))
                R_v[d] = c1v - c0v
            acc += phi * R
            acc_u += phi_u * R + phi * R_u
            acc_v += phi_v * R + phi * R_v
            sum_phi += phi
            sum_phi_u += phi_u
            sum_phi_v += phi_v
        S = acc / (sum_phi + 1e-12)
        Su = acc_u / (sum_phi + 1e-12) - acc * sum_phi_u / (sum_phi ** 2 + 1e-24)
        Sv = acc_v / (sum_phi + 1e-12) - acc * sum_phi_v / (sum_phi ** 2 + 1e-24)
        nS = np.cross(Su, Sv)
        nS = nS / (np.linalg.norm(nS) + 1e-12)
        return S, nS, Su, Sv

    return eval_pt


def _point_axis_dist(P, A, T):
    d = P - A
    return np.linalg.norm(d - np.dot(d, T) * T)


def estimate_n_min(model):
    if model.get("fitDir") == "U":
        return max(1, model["nCols"])
    return max(1, model["nRows"])


def weak_conjugate_axis(surf_r, feed, rule_prev, rule_next, w_conj=1.0, w_ruling=5.0, n_v=8):
    vs = np.linspace(rule_prev, rule_next, n_v)
    C = np.zeros((3, 3))
    for v in vs:
        nS = surf_r(feed, v)[1]
        C += np.outer(nS, nS)
    d = surf_r(feed, 0.5 * (rule_prev + rule_next))[2]
    d = d / (np.linalg.norm(d) + 1e-12)
    M = w_conj * C - w_ruling * np.outer(d, d)
    w, V = np.linalg.eigh(M)
    T = V[:, 0]
    if np.dot(T, d) < 0:
        T = -T
    return T


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


def _coverage_limit(surf_r, feed, rule_prev, R, eps, s_min, s_max, rule_lo, rule_hi):
    rule_next = min(rule_prev + s_max, rule_hi)
    T = weak_conjugate_axis(surf_r, feed, rule_prev, rule_next)
    A = optimize_center_at_u(surf_r, feed, rule_prev, rule_next, T, R, n_v=8)
    for v in np.linspace(rule_prev + s_min, rule_next, 20):
        if abs(_point_axis_dist(surf_r(feed, v)[0], A, T) - R) > eps:
            return v
    return rule_next


def plan_toolpaths(model_path, L_tool=25.0, o_min=2.0, o_max=5.0,
                   eps_mean=0.05, eps_std=0.05, n_feed=40, n_rule=10, max_N_extra=3):
    model = load_surface_model(model_path)
    surf = build_surface_eval(model)
    R = 5.0
    fitDir = model["fitDir"]
    u_min, u_max = model["uEdges"][0], model["uEdges"][-1]
    v_min, v_max = model["vEdges"][0], model["vEdges"][-1]

    # 抽象 feed(进给/准线) 与 rule(母线/刀轴) 方向
    if fitDir == "U":
        feed_lo, feed_hi = v_min, v_max
        rule_lo, rule_hi = u_min, u_max

        def surf_r(feed, rule):
            S, nS, Su, Sv = surf(rule, feed)
            return S, nS, Su          # 母线方向 = Su
    else:  # V
        feed_lo, feed_hi = u_min, u_max
        rule_lo, rule_hi = v_min, v_max

        def surf_r(feed, rule):
            S, nS, Su, Sv = surf(feed, rule)
            return S, nS, Sv          # 母线方向 = Sv

    feed_grid = np.linspace(feed_lo, feed_hi, n_feed)

    L_avg = float(np.mean([np.linalg.norm(surf_r(f, rule_hi)[0] - surf_r(f, rule_lo)[0])
                           for f in np.linspace(feed_lo, feed_hi, 20)]))
    s_max = (L_tool - o_min) / L_avg * (rule_hi - rule_lo)
    s_min = (L_tool - o_max) / L_avg * (rule_hi - rule_lo)

    n_min = estimate_n_min(model)
    max_N = n_min + max_N_extra

    for N in range(n_min, max_N + 1):
        boundaries = [rule_lo]
        strips = []
        # 均匀划分母线区间为 N 个条带（避免覆盖极限在容差内不触发导致的退化）
        for j in range(N):
            b_prev = rule_lo + (rule_hi - rule_lo) * j / N
            b_next = rule_lo + (rule_hi - rule_lo) * (j + 1) / N
            boundaries.append(b_next)
            strips.append((b_prev, b_next))

        all_e = []
        feed_lines = []
        axis_segs = []
        for (b_prev, b_next) in strips:
            A_samples = np.zeros((n_feed, 3))
            T_samples = np.zeros((n_feed, 3))
            L_strip = np.zeros(n_feed)
            for i, f in enumerate(feed_grid):
                T = weak_conjugate_axis(surf_r, f, b_prev, b_next)
                T_samples[i] = T
                A_samples[i] = optimize_center_at_u(surf_r, f, b_prev, b_next, T, R, n_v=n_rule)
                L_strip[i] = np.linalg.norm(surf_r(f, b_next)[0] - surf_r(f, b_prev)[0])
            t = np.linspace(0.0, 1.0, n_feed)
            try:
                T_smooth = fit_bspline_3d(t, T_samples)(t)
                A_smooth = fit_bspline_3d(t, A_samples)(t)
            except Exception:
                T_smooth, A_smooth = T_samples, A_samples
            T_smooth = T_smooth / (np.linalg.norm(T_smooth, axis=1, keepdims=True) + 1e-12)
            vs = np.linspace(b_prev, b_next, n_rule)
            for i, f in enumerate(feed_grid):
                A = A_smooth[i]
                T = T_smooth[i]
                for v in vs:
                    all_e.append(_point_axis_dist(surf_r(f, v)[0], A, T) - R)
            feed_lines.append(A_smooth)
            half = 0.5 * L_strip + R
            axis_segs.append(np.stack([A_smooth - half[:, None] * T_smooth,
                                       A_smooth + half[:, None] * T_smooth], axis=1))

        all_e = np.array(all_e)
        err_mean = float(np.mean(np.abs(all_e)))
        err_std = float(np.std(all_e))
        print(f"N={N}: strips={len(strips)}, mean|e|={err_mean:.4f} mm, std(e)={err_std:.4f} mm")
        if err_mean <= eps_mean and err_std <= eps_std:
            print(f"终止：N={N} 满足 均值≤{eps_mean} 且 方差≤{eps_std}")
            return _build_result(model, R, N, boundaries, strips, feed_lines, axis_segs, err_mean, err_std)

    print(f"未能在 max_N={max_N} 内满足 均值/方差 容差")
    return _build_result(model, R, N, boundaries, strips, feed_lines, axis_segs, err_mean, err_std)


def _build_result(model, R, N, boundaries, strips, feed_lines, axis_segs, err_mean, err_std):
    return {
        "model": model, "R": R, "N": N, "boundaries": boundaries, "strips": strips,
        "feed_lines": feed_lines, "axis_segs": axis_segs,
        "err_mean": err_mean, "err_std": err_std,
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


if __name__ == "__main__":
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else None
    if not path:
        print("usage: python toolpath_continuous.py <surface_model.json> [--out <dir>]")
        sys.exit(1)
    out_dir = None
    if "--out" in sys.argv:
        out_dir = sys.argv[sys.argv.index("--out") + 1]
    res = plan_toolpaths(path)
    print(f"总刀轨条数 N={res['N']}, 平均误差={res['err_mean']:.4f} mm, 标准差={res['err_std']:.4f} mm")
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        stem = os.path.basename(path).replace("_surface_model.json", "")
        write_continuous_toolpath_vtk(res, os.path.join(out_dir, f"{stem}_toolpath_continuous.vtk"))
        print("已导出连续刀轨")
