# -*- coding: utf-8 -*-
"""在完整 C^k 连续曲面上做可微刀位规划（CasADi + IPOPT）。

关键：不把曲面降维成点云。曲面 S 按解析的 C^k 单位分解公式精确重建，
法矢 nS 用解析导数（B 样条导数 + smootherstep 导数）在求积点上求值；
刀轴场 T(u) 用 B 样条参数化（控制点可学习），目标泛函
  J = ∫[ (nS·T)^2 + w·‖T'‖^2 ] du
在 Gauss-Legendre 求积点上求值，用 IPOPT 梯度下降。
连续性体现在：S∈C^k ⇒ nS∈C^{k-1} 连续 ⇒ 目标泛函连续可微 ⇒ 精确梯度。
"""
import json
import math
import numpy as np
from scipy.interpolate import BSpline
import casadi as ca


def load_surface_model(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


# ── 分解函数（smootherstep 及其导数，numpy） ─────────────────
def _smootherstep_poly(k):
    """返回 s_k(t) 在 t 的多项式系数（升幂），以及其导数系数。"""
    if k <= 0:
        return np.array([0.0, 1.0]), np.array([1.0])
    # s_k(t) = t^{k+1} * sum_i C(k+i,i) (1-t)^i
    # 展开成多项式：先用 numpy 多项式乘法
    coeff = np.array([1.0])
    for i in range(k + 1):
        # (1-t)^i 的系数（升幂）
        binom = math.comb(k + i, i)
        poly = np.poly1d([1.0, -1.0]) ** i
        coeff = np.polynomial.polynomial.polyadd(
            coeff, binom * np.polynomial.polynomial.polyfromroots(np.zeros(1)) * 0)
        # 累加：直接构造 (1-t)^i 系数向量
    # 用多项式对象更简单
    p = np.polynomial.Polynomial([0.0])
    for i in range(k + 1):
        binom = math.comb(k + i, i)
        omt = np.polynomial.Polynomial([1.0, -1.0]) ** i
        p = p + binom * omt
    # 乘 t^{k+1}
    t_pow = np.polynomial.Polynomial([0.0] * (k + 1) + [1.0])
    p = p * t_pow
    c = p.coef
    dc = np.polynomial.Polynomial(c).deriv().coef
    return c, dc


def _chi1_np(u, uL, uR, a, b, k, deriv=False):
    """一维单位分解基函数 chi 及其导数（numpy 标量）。"""
    c, dc = _smootherstep_poly(k)

    def polyval(coeff, x):
        return np.polynomial.Polynomial(coeff)(x)

    if u <= uL - a:
        return (0.0, 0.0) if deriv else 0.0
    if u < uL:
        t = (u - (uL - a)) / a
        val = polyval(c, t)
        dval = polyval(dc, t) / a
        return (val, dval) if deriv else val
    if u <= uR:
        return (1.0, 0.0) if deriv else 1.0
    if u < uR + b:
        t = (u - uR) / b
        val = 1.0 - polyval(c, t)
        dval = -polyval(dc, t) / b
        return (val, dval) if deriv else val
    return (0.0, 0.0) if deriv else 0.0


# ── 曲面数值重建（解析，含导数） ─────────────────────────────
def build_surface_eval(model):
    """返回 eval(u, v) -> (S, nS)，numpy 数值（S、nS 均为解析求值，非点云）。"""
    k = model["partition"]["continuity"]
    band_w = model["partition"]["bandWidth"]
    cells = model["cells"]
    min_u = min(c["u1"] - c["u0"] for c in cells)
    min_v = min(c["v1"] - c["v0"] for c in cells)
    a_u = band_w * min_u
    a_v = band_w * min_v

    cell_data = []
    for c in cells:
        c0 = [BSpline(c["c0"]["knots"], np.array(c["c0"]["ctrl"])[:, d], c["c0"]["degree"])
              for d in range(3)]
        c1 = [BSpline(c["c1"]["knots"], np.array(c["c1"]["ctrl"])[:, d], c["c1"]["degree"])
              for d in range(3)]
        cell_data.append((c, c0, c1))

    def eval_pt(u, v):
        acc = np.zeros(3)
        acc_u = np.zeros(3)
        acc_v = np.zeros(3)
        sum_phi = 0.0
        sum_phi_u = 0.0
        sum_phi_v = 0.0
        for (c, c0, c1) in cell_data:
            cu, cu_p = _chi1_np(u, c["u0"], c["u1"], a_u, a_u, k, deriv=True)
            cv, cv_p = _chi1_np(v, c["v0"], c["v1"], a_v, a_v, k, deriv=True)
            phi = cu * cv
            phi_u = cu_p * cv
            phi_v = cu * cv_p
            R = np.zeros(3)
            R_u = np.zeros(3)
            R_v = np.zeros(3)
            for d in range(3):
                c0v = float(c0[d](u))
                c1v = float(c1[d](u))
                c0p = float(c0[d].derivative()(u))
                c1p = float(c1[d].derivative()(u))
                R[d] = (1.0 - v) * c0v + v * c1v
                R_u[d] = (1.0 - v) * c0p + v * c1p
                R_v[d] = c1v - c0v
            acc += phi * R
            acc_u += phi_u * R + phi * R_u
            acc_v += phi_v * R + phi * R_v
            sum_phi += phi
            sum_phi_u += phi_u
            sum_phi_v += phi_v
        S = acc / sum_phi
        Su = acc_u / sum_phi - acc * sum_phi_u / (sum_phi ** 2)
        Sv = acc_v / sum_phi - acc * sum_phi_v / (sum_phi ** 2)
        nS = np.cross(Su, Sv)
        norm = np.linalg.norm(nS)
        nS = nS / (norm + 1e-12)
        return S, nS

    return eval_pt


# ── Cox-de Boor B 样条（CasADi MX，可微；规避 ca.bspline 对系数的 AD 缺陷） ──
def _bspline_mx(u, P, knots, degree):
    """3D B 样条曲线 C(u)=Σ N_{i,p}(u) P_i，P 为 (n_ctrl,3) MX。"""
    from functools import lru_cache
    n = len(knots) - degree - 1

    @lru_cache(maxsize=None)
    def basis(i, p):
        if p == 0:
            return ca.if_else(ca.logic_and(u >= knots[i], u < knots[i + 1]), 1.0, 0.0)
        d1 = knots[i + p] - knots[i]
        d2 = knots[i + p + 1] - knots[i + 1]
        a = 0.0 if d1 <= 1e-12 else (u - knots[i]) / d1 * basis(i, p - 1)
        b = 0.0 if d2 <= 1e-12 else (knots[i + p + 1] - u) / d2 * basis(i + 1, p - 1)
        return a + b

    N = [basis(i, degree) for i in range(n)]
    return ca.vertcat(*[sum(N[i] * P[i, d] for i in range(n)) for d in range(3)])


# ── 刀轴场优化（CasADi） ─────────────────────────────────────
def optimize_tool_axis_field(model_path, n_ctrl=24, w_smooth=0.5,
                             n_quad=100, solver="ipopt", continuity=None,
                             tool_r=5.0):
    model = load_surface_model(model_path)
    surf = build_surface_eval(model)
    k = continuity if continuity is not None else model["partition"]["continuity"]

    u_min = model["uEdges"][0]
    u_max = model["uEdges"][-1]
    v_mid = 0.5 * (model["vEdges"][0] + model["vEdges"][-1])

    # Gauss-Legendre 求积点
    pts, wts = np.polynomial.legendre.leggauss(n_quad)
    us = 0.5 * (u_max - u_min) * pts + 0.5 * (u_max + u_min)
    wts = wts * 0.5 * (u_max - u_min)

    # 曲面法矢在求积点的值（数值，解析求值）
    nS_grid = np.array([surf(ui, v_mid)[1] for ui in us])   # (n_quad, 3)

    # 刀轴场 B 样条（degree = k+1，C^k 连续）
    degree = k + 1
    n_ctrl = max(degree + 1, n_ctrl)
    inner = n_ctrl - degree - 1
    knots = [u_min] * (degree + 1) + \
            [u_min + (u_max - u_min) * (i + 1) / (inner + 1) for i in range(inner)] + \
            [u_max] * (degree + 1)

    u = ca.MX.sym("u")
    pvec = ca.MX.sym("p", 3 * n_ctrl)
    P = ca.reshape(pvec, n_ctrl, 3)   # (n_ctrl,3)，列主序（CasADi reshape）

    T_raw = _bspline_mx(u, P, knots, degree)   # (3,1)

    # 解析 B 样条导数 T'(u)（B 样条导数公式，避免对参数 u 求导的 AD）
    if degree >= 1:
        inner_knots = knots[1:-1]          # 去掉首尾（导数节点）
        Pderiv = ca.vertcat(*[
            (P[j + 1, :] - P[j, :]) * (degree / (knots[j + degree + 1] - knots[j + 1])
                                        if knots[j + degree + 1] - knots[j + 1] > 1e-12 else 0.0)
            for j in range(n_ctrl - 1)
        ])   # (n_ctrl-1, 3)
        Tp_raw = _bspline_mx(u, Pderiv, inner_knots, degree - 1)
    else:
        Tp_raw = ca.MX.zeros(3, 1)

    # 目标泛函在求积点上的离散和（T 不归一化，用罚项约束 ‖T‖≈1）
    conj = 0.0
    smooth = 0.0
    norm_pen = 0.0
    for i in range(n_quad):
        ui = float(us[i])
        wi = float(wts[i])
        Ti = ca.substitute(T_raw, u, ui)          # (3,1)
        nSi = ca.DM(nS_grid[i])                    # (3,)
        conj += wi * (ca.dot(Ti, nSi)) ** 2
        Tpi = ca.substitute(Tp_raw, u, ui)
        smooth += wi * ca.dot(Tpi, Tpi)
        norm_pen += wi * (ca.dot(Ti, Ti) - 1.0) ** 2

    obj = conj + w_smooth * smooth + 1.0 * norm_pen

    # 初值：法矢协方差最小特征向量（共轭最优方向，⊥ 平均法矢）
    C = np.zeros((3, 3))
    for ni in nS_grid:
        C += np.outer(ni, ni)
    _, V = np.linalg.eigh(C)
    T_init = V[:, 0]
    if np.dot(T_init, nS_grid.mean(axis=0)) > 0:
        T_init = -T_init
    # 列主序初值：pvec[i + j*n_ctrl] = T_init[j]（所有控制点同向 T_init）
    p0 = np.concatenate([np.full(n_ctrl, T_init[0]),
                         np.full(n_ctrl, T_init[1]),
                         np.full(n_ctrl, T_init[2])])

    nlp = {"x": pvec, "f": obj}
    opts = {"ipopt.print_level": 0, "print_time": 0}
    sol = ca.nlpsol("S", "ipopt", nlp, opts)(x0=p0)
    P_opt = np.array(sol["x"]).reshape(n_ctrl, 3, order="F")
    P_flat = P_opt.ravel(order="F")

    # 输出刀轴场 + 光顺度
    T_func = ca.Function("T", [u, pvec], [T_raw])
    Tp_func = ca.Function("Tp", [u, pvec], [Tp_raw])
    u_grid = np.linspace(u_min, u_max, 201)
    T_grid = np.array([T_func(ui, P_flat).full().ravel() for ui in u_grid])
    T_grid = T_grid / (np.linalg.norm(T_grid, axis=1, keepdims=True) + 1e-12)
    Tp_grid = np.array([Tp_func(ui, P_flat).full().ravel() for ui in u_grid])

    conj_vals = []
    for i in range(n_quad):
        nSi = nS_grid[i]
        Ti = np.array(T_func(float(us[i]), P_flat).full().ravel())
        Ti = Ti / (np.linalg.norm(Ti) + 1e-12)
        conj_vals.append(float(np.dot(Ti, nSi)) ** 2)

    # 刀位（CL）：刀心 A(u) = S(u,v_mid) + tool_r·nS；刀轴线段 A ± (L/2+tool_r)·T
    S_grid = np.array([surf(ui, v_mid)[0] for ui in u_grid])
    nS_dense = np.array([surf(ui, v_mid)[1] for ui in u_grid])
    A_grid = S_grid + tool_r * nS_dense
    # 母线长度（跨 v 两端距离的均值）
    edge_lens = []
    for ui in u_grid[::20]:
        edge_lens.append(np.linalg.norm(surf(ui, model["vEdges"][-1])[0] - surf(ui, model["vEdges"][0])[0]))
    L = float(np.mean(edge_lens)) if edge_lens else 0.0
    half = L * 0.5 + tool_r
    axis_segs = np.stack([A_grid - half * T_grid, A_grid + half * T_grid], axis=1)

    return {
        "u_grid": u_grid,
        "T_grid": T_grid,
        "Tp_grid": Tp_grid,
        "A_grid": A_grid,
        "axis_segs": axis_segs,
        "conjugate": {"mean": float(np.mean(conj_vals)), "rms": float(np.sqrt(np.mean(np.array(conj_vals) ** 2)))},
        "smoothness": {"mean": float(np.mean(np.linalg.norm(Tp_grid, axis=1))),
                       "max": float(np.max(np.linalg.norm(Tp_grid, axis=1)))},
        "meta": {
            "u_min": u_min, "u_max": u_max,
            "v_min": model["vEdges"][0], "v_max": model["vEdges"][-1],
            "fitDir": model["fitDir"], "continuity": k,
            "bandWidth": model["partition"]["bandWidth"],
            "n_cells": len(model["cells"]),
            "tool_r": tool_r,
        },
    }


def write_continuous_toolpath_vtk(res, out_path):
    """导出连续刀轨（进给折线 + 刀轴线段）到 VTK，供 UI 可视化。"""
    feed = res["A_grid"]
    axes = res["axis_segs"]
    lines = [feed]
    lines.extend([[seg[0], seg[1]] for seg in axes])
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


if __name__ == "__main__":
    import sys
    import os
    path = sys.argv[1] if len(sys.argv) > 1 else None
    if not path:
        print("usage: python toolpath_continuous.py <surface_model.json> [--out <dir>]")
        sys.exit(1)
    out_dir = None
    if "--out" in sys.argv:
        out_dir = sys.argv[sys.argv.index("--out") + 1]
    res = optimize_tool_axis_field(path)
    print(f"曲面: {res['meta']['n_cells']} 格, continuity=C{res['meta']['continuity']}, "
          f"bandWidth={res['meta']['bandWidth']}")
    print(f"共轭残差: mean={res['conjugate']['mean']:.6f} rms={res['conjugate']['rms']:.6f}")
    print(f"刀轴光顺度: mean={res['smoothness']['mean']:.6f} max={res['smoothness']['max']:.6f}")
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        stem = os.path.basename(path).replace("_surface_model.json", "")
        vtk_path = os.path.join(out_dir, f"{stem}_toolpath_continuous.vtk")
        write_continuous_toolpath_vtk(res, vtk_path)
        print(f"已导出连续刀轨: {vtk_path}")
