# -*- coding: utf-8 -*-
"""在完整 C^k 连续曲面上做可微刀位规划（CasADi + IPOPT）——完整版。

不把曲面降维成点云：曲面 S 按解析的 C^k 单位分解公式用 CasADi 符号重建，
法矢 nS 用自动微分解析求得；刀轴场 T(u)、刀心 A(u)、接触点 v(u) 都用 B 样条
参数化（控制点可学习）。目标泛函
  J = ∫[ (ρ(S(u,v(u)), axis(A,T)) − R)^2 + w1·(nS·T)^2
         + w2·‖T'‖^2 + w3·‖A'‖^2 + w4·‖v'‖^2 ] du
在 Gauss-Legendre 求积点上求值，用 IPOPT 梯度下降。
连续性体现在：S∈C^k ⇒ nS∈C^{k-1} 连续 ⇒ 目标泛函连续可微 ⇒ 精确梯度。
"""
import json
import math
import os
import numpy as np
from scipy.interpolate import BSpline
import casadi as ca


def load_surface_model(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


# ════════════════════════════════════════════════════════════
# 数值曲面（用于刀位输出 VTK + 初值；非优化路径）
# ════════════════════════════════════════════════════════════
def _smootherstep_poly(k):
    if k <= 0:
        return np.array([0.0, 1.0]), np.array([1.0])
    p = np.polynomial.Polynomial([0.0])
    for i in range(k + 1):
        omt = np.polynomial.Polynomial([1.0, -1.0]) ** i
        p = p + math.comb(k + i, i) * omt
    p = p * np.polynomial.Polynomial([0.0] * (k + 1) + [1.0])
    c = p.coef
    dc = np.polynomial.Polynomial(c).deriv().coef
    return c, dc


def _chi1_np(u, uL, uR, a, b, k, deriv=False):
    c, dc = _smootherstep_poly(k)

    def val(coef, x):
        return np.polynomial.Polynomial(coef)(x)

    if u <= uL - a:
        return (0.0, 0.0) if deriv else 0.0
    if u < uL:
        t = (u - (uL - a)) / a
        return (val(c, t), val(dc, t) / a) if deriv else val(c, t)
    if u <= uR:
        return (1.0, 0.0) if deriv else 1.0
    if u < uR + b:
        t = (u - uR) / b
        return (1.0 - val(c, t), -val(dc, t) / b) if deriv else 1.0 - val(c, t)
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
                c0v = float(c0[d](u)); c1v = float(c1[d](u))
                c0p = float(c0[d].derivative()(u)); c1p = float(c1[d].derivative()(u))
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
        nS = nS / (np.linalg.norm(nS) + 1e-12)
        return S, nS

    return eval_pt


# ════════════════════════════════════════════════════════════
# Cox-de Boor B 样条（CasADi MX，可微；规避 ca.bspline 的 AD 缺陷）
# ════════════════════════════════════════════════════════════
def _bspline_mx(u, P, knots, degree):
    """B 样条 C(u)=Σ N_{i,p}(u) P_i。P 为 (n_ctrl, dim) MX，返回 (dim,1)。"""
    from functools import lru_cache
    n = len(knots) - degree - 1
    dim = P.shape[1]

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
    return ca.vertcat(*[sum(N[i] * P[i, d] for i in range(n)) for d in range(dim)])


def _bspline_deriv_ctrl(P, knots, degree):
    """B 样条导数控制点（线性于 P），返回 (n_ctrl-1, dim)。"""
    n = P.shape[0]
    rows = []
    for j in range(n - 1):
        denom = knots[j + degree + 1] - knots[j + 1]
        w = degree / denom if denom > 1e-12 else 0.0
        rows.append((P[j + 1, :] - P[j, :]) * w)
    return ca.vertcat(*rows)


# ════════════════════════════════════════════════════════════
# 符号曲面（CasADi MX，对 u、v 都可微）
# ════════════════════════════════════════════════════════════
def _smootherstep_mx(t, k):
    if k <= 0:
        return t
    omt = 1.0 - t
    s = 0
    for i in range(k + 1):
        s = s + math.comb(k + i, i) * omt ** i
    return t ** (k + 1) * s


def _chi1_mx(u, uL, uR, a, b, k):
    tL = (u - (uL - a)) / a
    tR = (u - uR) / b
    return ca.if_else(
        u <= uL - a, 0.0,
        ca.if_else(u < uL, _smootherstep_mx(tL, k),
        ca.if_else(u <= uR, 1.0,
        ca.if_else(u < uR + b, 1.0 - _smootherstep_mx(tR, k), 0.0))))


def build_surface_mx(model):
    """符号重建 C^k 曲面 S(u,v) 与法矢 nS(u,v)。返回 (S_func, nS_func, meta)。"""
    k = model["partition"]["continuity"]
    band_w = model["partition"]["bandWidth"]
    cells = model["cells"]
    min_u = min(c["u1"] - c["u0"] for c in cells)
    min_v = min(c["v1"] - c["v0"] for c in cells)
    a_u = band_w * min_u
    a_v = band_w * min_v

    u = ca.MX.sym("u")
    v = ca.MX.sym("v")
    acc = ca.MX.zeros(3, 1)
    sum_phi = 0.0
    for c in cells:
        cu = _chi1_mx(u, c["u0"], c["u1"], a_u, a_u, k)
        cv = _chi1_mx(v, c["v0"], c["v1"], a_v, a_v, k)
        phi = cu * cv
        c0 = _bspline_mx(u, ca.DM(np.array(c["c0"]["ctrl"], dtype=float)), c["c0"]["knots"], c["c0"]["degree"])
        c1 = _bspline_mx(u, ca.DM(np.array(c["c1"]["ctrl"], dtype=float)), c["c1"]["knots"], c["c1"]["degree"])
        R = (1.0 - v) * c0 + v * c1
        acc = acc + phi * R
        sum_phi = sum_phi + phi
    S = acc / sum_phi

    Su = ca.jacobian(S, u)
    Sv = ca.jacobian(S, v)
    nS = ca.cross(Su, Sv)
    nS = nS / (ca.norm_2(nS) + 1e-12)

    S_func = ca.Function("S", [u, v], [S])
    nS_func = ca.Function("nS", [u, v], [nS])
    meta = {
        "u_min": model["uEdges"][0], "u_max": model["uEdges"][-1],
        "v_min": model["vEdges"][0], "v_max": model["vEdges"][-1],
        "fitDir": model["fitDir"], "continuity": k,
        "bandWidth": band_w, "n_cells": len(cells),
    }
    return S_func, nS_func, meta


# ════════════════════════════════════════════════════════════
# 完整版刀轴场优化（CasADi + IPOPT）
# ════════════════════════════════════════════════════════════
def optimize_tool_axis_field(model_path, n_ctrl=16, n_quad=40, tool_r=5.0,
                             w_conj=1.0, w_smooth_T=0.5, w_smooth_A=0.05,
                             w_smooth_v=0.05, continuity=None, max_iter=300):
    model = load_surface_model(model_path)
    k = continuity if continuity is not None else model["partition"]["continuity"]
    u_min, u_max = model["uEdges"][0], model["uEdges"][-1]
    v_min, v_max = model["vEdges"][0], model["vEdges"][-1]
    v_mid = 0.5 * (v_min + v_max)

    # Gauss-Legendre 求积点
    pts, wts = np.polynomial.legendre.leggauss(n_quad)
    us = 0.5 * (u_max - u_min) * pts + 0.5 * (u_max + u_min)
    wts = wts * 0.5 * (u_max - u_min)

    # 数值预计算曲面在求积点 (u_g, v_mid) 的 S、nS、Sv（解析求值，非点云）
    # v 向用线性 Taylor 模型 P(v)=S0+Sv0·(v-v_mid)，接触点 v 作为变量仍可微。
    surf = build_surface_eval(model)
    eps = 1e-4 * (v_max - v_min)
    S0 = np.array([surf(ui, v_mid)[0] for ui in us])
    nS0 = np.array([surf(ui, v_mid)[1] for ui in us])
    Sv0 = np.array([(surf(ui, v_mid + eps)[0] - surf(ui, v_mid - eps)[0]) / (2 * eps) for ui in us])

    # 刀轴/刀心/接触点 B 样条（degree = k+1，C^k 连续）
    degree = k + 1
    n_ctrl = max(degree + 1, n_ctrl)
    inner = n_ctrl - degree - 1
    knots = [u_min] * (degree + 1) + \
            [u_min + (u_max - u_min) * (i + 1) / (inner + 1) for i in range(inner)] + \
            [u_max] * (degree + 1)

    u = ca.MX.sym("u")
    nT = 3 * n_ctrl
    nA = 3 * n_ctrl
    nv = n_ctrl
    xvec = ca.MX.sym("x", nT + nA + nv)
    pT = ca.reshape(xvec[0:nT], n_ctrl, 3)
    pA = ca.reshape(xvec[nT:nT + nA], n_ctrl, 3)
    pv = ca.reshape(xvec[nT + nA:], n_ctrl, 1)

    T_raw = _bspline_mx(u, pT, knots, degree)          # (3,1)
    A_raw = _bspline_mx(u, pA, knots, degree)          # (3,1)
    v_raw = _bspline_mx(u, pv, knots, degree)          # (1,1)
    v_sig = v_min + (v_max - v_min) / (1.0 + ca.exp(-v_raw))  # 映射到 [v_min,v_max]

    Tp_raw = _bspline_mx(u, _bspline_deriv_ctrl(pT, knots, degree), knots[1:-1], degree - 1)
    Ap_raw = _bspline_mx(u, _bspline_deriv_ctrl(pA, knots, degree), knots[1:-1], degree - 1)
    vp_raw = _bspline_mx(u, _bspline_deriv_ctrl(pv, knots, degree), knots[1:-1], degree - 1)

    obj = 0.0
    for i in range(n_quad):
        ui = float(us[i])
        wi = float(wts[i])
        Ti = ca.substitute(T_raw, u, ui)               # (3,1)
        Ai = ca.substitute(A_raw, u, ui)
        vi = ca.substitute(v_sig, u, ui)               # (1,1)
        Pi = ca.DM(S0[i]) + ca.DM(Sv0[i]) * (vi - v_mid)   # 线性 Taylor 接触点 (3,1)
        nSi = ca.DM(nS0[i])                            # (3,1)
        d = Pi - Ai
        proj = ca.dot(d, Ti)
        rho = ca.norm_2(d - Ti * proj)
        e = rho - tool_r
        conj = ca.dot(nSi, Ti)
        Tpi = ca.substitute(Tp_raw, u, ui)
        Api = ca.substitute(Ap_raw, u, ui)
        vpi = ca.substitute(vp_raw, u, ui)
        obj += wi * (e ** 2 + w_conj * conj ** 2
                     + w_smooth_T * ca.dot(Tpi, Tpi)
                     + w_smooth_A * ca.dot(Api, Api)
                     + w_smooth_v * ca.dot(vpi, vpi)
                     + 1.0 * (ca.dot(Ti, Ti) - 1.0) ** 2)

    # 初值
    nS_init = nS0
    C = np.zeros((3, 3))
    for ni in nS_init:
        C += np.outer(ni, ni)
    _, V = np.linalg.eigh(C)
    T_init = V[:, 0]
    if np.dot(T_init, nS_init.mean(axis=0)) > 0:
        T_init = -T_init

    p0_T = np.tile(T_init, (n_ctrl, 1))                       # (n_ctrl,3)
    p0_A = np.array([surf(ui, v_mid)[0] + tool_r * surf(ui, v_mid)[1]
                     for ui in np.linspace(u_min, u_max, n_ctrl)])  # (n_ctrl,3)
    p0_v = np.full((n_ctrl, 1), 0.0)                          # v_raw=0 -> v=v_mid

    # 决策变量拼接（列主序展平）
    x0 = np.concatenate([p0_T.ravel(order="F"), p0_A.ravel(order="F"), p0_v.ravel()])

    nlp = {"x": xvec, "f": obj}
    opts = {"ipopt.print_level": 0, "print_time": 0, "ipopt.max_iter": max_iter}
    print("[optimize] IPOPT solving...")
    sol = ca.nlpsol("S", "ipopt", nlp, opts)(x0=x0)
    x_opt = np.array(sol["x"]).ravel()
    T_opt = x_opt[0:nT].reshape(n_ctrl, 3, order="F")
    A_opt = x_opt[nT:nT + nA].reshape(n_ctrl, 3, order="F")
    v_opt = x_opt[nT + nA:].reshape(n_ctrl, 1)

    # 输出刀轴场
    T_func = ca.Function("T", [u, xvec], [T_raw])
    A_func = ca.Function("A", [u, xvec], [A_raw])
    v_func = ca.Function("v", [u, xvec], [v_sig])
    Tp_func = ca.Function("Tp", [u, xvec], [Tp_raw])
    u_grid = np.linspace(u_min, u_max, 201)
    T_grid = np.array([T_func(ui, x_opt).full().ravel() for ui in u_grid])
    T_grid = T_grid / (np.linalg.norm(T_grid, axis=1, keepdims=True) + 1e-12)
    A_grid = np.array([A_func(ui, x_opt).full().ravel() for ui in u_grid])
    v_grid = np.array([v_func(ui, x_opt).full().ravel() for ui in u_grid])
    Tp_grid = np.array([Tp_func(ui, x_opt).full().ravel() for ui in u_grid])

    # 点轴残差与共轭统计（数值求值）
    resid = []
    conjv = []
    for i in range(n_quad):
        ui = float(us[i])
        Ti = T_func(ui, x_opt).full().ravel()
        Ti = Ti / (np.linalg.norm(Ti) + 1e-12)
        Ai = A_func(ui, x_opt).full().ravel()
        vi = float(v_func(ui, x_opt))
        Pi = surf(ui, vi)[0]
        nSi = surf(ui, vi)[1]
        d = Pi - Ai
        rho = np.linalg.norm(d - np.dot(d, Ti) * Ti)
        resid.append(abs(rho - tool_r))
        conjv.append(float(np.dot(nSi, Ti)) ** 2)

    # 刀轴线段（母线长度近似）
    edge_lens = [np.linalg.norm(surf(ui, v_max)[0] - surf(ui, v_min)[0]) for ui in u_grid[::20]]
    L = float(np.mean(edge_lens)) if edge_lens else 0.0
    half = L * 0.5 + tool_r
    axis_segs = np.stack([A_grid - half * T_grid, A_grid + half * T_grid], axis=1)

    return {
        "u_grid": u_grid, "T_grid": T_grid, "Tp_grid": Tp_grid,
        "A_grid": A_grid, "v_grid": v_grid, "axis_segs": axis_segs,
        "residual": {"mean": float(np.mean(resid)), "rms": float(np.sqrt(np.mean(np.array(resid) ** 2))),
                     "max": float(np.max(resid))},
        "conjugate": {"mean": float(np.mean(conjv)), "rms": float(np.sqrt(np.mean(np.array(conjv) ** 2)))},
        "smoothness": {"mean": float(np.mean(np.linalg.norm(Tp_grid, axis=1))),
                       "max": float(np.max(np.linalg.norm(Tp_grid, axis=1)))},
        "meta": {
            "u_min": u_min, "u_max": u_max, "v_min": v_min, "v_max": v_max,
            "fitDir": model["fitDir"], "continuity": k,
            "bandWidth": model["partition"]["bandWidth"],
            "n_cells": len(model["cells"]), "tool_r": tool_r,
        },
    }


def write_continuous_toolpath_vtk(res, out_path):
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
    print(f"点轴残差(mm): max={res['residual']['max']:.4f} mean={res['residual']['mean']:.4f} "
          f"rms={res['residual']['rms']:.4f}")
    print(f"共轭残差: mean={res['conjugate']['mean']:.6f}")
    print(f"刀轴光顺度: mean={res['smoothness']['mean']:.6f} max={res['smoothness']['max']:.6f}")
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        stem = os.path.basename(path).replace("_surface_model.json", "")
        vtk_path = os.path.join(out_dir, f"{stem}_toolpath_continuous.vtk")
        write_continuous_toolpath_vtk(res, vtk_path)
        print(f"已导出连续刀轨: {vtk_path}")
