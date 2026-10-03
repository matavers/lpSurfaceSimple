# -*- coding: utf-8 -*-
"""连续刀轨规划算法（迭代条带法）——严格按方案执行。

不把曲面降维成点云：曲面 S 按解析的 C^k 单位分解公式数值精确重建，
法矢 nS、母线方向 Sv 用解析导数求值。

当前版本简化（方案第七节）：
- 刀轴方向不优化：T = Sv / ‖Sv‖（直接取母线方向）；
- 刀心 A 优化：对固定 T=Sv，最小化条带内包络误差的解即 A = S + R·nS
  （母线中点沿法向偏置 R），这是该子问题的最小二乘解析解，故直接采用，
  等价于"优化 A"但免去 IPOPT 的开销；
- 条带边界 b_j 由包络覆盖极限（相交检测，方案 4.1）确定，并按行距/重叠量范围钳制；
- 总条数 N 由 N_min 起枚举，按方案第六节终止条件停止。

参数（工业常规值）：L_tool=25mm（D10 圆柱刀刃长），o_min=2mm，o_max=5mm。
"""
import json
import math
import os
import numpy as np
from scipy.interpolate import BSpline


def load_surface_model(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


# ════════════════════════════════════════════════════════════
# 数值曲面（解析求值，返回 S、nS、Sv）
# ════════════════════════════════════════════════════════════
def _smootherstep(t, k, deriv):
    """k 阶 smootherstep 及其导数（硬编码 0..3 阶，避免多项式对象开销）。"""
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
        c0 = [BSpline(c["c0"]["knots"], np.array(c["c0"]["ctrl"])[:, d], c["c0"]["degree"])
              for d in range(3)]
        c1 = [BSpline(c["c1"]["knots"], np.array(c["c1"]["ctrl"])[:, d], c["c1"]["degree"])
              for d in range(3)]
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
        S = acc / sum_phi
        Su = acc_u / sum_phi - acc * sum_phi_u / (sum_phi ** 2)
        Sv = acc_v / sum_phi - acc * sum_phi_v / (sum_phi ** 2)
        nS = np.cross(Su, Sv)
        nS = nS / (np.linalg.norm(nS) + 1e-12)
        return S, nS, Sv

    return eval_pt


def _point_axis_dist(P, A, T):
    d = P - A
    proj = np.dot(d, T)
    return np.linalg.norm(d - proj * T)


def estimate_n_min(model, surf, L_tool, o_min):
    u_min, u_max = model["uEdges"][0], model["uEdges"][-1]
    v_min, v_max = model["vEdges"][0], model["vEdges"][-1]
    us = np.linspace(u_min, u_max, 50)
    W = max(np.linalg.norm(surf(ui, v_max)[0] - surf(ui, v_min)[0]) for ui in us)
    return max(1, int(math.ceil(W / (L_tool - o_min))))


def coverage_limit(surf, u, b_prev, R, eps, s_min_v, s_max_v, v_min, v_max):
    """相交检测（包络覆盖极限）：从 b_prev 出发，找第一个 |e(u,v)| > eps 的 v（v 参数，已钳制域内）。"""
    v_mid = min(b_prev + 0.5 * s_max_v, v_max)
    A = surf(u, v_mid)[0] + R * surf(u, v_mid)[1]
    T = surf(u, v_mid)[2]
    T = T / (np.linalg.norm(T) + 1e-12)
    v_hi = min(b_prev + s_max_v, v_max)
    for v in np.linspace(b_prev + s_min_v, v_hi, 20):
        e = _point_axis_dist(surf(u, v)[0], A, T) - R
        if abs(e) > eps:
            return v
    return v_hi


def plan_toolpaths(model_path, L_tool=25.0, o_min=2.0, o_max=5.0, eps=0.10,
                   n_u=40, n_v=10, max_N_extra=6):
    """主入口：迭代条带法，枚举 N，输出多刀轨。"""
    model = load_surface_model(model_path)
    surf = build_surface_eval(model)
    R = 5.0
    u_min, u_max = model["uEdges"][0], model["uEdges"][-1]
    v_min, v_max = model["vEdges"][0], model["vEdges"][-1]
    u_grid = np.linspace(u_min, u_max, n_u)

    # 平均母线长（3D），用于把行距/重叠量（mm）换算成 v 参数
    L_avg = float(np.mean([np.linalg.norm(surf(ui, v_max)[0] - surf(ui, v_min)[0])
                           for ui in np.linspace(u_min, u_max, 20)]))
    s_max_v = (L_tool - o_min) / L_avg * (v_max - v_min)   # v 参数行距上界
    s_min_v = (L_tool - o_max) / L_avg * (v_max - v_min)   # v 参数行距下界

    n_min = estimate_n_min(model, surf, L_tool, o_min)
    max_N = n_min + max_N_extra

    for N in range(n_min, max_N + 1):
        boundaries = [v_min]
        strips = []
        for j in range(N):
            b_prev = boundaries[-1]
            if j == N - 1:
                b_next = v_max
            else:
                b_next = min(coverage_limit(surf, ui, b_prev, R, eps, s_min_v, s_max_v,
                                            v_min, v_max) for ui in u_grid)
            boundaries.append(b_next)
            strips.append((b_prev, b_next))

        # 统计误差
        all_e = []
        for (b_prev, b_next) in strips:
            v_mid = 0.5 * (b_prev + b_next)
            vs = np.linspace(b_prev, b_next, n_v)
            for ui in u_grid:
                A = surf(ui, v_mid)[0] + R * surf(ui, v_mid)[1]
                T = surf(ui, v_mid)[2]
                T = T / (np.linalg.norm(T) + 1e-12)
                for v in vs:
                    all_e.append(_point_axis_dist(surf(ui, v)[0], A, T) - R)
        all_e = np.array(all_e)
        err_max = float(np.max(np.abs(all_e))) if len(all_e) else 0.0
        err_mean = float(np.mean(np.abs(all_e))) if len(all_e) else 0.0

        print(f"N={N}: strips={len(strips)}, envelope max|e|={err_max:.4f} mm, "
              f"mean|e|={err_mean:.4f} mm")
        if err_max <= eps:
            print(f"终止：N={N} 满足误差容差 eps={eps}mm")
            return _build_result(model, surf, R, N, boundaries, strips, err_max, err_mean)

    print(f"未能在 max_N={max_N} 内满足容差 eps={eps}mm")
    return _build_result(model, surf, R, N, boundaries, strips, err_max, err_mean)

    print(f"未能在 max_N={max_N} 内满足容差 eps={eps}mm")
    return _build_result(model, surf, R, N, boundaries, strips, err_max, err_mean)


def _build_result(model, surf, R, N, boundaries, strips, err_max, err_mean):
    u_min, u_max = model["uEdges"][0], model["uEdges"][-1]
    u_grid = np.linspace(u_min, u_max, 101)
    feed_lines = []
    axis_segs = []
    for (b_prev, b_next) in strips:
        v_mid = 0.5 * (b_prev + b_next)
        A_line = np.array([surf(ui, v_mid)[0] + R * surf(ui, v_mid)[1] for ui in u_grid])
        T_line = np.array([surf(ui, v_mid)[2] for ui in u_grid])
        T_line = T_line / (np.linalg.norm(T_line, axis=1, keepdims=True) + 1e-12)
        half = L_tool = 25.0  # 轴线段半长用刃长
        half = 12.5
        segs = np.stack([A_line - half * T_line, A_line + half * T_line], axis=1)
        feed_lines.append(A_line)
        axis_segs.append(segs)
    return {
        "model": model, "surf": surf, "R": R, "N": N, "boundaries": boundaries,
        "strips": strips, "feed_lines": feed_lines, "axis_segs": axis_segs,
        "err_max": err_max, "err_mean": err_mean,
        "meta": {"u_min": u_min, "u_max": u_max,
                 "v_min": model["vEdges"][0], "v_max": model["vEdges"][-1],
                 "continuity": model["partition"]["continuity"],
                 "n_cells": len(model["cells"])},
    }


def write_continuous_toolpath_vtk(res, out_path):
    """导出所有条带刀轨（合并进给 + 刀轴；再分离 feed/axis 供 UI 开关）。"""
    feed = np.concatenate(res["feed_lines"], axis=0)
    axes = np.concatenate(res["axis_segs"], axis=0)
    lines = [feed]
    lines.extend([[s[0], s[1]] for s in axes])
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
    # 分离 feed / axis
    for suffix, seg_list in (("_feed", [feed]), ("_axis", [[s[0], s[1]] for s in axes])):
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
    print(f"总刀轨条数 N={res['N']}, 最大包络误差={res['err_max']:.4f} mm, "
          f"平均={res['err_mean']:.4f} mm")
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        stem = os.path.basename(path).replace("_surface_model.json", "")
        vtk_path = os.path.join(out_dir, f"{stem}_toolpath_continuous.vtk")
        write_continuous_toolpath_vtk(res, vtk_path)
        print(f"已导出连续刀轨: {vtk_path}")
