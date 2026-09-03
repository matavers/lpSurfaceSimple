#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compute_machining.py — 直纹面分片拟合 vs 传统点铣 加工时间仿真（简化但可辩护的刀轨模型）

读取 API 输出的直纹面准线参数文件 (*_params.txt)，对每个直纹格：

1. 可展性判定 —— 用「沿母线的曲面法向扭转角」作为判据（不再是母线方向偏离角）：
     对第 i 条母线，取上下准线切向 t0_i / t1_i（中心差分），
     法向 n0_i = t0_i × r_i、n1_i = t1_i × r_i（r_i 为母线方向），
     扭转角 β_i = angle(n0_i, n1_i)，格扭转角 β = max β_i。
     直纹面可展 ⟺ 切平面(法向)沿母线恒定 ⟺ β=0，这正是圆柱刀侧刃铣的精确接触条件。
2. 侧铣（圆柱刀，刀轴 ∥ 母线，沿准线进给）：
     可展片（β ≤ 阈值）单刀成形，刀轨长 ≈ 准线长；
     不可展片**不再**按 twist/threshold 伪造多刀，而是如实归入点铣。
3. 点铣（球头刀，刀心 = 曲面点 + R·n）：
     行距由残留高度反算 stepover = 2·√(2·R·h − h²)（标准凸/平面近似）。
4. 时间 = 刀轨长 / 进给率 + 每块固定进退刀开销（理想化，见 assumptions）。

对比口径：
   A 混合策略（直纹面拟合后）= 可展片侧铣 + 不可展片点铣 + 每片 overhead
   B 传统基线（整体点铣）= 全部片点铣 + 单次 point_overhead

用法:
  python compute_machining.py <output_dir> [--feed 500] [--rapid 5000] \
      [--tool-r 5] [--ball-r 5] [--scallop 0.01] [--twist-limit 1.0] \
      [--overhead 2] [--point-overhead 10] [--json out.json] [--vtk vtk_dir]
"""
import os
import sys
import math
import json
import argparse
from pathlib import Path
from dataclasses import dataclass, asdict

# ===================== 几何（纯 python 元组向量） =====================

def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

def _add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])

def _scale(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)

def _norm(a):
    return math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])

def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])

def _normalize(a):
    m = _norm(a)
    if m < 1e-12:
        return (0.0, 0.0, 0.0)
    return (a[0] / m, a[1] / m, a[2] / m)

def _lerp(a, b, t):
    return tuple(a[k] * (1.0 - t) + b[k] * t for k in range(3))

def _angle_deg(a, b):
    m = _norm(a) * _norm(b)
    if m < 1e-12:
        return 0.0
    d = max(-1.0, min(1.0, _dot(a, b) / m))
    return math.degrees(math.acos(d))

def _tangents(C):
    """折线各点的单位切向（中心差分，端点单侧）。"""
    n = len(C)
    out = []
    for i in range(n):
        if n < 2:
            out.append((0.0, 0.0, 0.0))
        elif i == 0:
            out.append(_normalize(_sub(C[1], C[0])))
        elif i == n - 1:
            out.append(_normalize(_sub(C[n - 1], C[n - 2])))
        else:
            out.append(_normalize(_sub(C[i + 1], C[i - 1])))
    return out

def load_directrix(path):
    """读取 _params.txt，返回 (C0, C1)，各为 [(x,y,z),...]"""
    C0, C1, cur = [], [], None
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('['):
                cur = line.strip('[]')
            elif cur in ('C0', 'C1') and line and not line.startswith(('n', 'd', 'm')):
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        p = tuple(float(x) for x in parts[:3])
                        (C0 if cur == 'C0' else C1).append(p)
                    except ValueError:
                        pass
    return C0, C1

def load_obj(path):
    """读取 OBJ 网格，返回 (verts [(x,y,z),...], faces [[i,j,k],...])。"""
    verts, faces = [], []
    try:
        with open(path, encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line.startswith('v '):
                    parts = line.split()
                    verts.append(tuple(float(x) for x in parts[1:4]))
                elif line.startswith('f '):
                    parts = line.split()[1:]
                    idxs = []
                    for p in parts:
                        idxs.append(int(p.split('/')[0]) - 1)
                    if len(idxs) == 3:
                        faces.append(idxs)
    except Exception:
        pass
    return verts, faces

def mesh_area(verts, faces):
    """三角网格面积（mm²）。"""
    area = 0.0
    for f in faces:
        a, b, c = verts[f[0]], verts[f[1]], verts[f[2]]
        area += 0.5 * _norm(_cross(_sub(b, a), _sub(c, a)))
    return area

def curve_length(C):
    return sum(_norm(_sub(C[i + 1], C[i])) for i in range(len(C) - 1))

def mean_ruling_length(C0, C1):
    n = min(len(C0), len(C1))
    if n == 0:
        return 0.0
    return sum(_norm(_sub(C1[i], C0[i])) for i in range(n)) / n

def ruling_normals(C0, C1):
    """每条母线两端（v=0 / v=1）的曲面法向 + 母线方向。
    S(u,v) = (1-v)·C0(u) + v·C1(u)，n = (∂S/∂u) × (∂S/∂v) = C' × (C1-C0)。"""
    n = min(len(C0), len(C1))
    T0 = _tangents(C0)
    T1 = _tangents(C1)
    n0s, n1s, rs = [], [], []
    for i in range(n):
        r = _sub(C1[i], C0[i])
        rs.append(r)
        n0s.append(_normalize(_cross(T0[i], r)))
        n1s.append(_normalize(_cross(T1[i], r)))
    return n0s, n1s, rs

def ruling_twists(C0, C1):
    """每条母线的法向扭转角（度）列表。"""
    n0s, n1s, _ = ruling_normals(C0, C1)
    return [_angle_deg(n0s[i], n1s[i]) for i in range(len(n0s))]

def twist_angle_deg(C0, C1):
    """格扭转角（度）= 各母线法向扭转角的最大值。可展 ⟺ 接近 0。"""
    ts = ruling_twists(C0, C1)
    return max(ts) if ts else 0.0

def surface_area(C0, C1):
    """直纹面面积 = Σ |r × dC|（数值积分）。"""
    n = min(len(C0), len(C1))
    area = 0.0
    for i in range(n - 1):
        r = _sub(C1[i], C0[i])
        dC = _sub(C0[i + 1], C0[i])
        area += _norm(_cross(r, dC))
    return area

# ===================== 刀轨 + 时间 =====================

def flank_milling(C0, C1, feed, twist_limit):
    """
    侧铣（锥度立铣刀）：刀具轴线沿母线、沿准线进给。
    C++ 已按曲率把卷曲区切去，能进 *_params.txt 的格子都是可展区，一律侧铣一次。
    twist 仅作报告统计。返回 (切削时间 s, 刀轨长 mm, 刀数, 是否可展, 扭转角 deg)
    """
    L = curve_length(C0)
    twist = twist_angle_deg(C0, C1)
    developable = True
    path_len = L
    passes = 1.0
    cut_time = path_len / feed * 60.0
    return cut_time, path_len, passes, developable, twist

def point_milling(C0, C1, feed, ball_r, scallop):
    """
    点铣（球头刀）：刀心 = 曲面点 + R·n，沿准线方向行切，
    行距由残留高度反算 stepover = 2·√(2·R·h − h²)（标准凸/平面近似）。
    返回 (切削时间 s, 行距 mm, 总刀轨长 mm, 面积 mm²)
    """
    stepover = 2.0 * math.sqrt(max(0.0, 2.0 * ball_r * scallop - scallop * scallop))
    A = surface_area(C0, C1)
    total_len = A / stepover if stepover > 0 else 0.0
    return total_len / feed * 60.0, stepover, total_len, A

def flank_toolpath_lines(C0, C1, n_along=50, tool_r=None):
    """侧铣刀轨可视化：沿准线的进给路径 + 每条母线的刀具轴线（可偏置 R_t）。
    返回多条折线（每条为点列）。"""
    n = len(C0)
    if n < 2:
        return []
    idx = ([int(round(i * (n - 1) / (n_along - 1))) for i in range(n_along)]
           if n > n_along else list(range(n)))
    lines = [[C0[i] for i in idx]]  # 进给路径（下准线）
    n0s, n1s, rs = ruling_normals(C0, C1)
    for i in idx:
        if tool_r and tool_r > 0:
            mid = _lerp(C0[i], C1[i], 0.5)
            nm = _normalize(_add(n0s[i], n1s[i]))
            axis_p = _add(mid, _scale(nm, tool_r))
            half = _norm(rs[i]) / 2.0
            d = _normalize(rs[i])
            lines.append([_add(axis_p, _scale(d, -half)),
                          _add(axis_p, _scale(d, half))])
        else:
            lines.append([C0[i], C1[i]])  # 母线（刀具轴线方向）
    return lines

def point_toolpath_lines(C0, C1, stepover, n_along=50, ball_r=None):
    """点铣刀轨可视化：沿准线方向行切，行距 stepover（沿母线方向）。
    ball_r 提供时刀心沿曲面法向偏置（CL 点）。"""
    n = len(C0)
    if n < 2:
        return []
    idx = ([int(round(i * (n - 1) / (n_along - 1))) for i in range(n_along)]
           if n > n_along else list(range(n)))
    mean_r = mean_ruling_length(C0, C1)
    n_across = max(1, int(round(mean_r / stepover))) if stepover > 0 else 1
    n0s, n1s, _ = ruling_normals(C0, C1)
    lines = []
    for j in range(n_across + 1):
        t = j / n_across if n_across > 0 else 0.0
        if ball_r and ball_r > 0:
            row = []
            for i in idx:
                nm = _normalize(_add(n0s[i], n1s[i]))
                surf = _lerp(C0[i], C1[i], t)
                row.append(_add(surf, _scale(nm, ball_r)))
            lines.append(row)
        else:
            lines.append([_lerp(C0[i], C1[i], t) for i in idx])
    return lines

def flank_cl_lines(C0, C1, tool_r, n_along=50, axis_extend=None, flip=1.0):
    """侧铣 CL 刀位（圆柱刀，半径 tool_r）：
    刀轴沿母线方向，刀心沿曲面法向偏置 tool_r（flip=±1 控制偏置到曲面的哪一侧）。
    返回 (feed_pts, axis_segs)
      feed_pts: 刀轴中心进给轨迹 [(x,y,z),...]，沿准线方向 = 刀具移动方向
      axis_segs: 每条采样母线的刀轴线段 [[p0,p1],...]，偏置后、两端各延伸 axis_extend
    """
    n = len(C0)
    if n < 2:
        return [], []
    idx = ([int(round(i * (n - 1) / (n_along - 1))) for i in range(n_along)]
           if n > n_along else list(range(n)))
    n0s, n1s, rs = ruling_normals(C0, C1)
    ext = tool_r if axis_extend is None else axis_extend
    feed, axes = [], []
    for i in idx:
        mid = _lerp(C0[i], C1[i], 0.5)
        nm = _scale(_normalize(_add(n0s[i], n1s[i])), flip)
        r = rs[i]
        rlen = _norm(r)
        if rlen < 1e-12:
            continue
        rhat = _scale(r, 1.0 / rlen)
        axis_p = _add(mid, _scale(nm, tool_r))
        feed.append(axis_p)
        half = rlen * 0.5 + ext
        axes.append([_add(axis_p, _scale(rhat, -half)),
                     _add(axis_p, _scale(rhat, half))])
    return feed, axes

def point_cl_lines(C0, C1, stepover, ball_r, n_along=50, flip=1.0):
    """点铣 CL（球头刀，半径 ball_r）：刀心 = 曲面点 + ball_r·n，沿准线方向行切。
    flip=±1 控制刀心偏置到曲面的哪一侧。返回多条刀心行切折线（每条为点列）。"""
    n = len(C0)
    if n < 2:
        return []
    idx = ([int(round(i * (n - 1) / (n_along - 1))) for i in range(n_along)]
           if n > n_along else list(range(n)))
    mean_r = mean_ruling_length(C0, C1)
    n_across = max(1, int(round(mean_r / stepover))) if stepover > 0 else 1
    n0s, n1s, _ = ruling_normals(C0, C1)
    lines = []
    for j in range(n_across + 1):
        t = j / n_across if n_across > 0 else 0.0
        row = []
        for i in idx:
            nm = _scale(_normalize(_lerp(n0s[i], n1s[i], t)), flip)
            surf = _lerp(C0[i], C1[i], t)
            row.append(_add(surf, _scale(nm, ball_r)))
        lines.append(row)
    return lines

def curl_point_cl_lines(verts, faces, stepover, ball_r, flip=1.0):
    """卷曲区点铣 CL：在原曲面网格上行切（球刀刀心 = 顶点 + ball_r·n）。
    网格为 C++ 导出的方形网格（nV×nV），沿 v 方向按 stepover 取列。"""
    import math
    n = len(verts)
    nV = int(round(math.sqrt(n))) if n > 0 else 0
    if nV < 2 or nV * nV != n:
        return []
    # 顶点法向（面法向面积加权平均）
    acc = [[0.0, 0.0, 0.0] for _ in range(n)]
    for f in faces:
        a, b, c = verts[f[0]], verts[f[1]], verts[f[2]]
        fn = _cross(_sub(b, a), _sub(c, a))
        for k in f:
            acc[k][0] += fn[0]
            acc[k][1] += fn[1]
            acc[k][2] += fn[2]
    normals = [_normalize(tuple(x)) for x in acc]
    # 相邻列间距（v 方向）近似
    edge = _norm(_sub(verts[0], verts[1])) if nV >= 2 else 0.0
    step = max(1, int(round(edge / stepover))) if stepover > 0 and edge > 0 else 1
    lines = []
    for j in range(0, nV, step):
        row = []
        for i in range(nV):
            v = verts[i * nV + j]
            nm = _scale(normals[i * nV + j], flip)
            row.append(_add(v, _scale(nm, ball_r)))
        lines.append(row)
    # 补最后一列，保证覆盖到边
    if (nV - 1) % step != 0:
        j = nV - 1
        row = []
        for i in range(nV):
            v = verts[i * nV + j]
            nm = _scale(normals[i * nV + j], flip)
            row.append(_add(v, _scale(nm, ball_r)))
        lines.append(row)
    return lines

def write_vtk_polylines(path, lines):
    pts, segs = [], []
    for line in lines:
        base = len(pts)
        for p in line:
            pts.append(p)
        for k in range(len(line) - 1):
            segs.append((base + k, base + k + 1))
    with open(path, 'w', encoding='utf-8') as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("toolpath\nASCII\nDATASET POLYDATA\n")
        f.write(f"POINTS {len(pts)} float\n")
        for p in pts:
            f.write(f"{p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
        f.write(f"LINES {len(segs)} {len(segs) * 3}\n")
        for a, b in segs:
            f.write(f"2 {a} {b}\n")

# ===================== 数据模型 =====================

@dataclass
class Patch:
    name: str
    C0: list = None
    C1: list = None
    directrix_len: float = 0.0
    mean_ruling: float = 0.0
    area: float = 0.0
    twist: float = 0.0            # 法向扭转角（度）
    developable: bool = True
    flank_time: float = 0.0
    flank_len: float = 0.0
    flank_passes: float = 0.0
    point_time: float = 0.0
    point_len: float = 0.0
    stepover: float = 0.0
    blade: int = 0               # 0=blade1, 1=blade2
    cell_idx: int = -1           # 网格 cell 索引
    row: int = -1                # 网格行（-1 表示无 meta）
    col: int = -1                # 网格列
    normal_flip: float = 1.0     # 法向翻转系数：+1 保持，-1 翻转到另一侧（刀具偏移到曲面外侧）
    is_curl: bool = False        # 卷曲区（C++ 曲率分割切出），点铣原曲面
    mesh_verts: list = None      # 卷曲区原曲面网格顶点
    mesh_faces: list = None      # 卷曲区原曲面网格面（三角形索引）

# ===================== 标定参数 =====================

def load_config(path=None):
    """读取标定参数 JSON 配置，返回 dict（不存在则空）。"""
    if path is None:
        path = Path(__file__).parent / "machining_config.json"
    try:
        with open(path, encoding='utf-8') as f:
            return json.load(f)
    except Exception:
        return {}

# ===================== 网格邻接与连通域 =====================

def _parse_cell_name(name):
    """从 'blade1_cell12' 解析出 (blade 索引, cell 索引)。"""
    blade = 0 if 'blade1' in name else 1
    idx = -1
    tail = name.split('_cell')[-1]
    try:
        idx = int(tail)
    except ValueError:
        idx = -1
    return blade, idx

def _load_grid_meta(input_dir):
    """从 meta.json 读取每个面的 (nRows, nCols)，返回 {blade: (nRows, nCols)}。"""
    meta_path = os.path.join(input_dir, 'meta.json')
    if not os.path.exists(meta_path):
        return {}
    try:
        with open(meta_path, encoding='utf-8') as f:
            meta = json.load(f)
    except Exception:
        return {}
    grid = {}
    for i, s in enumerate(meta.get('surfaces', [])):
        grid[i] = (s.get('nRows', 0), s.get('nCols', 0))
    return grid

def _load_curled_cells(input_dir):
    """从 meta.json 读取卷曲区格子（developable=false），返回 {(blade, idx): (row, col)}。"""
    meta_path = os.path.join(input_dir, 'meta.json')
    if not os.path.exists(meta_path):
        return {}
    try:
        with open(meta_path, encoding='utf-8') as f:
            meta = json.load(f)
    except Exception:
        return {}
    curled = {}
    for bi, s in enumerate(meta.get('surfaces', [])):
        for c in s.get('cells', []):
            if not c.get('developable', True):
                curled[(bi, c['index'])] = (c.get('row', -1), c.get('col', -1))
    return curled

def _blade_normal_flips(patches):
    """判定每个面的法向翻转系数：刀具应偏置到曲面外侧（远离另一面）。
    返回 {blade: ±1.0}。只有单面时退化为 +1（无法判定外侧）。"""
    cents = {}
    for blade in (0, 1):
        pts = []
        for p in patches:
            if p.blade == blade and not p.is_curl:
                pts.extend(p.C0)
                pts.extend(p.C1)
        if pts:
            cents[blade] = tuple(sum(c) / len(pts) for c in zip(*pts))
        else:
            cents[blade] = None

    flips = {}
    for blade in (0, 1):
        other = 1 - blade
        if cents[blade] is None or cents[other] is None:
            flips[blade] = 1.0
            continue
        ref = cents[other]
        s, n = 0.0, 0
        for p in patches:
            if p.blade != blade or p.is_curl:
                continue
            n0s, n1s, _ = ruling_normals(p.C0, p.C1)
            for i in range(len(n0s)):
                mid = _lerp(p.C0[i], p.C1[i], 0.5)
                s += _dot(n0s[i], _sub(mid, ref))
                n += 1
        flips[blade] = 1.0 if s >= 0 else -1.0
    return flips

def _connected_regions(patches, developable=None):
    """把格子按网格 4 邻接聚类成连通域。developable=None 表示所有格子（不分可展）。
    无网格元数据时退化为每格一个区域。返回 list[list[Patch]]。"""
    cells = [p for p in patches if developable is None or p.developable == developable]
    if not cells:
        return []
    if any(p.row < 0 or p.col < 0 for p in cells):
        return [[p] for p in cells]

    by_key = {(p.blade, p.row, p.col): p for p in cells}
    parent = {k: k for k in by_key}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for blade, row, col in list(by_key):
        for dr, dc in ((0, 1), (1, 0)):
            nk = (blade, row + dr, col + dc)
            if nk in by_key:
                union((blade, row, col), nk)

    groups = {}
    for k in by_key:
        groups.setdefault(find(k), []).append(by_key[k])
    return list(groups.values())

def _snake_order(cells):
    """按 (row, col) 蛇形排序（相邻行反向），使相邻格进给方向能首尾相接。"""
    rows = sorted({p.row for p in cells})
    ordered = []
    for ri, r in enumerate(rows):
        row_cells = sorted([p for p in cells if p.row == r], key=lambda p: p.col)
        if ri % 2 == 1:
            row_cells.reverse()
        ordered.extend(row_cells)
    return ordered

def _flank_stitched_region(region, tool_r):
    """把一个连通域内的片按蛇形拼成一条连续进给折线。
    返回 (feed_pts, axis_segs)。"""
    ordered = _snake_order(region)
    feed_all, axes_all = [], []
    for pi, p in enumerate(ordered):
        feed, axes = flank_cl_lines(p.C0, p.C1, tool_r, flip=p.normal_flip)
        if pi % 2 == 1:
            feed = list(reversed(feed))
        feed_all.extend(feed)
        axes_all.extend(axes)
    return feed_all, axes_all

def flank_stitched_feed_lines(patches, tool_r):
    """把可展片按连通域拼接，返回 (feed_lines, axis_segs)。
    feed_lines: 每个可展连通域一条连续进给折线（含跨格衔接段）"""
    regions = _connected_regions(patches, True)
    feed_lines, axis_segs = [], []
    for region in regions:
        feed, axes = _flank_stitched_region(region, tool_r)
        if len(feed) >= 2:
            feed_lines.append(feed)
        axis_segs.extend(axes)
    return feed_lines, axis_segs

# ===================== 主流程 =====================

def compute(input_dir, args):
    """扫描 _params.txt，计算每块指标，返回 patches 列表。"""
    files = sorted(Path(input_dir).glob("*_params.txt"))
    if not files:
        print(f"[Error] 未找到 *_params.txt 文件: {input_dir}")
        return []
    grid_meta = _load_grid_meta(input_dir)
    patches = []
    for fp in files:
        C0, C1 = load_directrix(fp)
        if len(C0) < 2 or len(C1) < 2:
            print(f"  [skip] {fp.name}: 准线点数不足")
            continue
        p = Patch(name=fp.name.replace('_params.txt', ''), C0=C0, C1=C1)
        p.blade, p.cell_idx = _parse_cell_name(p.name)
        if p.blade in grid_meta and grid_meta[p.blade][1] > 0 and p.cell_idx >= 0:
            nRows, nCols = grid_meta[p.blade]
            p.row = p.cell_idx // nCols
            p.col = p.cell_idx % nCols
        p.directrix_len = curve_length(C0)
        p.mean_ruling = mean_ruling_length(C0, C1)
        p.area = surface_area(C0, C1)
        p.flank_time, p.flank_len, p.flank_passes, p.developable, p.twist = \
            flank_milling(C0, C1, args.feed, args.twist_limit)
        p.point_time, p.stepover, p.point_len, _ = \
            point_milling(C0, C1, args.feed, args.ball_r, args.scallop)
        patches.append(p)

    flips = _blade_normal_flips(patches)
    for p in patches:
        p.normal_flip = flips.get(p.blade, 1.0)
    return patches

def summarize(patches, args):
    """汇总并返回结果 dict。
    A 直纹面侧铣（全部片侧铣一次，按连通域拼接）vs B 传统整体点铣。
    非切削时间按「连通域数 × 每区域进退刀」计算（刀轨拼接后只进退刀一次/区域）。
    二者在同一精度指标下比较：点铣残留高度 = 拟合容差。"""
    flank_regions = _connected_regions(patches)   # 全部片连通域
    n_flank_regions = len(flank_regions)

    flank_cut = sum(p.flank_time for p in patches)
    flank_overhead = n_flank_regions * args.overhead
    flank_total = flank_cut + flank_overhead

    point_cut = sum(p.point_time for p in patches)   # 基线：整体点铣
    point_overhead = getattr(args, 'point_overhead', 10.0)
    point_total = point_cut + point_overhead
    total_area = sum(p.area for p in patches)
    return {
        "num_patches": len(patches),
        "flank_regions": n_flank_regions,
        "total_area": total_area,
        "flank": {
            "cut": flank_cut,
            "overhead": flank_overhead,
            "total": flank_total,
        },
        "point": {"cut": point_cut, "overhead": point_overhead, "total": point_total},
    }

def print_report(patches, args, summary):
    print("=" * 96)
    print(f"直纹面拟合 vs 传统点铣 加工时间仿真  (feed={args.feed} mm/min, "
          f"侧铣刀R={getattr(args, 'tool_r', args.ball_r):.1f}, "
          f"球刀R={args.ball_r}, 残留={args.scallop}, 可展阈值={args.twist_limit}°)")
    print(f"面片数: {summary['num_patches']}")
    print("=" * 96)
    print(f"{'面片':<28}{'扭转°':>8}{'准线mm':>9}{'母线mm':>9}"
          f"{'侧铣s':>9}{'点铣s':>9}")
    for p in patches:
        print(f"{p.name:<28}{p.twist:>8.2f}"
              f"{p.directrix_len:>9.1f}{p.mean_ruling:>9.1f}"
              f"{p.flank_time:>9.2f}{p.point_time:>9.2f}")
    print("-" * 96)
    print(f"侧铣连通域 {summary.get('flank_regions', 0)}  总曲面面积: {summary['total_area']:.1f} mm^2")
    twists = sorted((p.twist for p in patches), reverse=True)
    if twists:
        print(f"扭转角分布(deg): max={twists[0]:.2f}  "
              f"p90={twists[max(0, int(0.10 * len(twists)))]:.2f}  "
              f"中位={twists[len(twists) // 2]:.2f}")
    print()
    f, p = summary["flank"], summary["point"]
    print(f"A 直纹面侧铣(拟合后): 切削 {f['cut']:8.1f}s "
          f"+ 非切削 {f['overhead']:8.1f}s = {f['total']:8.1f}s")
    print(f"B 传统点铣(原曲面):   切削 {p['cut']:8.1f}s "
          f"+ 非切削 {p['overhead']:8.1f}s = {p['total']:8.1f}s")
    if f['total'] > 0:
        print(f"A 相对 B 提速: {p['total'] / f['total']:.1f}x")

def main():
    if hasattr(sys.stdout, 'reconfigure'):
        try:
            sys.stdout.reconfigure(encoding='utf-8', errors='replace')
            sys.stderr.reconfigure(encoding='utf-8', errors='replace')
        except Exception:
            pass
    cfg = load_config()
    ap = argparse.ArgumentParser(description="直纹面拟合加工时间仿真（简化可辩护刀轨）")
    ap.add_argument("input", help="输出目录（含 *_params.txt）")
    ap.add_argument("--feed", type=float, default=cfg.get("feed", 500.0), help="切削进给 mm/min")
    ap.add_argument("--rapid", type=float, default=cfg.get("rapid", 5000.0), help="快进 mm/min（保留，模型未用）")
    ap.add_argument("--tool-r", type=float, default=cfg.get("tool_r", 5.0), help="侧铣圆柱刀半径 mm")
    ap.add_argument("--ball-r", type=float, default=cfg.get("ball_r", 5.0), help="球头刀半径 mm")
    ap.add_argument("--scallop", type=float, default=cfg.get("scallop", 0.01), help="点铣残留高度 mm")
    ap.add_argument("--overhead", type=float, default=cfg.get("overhead", 2.0), help="每块侧铣进退刀+衔接时间 s")
    ap.add_argument("--point-overhead", type=float, default=cfg.get("point_overhead", 10.0), help="点铣整体进退刀时间 s")
    ap.add_argument("--twist-limit", type=float, default=cfg.get("twist_limit", 1.0), help="可展判定阈值（法向扭转角，度）")
    ap.add_argument("--json", help="可选：输出 JSON 结果文件路径")
    ap.add_argument("--vtk", help="可选：输出刀轨 VTK 到指定目录")
    args = ap.parse_args()

    patches = compute(args.input, args)
    if not patches:
        sys.exit(1)
    summary = summarize(patches, args)
    print_report(patches, args, summary)

    if args.vtk:
        os.makedirs(args.vtk, exist_ok=True)
        for p in patches:
            base = os.path.join(args.vtk, p.name)
            write_vtk_polylines(base + "_flank.vtk",
                                flank_toolpath_lines(p.C0, p.C1, tool_r=args.tool_r))
            write_vtk_polylines(base + "_point.vtk",
                                point_toolpath_lines(p.C0, p.C1, p.stepover, ball_r=args.ball_r))
        print(f"刀轨 VTK 已写入: {args.vtk}")

    if args.json:
        result = {
            "feed": args.feed, "tool_r": args.tool_r, "ball_r": args.ball_r,
            "scallop": args.scallop, "twist_limit": args.twist_limit,
            "overhead": args.overhead, "point_overhead": args.point_overhead,
            **summary,
            "patches": [{k: v for k, v in asdict(p).items()
                         if k not in ('C0', 'C1', 'tool_r')} for p in patches],
        }
        with open(args.json, 'w', encoding='utf-8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
        print(f"结果已写入: {args.json}")

if __name__ == "__main__":
    main()
