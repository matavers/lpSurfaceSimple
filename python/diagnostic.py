"""Objective post-processing diagnostics for Simple ruled-surface fitting output.

Reads the OBJ files already exported by simple.exe and measures quality purely
from the data (it never re-runs the fitting algorithm):

  fit error    : composite -> original mesh distance (accuracy)
  coverage     : original mesh -> composite distance (gaps / missing patches)
  normal dev.  : composite normals vs original normals (G1 proxy)
  census       : component counts + area ratio

The composite is trim+strip+corner when blend is enabled, otherwise the raw
ruled cells. The original mesh (blade*_mesh.obj) is the ground truth.

Imported by gui_main.py as the "Diagnose" action; also runnable standalone:

    python diagnostic.py <output_dir>
"""

import os
import sys
import json
from dataclasses import dataclass, field
from typing import Dict, List

import numpy as np


def _pv():
    import pyvista as pv
    return pv


@dataclass
class BladeDiagnostic:
    blade_index: int = 0
    n_trim: int = 0
    n_strip: int = 0
    n_corner: int = 0
    n_cell: int = 0
    composite_points: int = 0
    composite_tris: int = 0
    composite_area: float = 0.0
    original_area: float = 0.0
    area_ratio: float = 0.0
    fit_max: float = 0.0
    fit_mean: float = 0.0
    fit_rms: float = 0.0
    coverage_max: float = 0.0
    coverage_mean: float = 0.0
    coverage_pct_beyond_tol: float = 0.0
    normal_dev_max: float = 0.0
    normal_dev_mean: float = 0.0
    shared_mesh: bool = False


@dataclass
class DiagnosticResult:
    blades: List[BladeDiagnostic] = field(default_factory=list)
    gap_tol: float = 0.05
    kink_angle: float = 20.0
    fit_meshes: Dict[int, object] = field(default_factory=dict)
    normal_dev_meshes: Dict[int, object] = field(default_factory=dict)
    coverage_meshes: Dict[int, object] = field(default_factory=dict)


def _merge_tri(meshes):
    """Concatenate triangulated PolyData into one PolyData."""
    pv = _pv()
    pts, tris, off = [], [], 0
    for m in meshes:
        t = m.triangulate()
        if t.n_points == 0:
            continue
        pts.append(np.asarray(t.points, dtype=np.float64))
        f = np.asarray(t.faces).reshape(-1, 4)
        tris.append(f[:, 1:4].astype(np.int64) + off)
        off += t.n_points
    if not pts:
        return pv.PolyData()
    points = np.vstack(pts)
    faces = np.hstack([
        np.full((np.vstack(tris).shape[0], 1), 3, dtype=np.int64),
        np.vstack(tris),
    ]).ravel()
    return pv.PolyData(points, faces)


def collect_components(out_dir):
    """Group exported OBJ files per blade (0/1) by role."""
    files = sorted(os.listdir(out_dir))
    blades = {0: {'mesh': None, 'trim': [], 'strip': [], 'corner': [], 'cell': []},
              1: {'mesh': None, 'trim': [], 'strip': [], 'corner': [], 'cell': []}}
    for fn in files:
        if not fn.endswith('.obj'):
            continue
        if fn.startswith('blade1'):
            bi = 0
        elif fn.startswith('blade2'):
            bi = 1
        else:
            continue
        b = blades[bi]
        if fn.endswith('_mesh.obj'):
            b['mesh'] = fn
        elif '_trim_' in fn:
            b['trim'].append(fn)
        elif '_strip' in fn:
            b['strip'].append(fn)
        elif '_corner_' in fn:
            b['corner'].append(fn)
        elif '_cell' in fn:
            b['cell'].append(fn)
    return blades


def _dist(mesh, surface):
    """Absolute point-to-surface distance of `mesh` vertices to `surface`."""
    out = mesh.compute_implicit_distance(surface, inplace=False)
    return np.abs(np.asarray(out['implicit_distance']))


def run_diagnostic(out_dir, gap_tol=0.05):
    pv = _pv()
    blades = collect_components(out_dir)
    result = DiagnosticResult(gap_tol=gap_tol)

    # Load ground-truth meshes first; detect the shared-mesh case
    # (faceIdx1 == faceIdx2 -> blade2_mesh.obj is a copy of blade1_mesh.obj).
    originals = {}
    shared = False
    for bi in (0, 1):
        if blades[bi]['mesh'] is None:
            continue
        m = pv.read(os.path.join(out_dir, blades[bi]['mesh']))
        originals[bi] = m
    if 0 in originals and 1 in originals:
        p0 = np.asarray(originals[0].points)
        p1 = np.asarray(originals[1].points)
        if p0.shape == p1.shape and np.allclose(p0, p1):
            shared = True

    for bi in (0, 1):
        b = blades[bi]
        if bi not in originals:
            continue
        original = originals[bi]
        if b['trim']:
            part_files = b['trim'] + b['strip'] + b['corner']
            census = (len(b['trim']), len(b['strip']), len(b['corner']), 0)
        else:
            part_files = b['cell']
            census = (0, 0, 0, len(b['cell']))
        if not part_files:
            continue

        parts = [pv.read(os.path.join(out_dir, f)) for f in part_files]
        composite = _merge_tri(parts)

        d = BladeDiagnostic(
            blade_index=bi,
            n_trim=census[0], n_strip=census[1],
            n_corner=census[2], n_cell=census[3],
            shared_mesh=(shared and bi == 1),
        )
        d.composite_points = composite.n_points
        d.composite_tris = composite.n_cells

        try:
            d.composite_area = composite.area
            d.original_area = original.area
            d.area_ratio = d.composite_area / d.original_area if d.original_area > 0 else 0.0
        except Exception:
            pass

        fit = _dist(composite, original)
        d.fit_max = float(fit.max())
        d.fit_mean = float(fit.mean())
        d.fit_rms = float(np.sqrt(np.mean(fit ** 2)))

        cov = _dist(original, composite)
        d.coverage_max = float(cov.max())
        d.coverage_mean = float(cov.mean())
        d.coverage_pct_beyond_tol = float(np.mean(cov > gap_tol) * 100.0)

        dev = None
        try:
            comp_n = composite.compute_normals(point_normals=True, cell_normals=False)
            orig_cell = original.compute_normals(point_normals=False, cell_normals=True)
            idx = original.find_closest_cell(composite.points)
            orig_n = np.asarray(orig_cell.cell_data['Normals'])[idx]
            dot = np.einsum('ij,ij->i', np.asarray(comp_n.point_data['Normals']), orig_n)
            dot = np.clip(np.abs(dot), 0.0, 1.0)
            dev = np.degrees(np.arccos(dot))
            d.normal_dev_max = float(dev.max())
            d.normal_dev_mean = float(dev.mean())
        except Exception:
            pass

        fit_mesh = composite.copy()
        fit_mesh['fit_error'] = fit
        result.fit_meshes[bi] = fit_mesh

        cov_mesh = original.copy()
        cov_mesh['coverage'] = cov
        result.coverage_meshes[bi] = cov_mesh

        if dev is not None:
            nd_mesh = composite.copy()
            nd_mesh['normal_dev'] = dev
            result.normal_dev_meshes[bi] = nd_mesh

        result.blades.append(d)

    return result


def report_text(result):
    lines = ["objective diagnostics (gap_tol=%.3f mm)" % result.gap_tol]
    for d in result.blades:
        lines.append(
            "Blade %d: trim=%d strip=%d corner=%d cell=%d"
            % (d.blade_index + 1, d.n_trim, d.n_strip, d.n_corner, d.n_cell))
        lines.append(
            "  composite %d pts / %d tris, area=%.3f vs original %.3f (ratio=%.4f)"
            % (d.composite_points, d.composite_tris,
               d.composite_area, d.original_area, d.area_ratio))
        lines.append(
            "  fit error   max=%.4f  mean=%.4f  rms=%.4f mm"
            % (d.fit_max, d.fit_mean, d.fit_rms))
        lines.append(
            "  coverage    max=%.4f  mean=%.4f mm  (%.2f%% points > %.3f mm)"
            % (d.coverage_max, d.coverage_mean,
               d.coverage_pct_beyond_tol, result.gap_tol))
        lines.append(
            "  normal dev  max=%.2f  mean=%.2f deg"
            % (d.normal_dev_max, d.normal_dev_mean))
        if d.shared_mesh:
            lines.append(
                "  [WARN] ground-truth mesh is a copy of blade 1's "
                "(same STEP face); comparison not meaningful")
    return "\n".join(lines)


def write_report(result, out_dir):
    data = {'gap_tol': result.gap_tol, 'blades': []}
    for d in result.blades:
        data['blades'].append({
            'blade': d.blade_index + 1,
            'census': {'trim': d.n_trim, 'strip': d.n_strip,
                       'corner': d.n_corner, 'cell': d.n_cell},
            'composite_points': d.composite_points,
            'composite_tris': d.composite_tris,
            'composite_area': d.composite_area,
            'original_area': d.original_area,
            'area_ratio': d.area_ratio,
            'fit_error': {'max': d.fit_max, 'mean': d.fit_mean, 'rms': d.fit_rms},
            'coverage': {'max': d.coverage_max, 'mean': d.coverage_mean,
                         'pct_beyond_tol': d.coverage_pct_beyond_tol},
            'normal_dev': {'max': d.normal_dev_max, 'mean': d.normal_dev_mean},
            'shared_mesh': d.shared_mesh,
        })
    path = os.path.join(out_dir, 'diagnostic.json')
    with open(path, 'w') as f:
        json.dump(data, f, indent=2)
    return path


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("usage: python diagnostic.py <output_dir> [gap_tol_mm]")
        sys.exit(1)
    out = sys.argv[1]
    tol = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
    res = run_diagnostic(out, gap_tol=tol)
    print(report_text(res))
    print("wrote", write_report(res, out))
