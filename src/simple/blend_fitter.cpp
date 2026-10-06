#include "simple/blend_fitter.hpp"

#include "simple/ruled_fitter.hpp"

#include <Geom_BSplineCurve.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>

#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

namespace simple {

namespace {

// ── 分解函数（单位分解基函数，与 分解函数选取.docx 一致） ──────────────
// 二项式系数 C(n,k)。
double binom(int n, int k) {
    if (k < 0 || k > n) return 0.0;
    double r = 1.0;
    for (int i = 1; i <= k; ++i) r = r * (n - k + i) / i;
    return r;
}

// 一维过渡函数 s_k(t)（最低次多项式解，即 smootherstep，分解函数选取.docx §2.3）：
//   满足 s(0)=0, s(1)=1, s^(m)(0)=s^(m)(1)=0 (m=1..k)，从而整体 C^k 连续。
//   s_k(t) = t^{k+1} · Σ_{i=0..k} C(k+i, i) (1-t)^i
//     k=0 → t                   （C⁰，线性）
//     k=1 → 3t²-2t³             （C¹）
//     k=2 → 6t⁵-15t⁴+10t³       （C²，五次 smootherstep，默认）
//     k=3 → 35t⁴-84t⁵+70t⁶-20t⁷（C³）
// 该 Bernstein 型写法系数全为正，数值稳定；是 B 样条构造的最低次解（唯一）。
double smootherstep(double t, int order) {
    if (t <= 0.0) return 0.0;
    if (t >= 1.0) return 1.0;
    if (order <= 0) return t;
    double omt = 1.0 - t;
    double sum = 0.0;
    for (int i = 0; i <= order; ++i)
        sum += binom(order + i, i) * std::pow(omt, i);
    return std::pow(t, order + 1) * sum;
}

// 一维基函数 chi(u; uL, uR, a, b, order)：
//   u ≤ uL-a 或 u ≥ uR+b  → 0
//   uL ≤ u ≤ uR           → 1（恒一区域）
//   两侧过渡带内用 order 阶 smootherstep 单调过渡。
double chi1(double u, double uL, double uR, double a, double b, int order) {
    if (u <= uL - a || u >= uR + b) return 0.0;
    if (u < uL) return smootherstep((u - (uL - a)) / a, order);
    if (u <= uR) return 1.0;
    return 1.0 - smootherstep((u - uR) / b, order);
}

// 单格直纹面求值：准线方向截断到格胞内（外推区常数延拓），
// 母线方向线性延拓。过渡带内 φ_i 平滑衰减到 0，乘积 φ_i·R_i 仍 C^order。
Vec3 evalRuledCellSurface(const GridCell& cell, double u, double v) {
    const Vec3Arr& c0 = cell.ruled.curveC0Samples;
    const Vec3Arr& c1 = cell.ruled.curveC1Samples;
    int n = static_cast<int>(c0.size());
    if (n < 2 || static_cast<int>(c1.size()) != n)
        return Vec3(0, 0, 0);

    if (cell.fitDir == ParamDir::V) {
        // 准线沿 u（索引由 u 决定，截断），母线沿 v（线性延拓）
        double uN = clamp((u - cell.u0) / (cell.u1 - cell.u0), 0.0, 1.0);
        double idx = uN * (n - 1);
        int i0 = static_cast<int>(std::floor(idx));
        int i1 = std::min(i0 + 1, n - 1);
        double fr = idx - i0;
        Vec3 C0 = c0[i0] + (c0[i1] - c0[i0]) * fr;
        Vec3 C1 = c1[i0] + (c1[i1] - c1[i0]) * fr;
        double t = (v - cell.v0) / (cell.v1 - cell.v0);
        return C0 + (C1 - C0) * t;
    } else {
        // 准线沿 v（截断），母线沿 u（线性延拓）
        double vN = clamp((v - cell.v0) / (cell.v1 - cell.v0), 0.0, 1.0);
        double idx = vN * (n - 1);
        int i0 = static_cast<int>(std::floor(idx));
        int i1 = std::min(i0 + 1, n - 1);
        double fr = idx - i0;
        Vec3 C0 = c0[i0] + (c0[i1] - c0[i0]) * fr;
        Vec3 C1 = c1[i0] + (c1[i1] - c1[i0]) * fr;
        double s = (u - cell.u0) / (cell.u1 - cell.u0);
        return C0 + (C1 - C0) * s;
    }
}

// 过渡带宽度：用最小格胞宽度计算全局统一的 aU/aV。
// 这样任意格胞的邻域都不会越过相邻格胞（满足外推条件）。
void bandWidths(const GridResult& gr, const BlendConfig& cfg, double& aU, double& aV) {
    double minU = std::numeric_limits<double>::max();
    double minV = std::numeric_limits<double>::max();
    for (const auto& cell : gr.cells) {
        minU = std::min(minU, cell.u1 - cell.u0);
        minV = std::min(minV, cell.v1 - cell.v0);
    }
    if (!std::isfinite(minU) || minU <= 0.0) minU = 1.0;
    if (!std::isfinite(minV) || minV <= 0.0) minV = 1.0;
    aU = cfg.bandWidth * minU;
    aV = cfg.bandWidth * minV;
}

struct ErrorStat {
    double maxError = 0.0;
    double sumSq = 0.0;
    int count = 0;
};

void addSample(ErrorStat& e, double d) {
    e.maxError = std::max(e.maxError, d);
    e.sumSq += d * d;
    e.count++;
}

double rmsOf(const ErrorStat& e) {
    return e.count ? std::sqrt(e.sumSq / e.count) : 0.0;
}

} // namespace

Vec3 evalBlend(const GridResult& gr, const BlendConfig& cfg, double u, double v) {
    double aU = 0.0, aV = 0.0;
    bandWidths(gr, cfg, aU, aV);

    double sumPhi = 0.0;
    Vec3 acc(0, 0, 0);
    for (const auto& cell : gr.cells) {
        double cu = chi1(u, cell.u0, cell.u1, aU, aU, cfg.continuity);
        if (cu <= 0.0) continue;
        double cv = chi1(v, cell.v0, cell.v1, aV, aV, cfg.continuity);
        if (cv <= 0.0) continue;
        double phi = cu * cv;
        acc += phi * evalRuledCellSurface(cell, u, v);
        sumPhi += phi;
    }

    if (sumPhi > 1e-12) return acc / sumPhi;

    // 兜底：采样点落在邻域覆盖之外（理论上不会发生），退化为最近格胞直纹面
    for (const auto& cell : gr.cells)
        if (u >= cell.u0 && u <= cell.u1 && v >= cell.v0 && v <= cell.v1)
            return evalRuledCellSurface(cell, u, v);
    return gr.cells.empty() ? Vec3(0, 0, 0) : evalRuledCellSurface(gr.cells.front(), u, v);
}

BlendResult computeBlendError(const SurfaceWrapper& surf, const GridResult& gr,
                              const BlendConfig& cfg) {
    BlendResult br;

    auto [uMin, uMax] = surf.paramDomainU();
    auto [vMin, vMax] = surf.paramDomainV();

    ErrorStat e;
    int sU = std::max(2, cfg.nSampleU);
    int sV = std::max(2, cfg.nSampleV);
    for (int i = 0; i < sU; ++i) {
        double u = uMin + (uMax - uMin) * i / (sU - 1.0);
        for (int j = 0; j < sV; ++j) {
            double v = vMin + (vMax - vMin) * j / (sV - 1.0);
            addSample(e, (evalBlend(gr, cfg, u, v) - surf.evaluate(u, v)).norm());
        }
    }
    br.maxError = e.maxError;
    br.rmsError = rmsOf(e);
    return br;
}

bool exportTransitionBandVTK(const std::string& path,
                             const GridResult& gr,
                             const BlendConfig& cfg,
                             int nPerEdge) {
    if (nPerEdge < 2) nPerEdge = 2;
    double aU = 0.0, aV = 0.0;
    bandWidths(gr, cfg, aU, aV);

    struct Line { std::vector<Vec3> pts; };
    std::vector<Line> lines;

    for (const auto& cell : gr.cells) {
        // 内缩核心区：内部边向内缩 aU/aV，外边界不缩。
        double u0 = (cell.col > 0)            ? cell.u0 + aU : cell.u0;
        double u1 = (cell.col < gr.nCols - 1) ? cell.u1 - aU : cell.u1;
        double v0 = (cell.row > 0)            ? cell.v0 + aV : cell.v0;
        double v1 = (cell.row < gr.nRows - 1) ? cell.v1 - aV : cell.v1;
        if (u1 <= u0 || v1 <= v0) continue;   // 过渡带过宽导致核心退化

        Line bottom, right, top, left;
        for (int i = 0; i < nPerEdge; ++i) {
            double u = u0 + (u1 - u0) * i / (nPerEdge - 1.0);
            bottom.pts.push_back(evalRuledCellSurface(cell, u, v0));
            top.pts.push_back(evalRuledCellSurface(cell, u1 - (u1 - u0) * i / (nPerEdge - 1.0), v1));
        }
        for (int j = 0; j < nPerEdge; ++j) {
            double v = v0 + (v1 - v0) * j / (nPerEdge - 1.0);
            right.pts.push_back(evalRuledCellSurface(cell, u1, v));
            left.pts.push_back(evalRuledCellSurface(cell, u0, v1 - (v1 - v0) * j / (nPerEdge - 1.0)));
        }
        lines.push_back(std::move(bottom));
        lines.push_back(std::move(right));
        lines.push_back(std::move(top));
        lines.push_back(std::move(left));
    }

    std::ofstream out(path);
    if (!out) return false;

    int nPts = 0, nSegs = 0;
    for (const auto& L : lines) {
        nPts += static_cast<int>(L.pts.size());
        if (L.pts.size() >= 2) nSegs += static_cast<int>(L.pts.size()) - 1;
    }

    out << "# vtk DataFile Version 3.0\n";
    out << "transition band (retracted core boundary)\n";
    out << "ASCII\nDATASET POLYDATA\n";
    out << "POINTS " << nPts << " float\n";
    out << std::fixed << std::setprecision(6);
    for (const auto& L : lines)
        for (const auto& p : L.pts)
            out << p.x() << " " << p.y() << " " << p.z() << "\n";
    out << "LINES " << nSegs << " " << nSegs * 3 << "\n";
    int base = 0;
    for (const auto& L : lines) {
        for (int i = 0; i + 1 < static_cast<int>(L.pts.size()); ++i)
            out << "2 " << (base + i) << " " << (base + i + 1) << "\n";
        base += static_cast<int>(L.pts.size());
    }
    return true;
}

namespace {

std::string bsplineToJson(const Handle(Geom_BSplineCurve)& c) {
    std::ostringstream o;
    o << std::fixed << std::setprecision(12);
    if (c.IsNull()) return "{\"degree\":0,\"knots\":[],\"ctrl\":[]}";
    int degree = c->Degree();
    int nCtrl = c->NbPoles();
    TColgp_Array1OfPnt poles(1, nCtrl);
    c->Poles(poles);
    int nKnots = c->NbKnots();
    TColStd_Array1OfReal knots(1, nKnots);
    c->Knots(knots);
    TColStd_Array1OfInteger mults(1, nKnots);
    c->Multiplicities(mults);
    o << "{\"degree\":" << degree << ",\"knots\":[";
    bool first = true;
    for (int k = 1; k <= nKnots; ++k) {
        for (int m = 0; m < mults(k); ++m) {
            if (!first) o << ",";
            first = false;
            o << knots(k);
        }
    }
    o << "],\"ctrl\":[";
    first = true;
    for (int i = 1; i <= nCtrl; ++i) {
        if (!first) o << ",";
        first = false;
        o << "[" << poles(i).X() << "," << poles(i).Y() << "," << poles(i).Z() << "]";
    }
    o << "]}";
    return o.str();
}

// 由优化后的准线采样点重建 B 样条（参数范围 [p0,p1] 与格胞参数域一致），
// 供 JSON 导出使用。这样导出的准线才是「拟合后」的准线，而非原始等参线，
// 与 exportGridOBJs / exportBlendSurfaceOBJ 使用的优化准线保持一致。
Handle(Geom_BSplineCurve) samplesToBSpline(const Vec3Arr& samples,
                                           double p0, double p1,
                                           int degree = 3) {
    int n = static_cast<int>(samples.size());
    if (n < 2) return Handle(Geom_BSplineCurve)();
    TColgp_Array1OfPnt pts(1, n);
    TColStd_Array1OfReal params(1, n);
    for (int i = 0; i < n; ++i) {
        pts.SetValue(i + 1, gp_Pnt(samples[i].x(), samples[i].y(), samples[i].z()));
        params.SetValue(i + 1, p0 + (p1 - p0) * i / (n - 1.0));
    }
    GeomAPI_PointsToBSpline fit(pts, params, degree, degree, GeomAbs_C2, 1e-3);
    if (!fit.IsDone()) return Handle(Geom_BSplineCurve)();
    return fit.Curve();
}

} // namespace

bool exportSurfaceModelJson(const std::string& path,
                            const GridResult& gr,
                            const BlendConfig& cfg) {
    std::ofstream o(path);
    if (!o) return false;
    o << std::fixed << std::setprecision(12);
    o << "{\"name\":\"" << gr.name << "\"";
    o << ",\"nRows\":" << gr.nRows << ",\"nCols\":" << gr.nCols;
    o << ",\"fitDir\":\"" << (gr.fitDir == ParamDir::U ? "U" : "V") << "\"";
    o << ",\"uEdges\":[";
    for (size_t i = 0; i < gr.uEdges.size(); ++i) { if (i) o << ","; o << gr.uEdges[i]; }
    o << "],\"vEdges\":[";
    for (size_t i = 0; i < gr.vEdges.size(); ++i) { if (i) o << ","; o << gr.vEdges[i]; }
    o << "],\"partition\":{\"type\":\"smootherstep\",\"bandWidth\":" << cfg.bandWidth
      << ",\"continuity\":" << cfg.continuity << "}";
    o << ",\"cells\":[";
    bool firstCell = true;
    for (const auto& c : gr.cells) {
        if (!firstCell) o << ",";
        firstCell = false;
        o << "{\"row\":" << c.row << ",\"col\":" << c.col;
        o << ",\"u0\":" << c.u0 << ",\"u1\":" << c.u1;
        o << ",\"v0\":" << c.v0 << ",\"v1\":" << c.v1;
        Handle(Geom_BSplineCurve) c0bs, c1bs;
        if (gr.fitDir == ParamDir::U) {
            c0bs = samplesToBSpline(c.ruled.curveC0Samples, c.v0, c.v1);
            c1bs = samplesToBSpline(c.ruled.curveC1Samples, c.v0, c.v1);
        } else {
            c0bs = samplesToBSpline(c.ruled.curveC0Samples, c.u0, c.u1);
            c1bs = samplesToBSpline(c.ruled.curveC1Samples, c.u0, c.u1);
        }
        o << ",\"c0\":" << bsplineToJson(c0bs.IsNull() ? c.ruled.curveC0 : c0bs);
        o << ",\"c1\":" << bsplineToJson(c1bs.IsNull() ? c.ruled.curveC1 : c1bs);
        o << "}";
    }
    o << "]}";
    return true;
}

bool exportBlendSurfaceOBJ(const std::string& path,
                           const GridResult& gr,
                           const BlendConfig& cfg,
                           int nU, int nV) {
    if (gr.uEdges.empty() || gr.vEdges.empty()) return false;
    double u0 = gr.uEdges.front(), u1 = gr.uEdges.back();
    double v0 = gr.vEdges.front(), v1 = gr.vEdges.back();
    nU = std::max(2, nU);
    nV = std::max(2, nV);
    Vec3Arr verts((nU + 1) * (nV + 1));
    FaceArr faces;
    for (int i = 0; i <= nU; ++i) {
        double u = u0 + (u1 - u0) * i / nU;
        for (int j = 0; j <= nV; ++j) {
            double v = v0 + (v1 - v0) * j / nV;
            verts[i * (nV + 1) + j] = evalBlend(gr, cfg, u, v);
        }
    }
    for (int i = 0; i < nU; ++i) {
        for (int j = 0; j < nV; ++j) {
            int a = i * (nV + 1) + j;
            faces.push_back({a, a + 1, a + nV + 1});
            faces.push_back({a + 1, a + nV + 2, a + nV + 1});
        }
    }
    return exportOBJ(path, verts, faces);
}

} // namespace simple
