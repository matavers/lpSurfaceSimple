#include "simple/blend_fitter.hpp"

#include "simple/ruled_fitter.hpp"

#include <cmath>
#include <algorithm>

namespace simple {

namespace {

// ── 工具 ──────────────────────────────────────────────────────

struct ErrorStat {
    double maxError = 0.0;
    double rmsError = 0.0;
    double sumSq = 0.0;
    int count = 0;
};

void addSample(ErrorStat& e, double d) {
    e.maxError = std::max(e.maxError, d);
    e.sumSq += d * d;
    e.count++;
}

void finalize(ErrorStat& e) {
    e.rmsError = e.count ? std::sqrt(e.sumSq / e.count) : 0.0;
}

Vec3 lerp3(const Vec3& a, const Vec3& b, double t) { return a + (b - a) * t; }

// Hermite 基函数
double h00(double t) { return 2 * t * t * t - 3 * t * t + 1; }
double h10(double t) { return -2 * t * t * t + 3 * t * t; }
double h01(double t) { return t * t * t - 2 * t * t + t; }
double h11(double t) { return t * t * t - t * t; }

// 在格子直纹面上求值 P 及其偏导 ∂P/∂u、∂P/∂v
void evalRuledCell(const GridCell& cell, double u, double v,
                   Vec3& P, Vec3& dU, Vec3& dV)
{
    const Vec3Arr& c0 = cell.ruled.curveC0Samples;
    const Vec3Arr& c1 = cell.ruled.curveC1Samples;
    int n = static_cast<int>(c0.size());
    if (n < 2) {
        P = Vec3(0, 0, 0); dU = Vec3(0, 0, 0); dV = Vec3(0, 0, 0);
        return;
    }

    if (cell.fitDir == ParamDir::V) {
        // 准线沿 u，母线沿 v
        double uN = (u - cell.u0) / (cell.u1 - cell.u0);
        uN = clamp(uN, 0.0, 1.0);
        double idx = uN * (n - 1);
        int i0 = (int)std::floor(idx);
        int i1 = std::min(i0 + 1, n - 1);
        double fr = idx - i0;
        Vec3 C0 = lerp3(c0[i0], c0[i1], fr);
        Vec3 C1 = lerp3(c1[i0], c1[i1], fr);
        double t = (v - cell.v0) / (cell.v1 - cell.v0);
        P = lerp3(C0, C1, t);

        double du = (cell.u1 - cell.u0) / (n - 1);
        Vec3 C0d = (c0[i1] - c0[i0]) / du;
        Vec3 C1d = (c1[i1] - c1[i0]) / du;
        dU = lerp3(C0d, C1d, t);
        dV = (C1 - C0) / (cell.v1 - cell.v0);
    } else {
        // 准线沿 v，母线沿 u
        double vN = (v - cell.v0) / (cell.v1 - cell.v0);
        vN = clamp(vN, 0.0, 1.0);
        double idx = vN * (n - 1);
        int i0 = (int)std::floor(idx);
        int i1 = std::min(i0 + 1, n - 1);
        double fr = idx - i0;
        Vec3 C0 = lerp3(c0[i0], c0[i1], fr);
        Vec3 C1 = lerp3(c1[i0], c1[i1], fr);
        double s = (u - cell.u0) / (cell.u1 - cell.u0);
        P = lerp3(C0, C1, s);

        double dv = (cell.v1 - cell.v0) / (n - 1);
        Vec3 C0d = (c0[i1] - c0[i0]) / dv;
        Vec3 C1d = (c1[i1] - c1[i0]) / dv;
        dV = lerp3(C0d, C1d, s);
        dU = (C1 - C0) / (cell.u1 - cell.u0);
    }
}

// 每格只内缩「内部边」，外边界不缩，保证复合面仍覆盖整个参数域
struct TrimmedBounds { double u0, u1, v0, v1; };

TrimmedBounds trimOf(const GridResult& gr, const GridCell& cell, double f) {
    TrimmedBounds t;
    double du = cell.u1 - cell.u0;
    double dv = cell.v1 - cell.v0;
    t.u0 = (cell.col > 0)             ? cell.u0 + f * du : cell.u0;
    t.u1 = (cell.col < gr.nCols - 1)  ? cell.u1 - f * du : cell.u1;
    t.v0 = (cell.row > 0)             ? cell.v0 + f * dv : cell.v0;
    t.v1 = (cell.row < gr.nRows - 1)  ? cell.v1 - f * dv : cell.v1;
    return t;
}

ErrorStat cellError(const SurfaceWrapper& surf, const GridResult& gr,
                    const GridCell& cell, double f, const BlendConfig& cfg) {
    ErrorStat e;
    TrimmedBounds t = trimOf(gr, cell, f);
    int n = cfg.nCellGrid;
    for (int i = 0; i < n; ++i) {
        double u = t.u0 + (t.u1 - t.u0) * i / (n - 1.0);
        for (int j = 0; j < n; ++j) {
            double v = t.v0 + (t.v1 - t.v0) * j / (n - 1.0);
            Vec3 P, dU, dV;
            evalRuledCell(cell, u, v, P, dU, dV);
            Vec3 S = surf.evaluate(u, v);
            addSample(e, (P - S).norm());
        }
    }
    finalize(e);
    return e;
}

// 生成 trimmed 格子渲染网格
void meshTrimmedCell(const GridResult& gr, const GridCell& cell, double f,
                     int res, Vec3Arr& verts, FaceArr& faces) {
    TrimmedBounds t = trimOf(gr, cell, f);
    verts.resize((res + 1) * (res + 1));
    for (int i = 0; i <= res; ++i) {
        double u = t.u0 + (t.u1 - t.u0) * i / res;
        for (int j = 0; j <= res; ++j) {
            double v = t.v0 + (t.v1 - t.v0) * j / res;
            Vec3 P, dU, dV;
            evalRuledCell(cell, u, v, P, dU, dV);
            verts[i * (res + 1) + j] = P;
        }
    }
    faces.clear();
    for (int i = 0; i < res; ++i)
        for (int j = 0; j < res; ++j) {
            int a = i * (res + 1) + j;
            faces.push_back({a, a + 1, a + res + 1});
            faces.push_back({a + 1, a + res + 2, a + res + 1});
        }
}

// 竖直缝（u = uEdges[c]）：左格 (r,c-1) 与右格 (r,c) 之间
BlendStrip buildVerticalStrip(const SurfaceWrapper& surf, const GridResult& gr,
                              const std::vector<double>& f,
                              int r, int c, const BlendConfig& cfg) {
    BlendStrip s;
    s.r = r; s.c = c; s.isVertical = true;

    const GridCell& L = gr.cells[r * gr.nCols + (c - 1)];
    const GridCell& R = gr.cells[r * gr.nCols + c];
    TrimmedBounds tL = trimOf(gr, L, f[r * gr.nCols + (c - 1)]);
    TrimmedBounds tR = trimOf(gr, R, f[r * gr.nCols + c]);

    double uL = tL.u1;   // 左格右边界
    double uR = tR.u0;   // 右格左边界
    double gapU = std::max(uR - uL, 1e-12);
    double vMin = gr.vEdges[r];
    double vMax = gr.vEdges[r + 1];

    int nA = cfg.nAlong, nB = cfg.nAcross;
    s.meshVerts.resize(nA * nB);
    ErrorStat e;

    for (int i = 0; i < nA; ++i) {
        double sv = (nA == 1) ? 0.5 : (double)i / (nA - 1);
        double v = vMin + (vMax - vMin) * sv;
        Vec3 PL, dUL, dVL;
        evalRuledCell(L, uL, v, PL, dUL, dVL);
        Vec3 PR, dUR, dVR;
        evalRuledCell(R, uR, v, PR, dUR, dVR);
        Vec3 TP = dUL * gapU;
        Vec3 TQ = dUR * gapU;

        for (int j = 0; j < nB; ++j) {
            double t = (nB == 1) ? 0.5 : (double)j / (nB - 1);
            Vec3 B = h00(t) * PL + h10(t) * PR + h01(t) * TP + h11(t) * TQ;
            s.meshVerts[i * nB + j] = B;
            double u = uL + gapU * t;
            Vec3 S = surf.evaluate(u, v);
            addSample(e, (B - S).norm());
        }
    }

    for (int i = 0; i < nA - 1; ++i)
        for (int j = 0; j < nB - 1; ++j) {
            int a = i * nB + j;
            s.meshFaces.push_back({a, a + 1, a + nB});
            s.meshFaces.push_back({a + 1, a + nB + 1, a + nB});
        }

    finalize(e);
    s.maxError = e.maxError;
    s.rmsError = e.rmsError;
    return s;
}

// 水平缝（v = vEdges[r]）：上格 (r-1,c) 与下格 (r,c) 之间
BlendStrip buildHorizontalStrip(const SurfaceWrapper& surf, const GridResult& gr,
                                const std::vector<double>& f,
                                int r, int c, const BlendConfig& cfg) {
    BlendStrip s;
    s.r = r; s.c = c; s.isVertical = false;

    const GridCell& T = gr.cells[(r - 1) * gr.nCols + c];
    const GridCell& B = gr.cells[r * gr.nCols + c];
    TrimmedBounds tT = trimOf(gr, T, f[(r - 1) * gr.nCols + c]);
    TrimmedBounds tB = trimOf(gr, B, f[r * gr.nCols + c]);

    double vT = tT.v1;   // 上格下边界
    double vB = tB.v0;   // 下格上边界
    double gapV = std::max(vB - vT, 1e-12);
    double uMin = gr.uEdges[c];
    double uMax = gr.uEdges[c + 1];

    int nA = cfg.nAlong, nB = cfg.nAcross;
    s.meshVerts.resize(nA * nB);
    ErrorStat e;

    for (int i = 0; i < nA; ++i) {
        double su = (nA == 1) ? 0.5 : (double)i / (nA - 1);
        double u = uMin + (uMax - uMin) * su;
        Vec3 PT, dUT, dVT;
        evalRuledCell(T, u, vT, PT, dUT, dVT);
        Vec3 PB, dUB, dVB;
        evalRuledCell(B, u, vB, PB, dUB, dVB);
        Vec3 TP = dVT * gapV;
        Vec3 TQ = dVB * gapV;

        for (int j = 0; j < nB; ++j) {
            double t = (nB == 1) ? 0.5 : (double)j / (nB - 1);
            Vec3 P = h00(t) * PT + h10(t) * PB + h01(t) * TP + h11(t) * TQ;
            s.meshVerts[i * nB + j] = P;
            double v = vT + gapV * t;
            Vec3 S = surf.evaluate(u, v);
            addSample(e, (P - S).norm());
        }
    }

    for (int i = 0; i < nA - 1; ++i)
        for (int j = 0; j < nB - 1; ++j) {
            int a = i * nB + j;
            s.meshFaces.push_back({a, a + 1, a + nB});
            s.meshFaces.push_back({a + 1, a + nB + 1, a + nB});
        }

    finalize(e);
    s.maxError = e.maxError;
    s.rmsError = e.rmsError;
    return s;
}

double stripErrorMax(const BlendStrip& s, bool useMaxError) {
    return useMaxError ? s.maxError : s.rmsError;
}

double cellErrorVal(const ErrorStat& e, bool useMaxError) {
    return useMaxError ? e.maxError : e.rmsError;
}

} // anon

BlendResult buildBlend(const SurfaceWrapper& surf, const GridResult& gr,
                       const std::vector<double>& f, const BlendConfig& cfg)
{
    BlendResult br;
    double totalSumSq = 0.0;
    int totalCnt = 0;

    br.cells.reserve(gr.cells.size());
    for (const auto& cell : gr.cells) {
        int idx = cell.row * gr.nCols + cell.col;
        double fi = (idx < (int)f.size()) ? f[idx] : 0.0;
        TrimmedCell tc;
        tc.row = cell.row;
        tc.col = cell.col;
        tc.f = fi;
        TrimmedBounds t = trimOf(gr, cell, fi);
        tc.u0 = t.u0; tc.u1 = t.u1; tc.v0 = t.v0; tc.v1 = t.v1;
        ErrorStat e = cellError(surf, gr, cell, fi, cfg);
        tc.maxError = e.maxError;
        tc.rmsError = e.rmsError;
        totalSumSq += e.sumSq;
        totalCnt += e.count;
        meshTrimmedCell(gr, cell, fi, cfg.nMeshRes, tc.meshVerts, tc.meshFaces);
        br.cells.push_back(std::move(tc));
    }

    for (int r = 0; r < gr.nRows; ++r)
        for (int c = 1; c < gr.nCols; ++c) {
            BlendStrip s = buildVerticalStrip(surf, gr, f, r, c, cfg);
            totalSumSq += s.rmsError * s.rmsError * (int)s.meshVerts.size();
            totalCnt += (int)s.meshVerts.size();
            br.strips.push_back(std::move(s));
        }
    for (int r = 1; r < gr.nRows; ++r)
        for (int c = 0; c < gr.nCols; ++c) {
            BlendStrip s = buildHorizontalStrip(surf, gr, f, r, c, cfg);
            totalSumSq += s.rmsError * s.rmsError * (int)s.meshVerts.size();
            totalCnt += (int)s.meshVerts.size();
            br.strips.push_back(std::move(s));
        }

    br.totalMaxError = 0.0;
    for (const auto& tc : br.cells) br.totalMaxError = std::max(br.totalMaxError, tc.maxError);
    for (const auto& s : br.strips) br.totalMaxError = std::max(br.totalMaxError, s.maxError);
    br.totalRmsError = totalCnt ? std::sqrt(totalSumSq / totalCnt) : 0.0;
    return br;
}

std::vector<double> optimizeTrim(const SurfaceWrapper& surf, const GridResult& gr,
                                 const BlendConfig& cfg)
{
    int nCells = (int)gr.cells.size();
    std::vector<double> f(nCells, 0.0);

    auto localObjective = [&](int i) -> double {
        const GridCell& cell = gr.cells[i];
        double e = cellErrorVal(cellError(surf, gr, cell, f[i], cfg), cfg.useMaxError);
        int r = cell.row, c = cell.col;
        if (c > 0)
            e = std::max(e, stripErrorMax(buildVerticalStrip(surf, gr, f, r, c, cfg), cfg.useMaxError));
        if (c + 1 < gr.nCols)
            e = std::max(e, stripErrorMax(buildVerticalStrip(surf, gr, f, r, c + 1, cfg), cfg.useMaxError));
        if (r > 0)
            e = std::max(e, stripErrorMax(buildHorizontalStrip(surf, gr, f, r, c, cfg), cfg.useMaxError));
        if (r + 1 < gr.nRows)
            e = std::max(e, stripErrorMax(buildHorizontalStrip(surf, gr, f, r + 1, c, cfg), cfg.useMaxError));
        return e;
    };

    for (int pass = 0; pass < cfg.maxPasses; ++pass) {
        bool changed = false;
        for (int i = 0; i < nCells; ++i) {
            double orig = f[i];
            double bestF = orig;
            double bestE = localObjective(i);
            for (int k = 1; k <= cfg.nSearch; ++k) {
                double cand = cfg.fMax * k / cfg.nSearch;
                f[i] = cand;
                double e = localObjective(i);
                if (e < bestE) { bestE = e; bestF = cand; }
            }
            f[i] = bestF;
            if (std::abs(bestF - orig) > 1e-6) changed = true;
        }
        if (!changed) break;
    }
    return f;
}

bool exportBlendOBJs(const std::string& outDir, const std::string& prefix,
                     const BlendResult& br)
{
    for (const auto& tc : br.cells) {
        std::string op = outDir + "/" + prefix + "_trim_" + std::to_string(tc.row)
                         + "_" + std::to_string(tc.col) + ".obj";
        exportOBJ(op, tc.meshVerts, tc.meshFaces);
    }
    for (size_t si = 0; si < br.strips.size(); ++si) {
        const auto& s = br.strips[si];
        std::string op = outDir + "/" + prefix + "_strip" + std::to_string(si) + ".obj";
        exportOBJ(op, s.meshVerts, s.meshFaces);
    }
    return true;
}

} // namespace simple
