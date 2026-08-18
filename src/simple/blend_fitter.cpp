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
    double gapRaw = uR - uL;
    if (gapRaw < 1e-9) { s.valid = false; return s; }   // 两侧 f=0，退化
    double gapU = gapRaw;
    double vMin = std::max(tL.v0, tR.v0);   // 缩到两侧 trimmed v 范围交集
    double vMax = std::min(tL.v1, tR.v1);
    if (vMax - vMin < 1e-12) { s.valid = false; return s; }

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
    double gapRaw = vB - vT;
    if (gapRaw < 1e-9) { s.valid = false; return s; }   // 两侧 f=0，退化
    double gapV = gapRaw;
    double uMin = std::max(tT.u0, tB.u0);   // 缩到两侧 trimmed u 范围交集
    double uMax = std::min(tT.u1, tB.u1);
    if (uMax - uMin < 1e-12) { s.valid = false; return s; }

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

// 角点过渡 patch：内部网格顶点 (c,r) 处的双三次 Hermite (Ferguson) 面片。
// 插值 4 角点位置 + 8 个一阶偏导（扭矢置 0），与周围 4 条带同源：
// 四条边即条带端边的 Hermite 曲线（G0 严丝合缝），跨界切线在角点处严格、
// 沿边近似 G1。
BlendCorner buildCornerPatch(const SurfaceWrapper& surf, const GridResult& gr,
                             const std::vector<double>& f,
                             int r, int c, const BlendConfig& cfg) {
    BlendCorner cp;
    cp.r = r; cp.c = c;

    const GridCell& TL = gr.cells[(r - 1) * gr.nCols + (c - 1)];
    const GridCell& TR = gr.cells[(r - 1) * gr.nCols + c];
    const GridCell& BL = gr.cells[r * gr.nCols + (c - 1)];
    const GridCell& BR = gr.cells[r * gr.nCols + c];
    TrimmedBounds tTL = trimOf(gr, TL, f[(r - 1) * gr.nCols + (c - 1)]);
    TrimmedBounds tTR = trimOf(gr, TR, f[(r - 1) * gr.nCols + c]);
    TrimmedBounds tBL = trimOf(gr, BL, f[r * gr.nCols + (c - 1)]);
    TrimmedBounds tBR = trimOf(gr, BR, f[r * gr.nCols + c]);

    double uL = std::max(tTL.u1, tBL.u1);
    double uR = std::min(tTR.u0, tBR.u0);
    double vB = std::max(tTL.v1, tTR.v1);
    double vT = std::min(tBL.v0, tBR.v0);
    if (uR - uL < 1e-12 || vT - vB < 1e-12) return cp;   // 退化角点

    int n = cfg.nMeshRes + 1;
    cp.meshVerts.resize(n * n);

    // 双三次 Hermite (Ferguson) 角点面片：
    // 插值 4 角点位置 + 8 个一阶偏导（扭矢置 0），四条边即为条带端边的
    // Hermite 曲线（G0 严丝合缝），跨界切线由角点偏导插值（角点处严格、
    // 沿边近似 G1）。
    double du = uR - uL, dv = vT - vB;

    Vec3 P00, P10, P01, P11;                 // 角点位置 (s,t)=(i,j)
    Vec3 Su00, Su10, Su01, Su11;             // ∂/∂s = ∂/∂u * du
    Vec3 Sv00, Sv10, Sv01, Sv11;             // ∂/∂t = ∂/∂v * dv
    {
        Vec3 dU, dV;
        evalRuledCell(TL, uL, vB, P00, dU, dV); Su00 = dU * du; Sv00 = dV * dv;
        evalRuledCell(TR, uR, vB, P10, dU, dV); Su10 = dU * du; Sv10 = dV * dv;
        evalRuledCell(BL, uL, vT, P01, dU, dV); Su01 = dU * du; Sv01 = dV * dv;
        evalRuledCell(BR, uR, vT, P11, dU, dV); Su11 = dU * du; Sv11 = dV * dv;
    }

    ErrorStat e;
    for (int i = 0; i < n; ++i) {
        double s = (n == 1) ? 0.5 : (double)i / (n - 1);
        for (int j = 0; j < n; ++j) {
            double t = (n == 1) ? 0.5 : (double)j / (n - 1);
            Vec3 P =
                  h00(s) * (h00(t) * P00 + h10(t) * P01 + h01(t) * Sv00 + h11(t) * Sv01)
                + h10(s) * (h00(t) * P10 + h10(t) * P11 + h01(t) * Sv10 + h11(t) * Sv11)
                + h01(s) * (h00(t) * Su00 + h10(t) * Su01)
                + h11(s) * (h00(t) * Su10 + h10(t) * Su11);
            cp.meshVerts[i * n + j] = P;
            double u = uL + s * (uR - uL);
            double v = vB + t * (vT - vB);
            Vec3 S = surf.evaluate(u, v);
            addSample(e, (P - S).norm());
        }
    }
    for (int i = 0; i < n - 1; ++i)
        for (int j = 0; j < n - 1; ++j) {
            int a = i * n + j;
            cp.meshFaces.push_back({a, a + 1, a + n});
            cp.meshFaces.push_back({a + 1, a + n + 1, a + n});
        }

    finalize(e);
    cp.maxError = e.maxError;
    cp.rmsError = e.rmsError;
    return cp;
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
            if (!s.valid) continue;
            totalSumSq += s.rmsError * s.rmsError * (int)s.meshVerts.size();
            totalCnt += (int)s.meshVerts.size();
            br.strips.push_back(std::move(s));
        }
    for (int r = 1; r < gr.nRows; ++r)
        for (int c = 0; c < gr.nCols; ++c) {
            BlendStrip s = buildHorizontalStrip(surf, gr, f, r, c, cfg);
            if (!s.valid) continue;
            totalSumSq += s.rmsError * s.rmsError * (int)s.meshVerts.size();
            totalCnt += (int)s.meshVerts.size();
            br.strips.push_back(std::move(s));
        }

    for (int r = 1; r < gr.nRows; ++r)
        for (int c = 1; c < gr.nCols; ++c) {
            BlendCorner cp = buildCornerPatch(surf, gr, f, r, c, cfg);
            if (cp.meshVerts.empty()) continue;
            totalSumSq += cp.rmsError * cp.rmsError * (int)cp.meshVerts.size();
            totalCnt += (int)cp.meshVerts.size();
            br.corners.push_back(std::move(cp));
        }

    br.totalMaxError = 0.0;
    for (const auto& tc : br.cells) br.totalMaxError = std::max(br.totalMaxError, tc.maxError);
    for (const auto& s : br.strips) br.totalMaxError = std::max(br.totalMaxError, s.maxError);
    for (const auto& cp : br.corners) br.totalMaxError = std::max(br.totalMaxError, cp.maxError);
    br.totalRmsError = totalCnt ? std::sqrt(totalSumSq / totalCnt) : 0.0;
    return br;
}

// 统一内缩量 f：所有格用同一个 f，保证过渡网格（条带+角点）无台阶、可 G0 拼接。
// 对 f ∈ [0, fMax] 做一维搜索，最小化复合面总误差。
std::vector<double> optimizeTrim(const SurfaceWrapper& surf, const GridResult& gr,
                                 const BlendConfig& cfg)
{
    int nCells = (int)gr.cells.size();
    double bestF = 0.0;
    double bestE = 1e30;

    for (int k = 0; k <= cfg.nSearch; ++k) {
        double cand = cfg.fMax * k / cfg.nSearch;
        std::vector<double> f(nCells, cand);
        BlendResult br = buildBlend(surf, gr, f, cfg);
        double e = cfg.useMaxError ? br.totalMaxError : br.totalRmsError;
        if (e < bestE) { bestE = e; bestF = cand; }
    }
    return std::vector<double>(nCells, bestF);
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
    for (const auto& cp : br.corners) {
        std::string op = outDir + "/" + prefix + "_corner_" + std::to_string(cp.r)
                         + "_" + std::to_string(cp.c) + ".obj";
        exportOBJ(op, cp.meshVerts, cp.meshFaces);
    }
    return true;
}

} // namespace simple
