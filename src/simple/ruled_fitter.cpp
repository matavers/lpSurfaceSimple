#include "simple/ruled_fitter.hpp"

#include <fstream>
#include <iomanip>

#include <Geom_BSplineCurve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <gp_Pnt.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColStd_Array1OfReal.hxx>

namespace simple {

static double evalU(double u0, double u1, double t, bool wrap) {
    if (!wrap || u0 <= u1)
        return u0 + (u1 - u0) * t;
    double rng = (1.0 - u0) + u1;
    double u = u0 + rng * t;
    if (u > 1.0) u -= 1.0;
    return u;
}

bool exportOBJ(const std::string& path,
               const Vec3Arr& verts, const FaceArr& faces)
{
    std::ofstream out(path);
    if (!out) return false;
    out << std::fixed << std::setprecision(6);
    for (const auto& v : verts)
        out << "v " << v.x() << " " << v.y() << " " << v.z() << "\n";
    for (const auto& f : faces)
        out << "f " << (f[0] + 1) << " " << (f[1] + 1) << " " << (f[2] + 1) << "\n";
    return true;
}

bool exportDirectrixTXT(const std::string& path,
                        const Vec3Arr& c0Samples, const Vec3Arr& c1Samples)
{
    std::ofstream out(path);
    if (!out) return false;
    out << std::fixed << std::setprecision(6);

    auto writePoints = [&](const char* name, const Vec3Arr& pts) {
        out << "[" << name << "]\n";
        out << "n = " << pts.size() << "\n";
        for (const auto& p : pts)
            out << p.x() << " " << p.y() << " " << p.z() << "\n";
    };

    writePoints("C0", c0Samples);
    out << "\n";
    writePoints("C1", c1Samples);
    out << "\n[mapping]\n";
    out << "description = identity, C0 and C1 share the same parameterization (i-th sample pair forms a ruling)\n";
    return true;
}

Vec3Arr generateRuledMesh(const Vec3Arr& c0Samples,
                           const Vec3Arr& c1Samples,
                           int nAlong, int nAcross,
                           FaceArr& faces)
{
    int n = static_cast<int>(c0Samples.size());
    Vec3Arr verts(n * nAcross);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < nAcross; ++j) {
            double t = j / (nAcross - 1.0);
            int idx = i * nAcross + j;
            verts[idx] = c0Samples[i] * (1.0 - t) + c1Samples[i] * t;
        }
    }

    faces.clear();
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < nAcross - 1; ++j) {
            int a = i * nAcross + j;
            faces.push_back({a, a + 1, a + nAcross});
            faces.push_back({a + 1, a + nAcross + 1, a + nAcross});
        }
    }

    return verts;
}

static std::pair<double, double> computeError(const SurfaceWrapper& surf,
                                               const RuledCellFit& seg,
                                               double uSeg0, double uSeg1,
                                               double vSeg0, double vSeg1,
                                               int nAlong, int nAcross)
{
    double totalDist = 0.0;
    double maxDist = 0.0;
    int count = 0;

    bool wrap = surf.isWrapU();
    for (int i = 0; i < nAlong; ++i) {
        for (int j = 0; j < nAcross; ++j) {
            double u, v;
            if (seg.fitDir == ParamDir::V) {
                u = evalU(uSeg0, uSeg1, i / (nAlong - 1.0), wrap);
                v = vSeg0 + (vSeg1 - vSeg0) * j / (nAcross - 1.0);
            } else {
                v = vSeg0 + (vSeg1 - vSeg0) * i / (nAlong - 1.0);
                u = evalU(uSeg0, uSeg1, j / (nAcross - 1.0), wrap);
            }
            Vec3 surfPt = surf.evaluate(u, v);
            double bestDist = std::numeric_limits<double>::max();
            int nR = static_cast<int>(seg.curveC0Samples.size());
            for (int ri = 0; ri < nR; ++ri) {
                for (int rj = 0; rj < 10; ++rj) {
                    double t = rj / 9.0;
                    Vec3 ruledPt = seg.curveC0Samples[ri] * (1.0 - t) +
                                   seg.curveC1Samples[ri] * t;
                    double d = (surfPt - ruledPt).norm();
                    if (d < bestDist) bestDist = d;
                }
            }
            totalDist += bestDist;
            if (bestDist > maxDist) maxDist = bestDist;
            ++count;
        }
    }

    double rms = count > 0 ? std::sqrt(totalDist / count) : 0.0;
    return {maxDist, rms};
}

void optimizeDirectrices(
    const SurfaceWrapper& surf,
    double uSeg0, double uSeg1,
    double vSeg0, double vSeg1,
    ParamDir directrixDir,
    Vec3Arr& c0Samples,
    Vec3Arr& c1Samples,
    int nRibs, double lambda)
{
    int n = static_cast<int>(c0Samples.size());
    Vec3Arr origC0 = c0Samples;
    Vec3Arr origC1 = c1Samples;

    bool wrap = surf.isWrapU();

    for (int i = 0; i < n; ++i) {
        double sumW = 0, sumZ = 0, sumWW = 0, sumZZ = 0, sumWZ = 0;
        Vec3 sumWR(0,0,0), sumZR(0,0,0);

        for (int j = 0; j < nRibs; ++j) {
            double t = j / (nRibs - 1.0);
            double w = 1.0 - t;
            double z = t;

            Vec3 R;
            if (directrixDir == ParamDir::V) {
                double u = evalU(uSeg0, uSeg1, i / (n - 1.0), wrap);
                double v = vSeg0 + (vSeg1 - vSeg0) * t;
                R = surf.evaluate(u, v);
            } else {
                double v = vSeg0 + (vSeg1 - vSeg0) * i / (n - 1.0);
                double u = evalU(uSeg0, uSeg1, t, wrap);
                R = surf.evaluate(u, v);
            }

            sumW += w;
            sumZ += z;
            sumWW += w * w;
            sumZZ += z * z;
            sumWZ += w * z;
            sumWR += w * R;
            sumZR += z * R;
        }

        double a = sumWW + lambda;
        double b = sumWZ;
        double c = sumZZ + lambda;
        double det = a * c - b * b;
        if (std::abs(det) < 1e-12) continue;

        Vec3 X0 = origC0[i];
        Vec3 Y0 = origC1[i];
        Vec3 p = sumWR + lambda * X0;
        Vec3 q = sumZR + lambda * Y0;
        double invDet = 1.0 / det;
        c0Samples[i] = Vec3(
            (c * p.x() - b * q.x()) * invDet,
            (c * p.y() - b * q.y()) * invDet,
            (c * p.z() - b * q.z()) * invDet);
        c1Samples[i] = Vec3(
            (a * q.x() - b * p.x()) * invDet,
            (a * q.y() - b * p.y()) * invDet,
            (a * q.z() - b * p.z()) * invDet);
    }
}

// 逐母线法向扭转角(度)：直纹面可展 ⟺ 切平面(法向)沿母线恒定 ⟺ β≈0。
// 与 compute_machining.py 的 ruling_twists 口径一致。
static VecDArr computeRuledTwists(const Vec3Arr& c0, const Vec3Arr& c1) {
    int n = static_cast<int>(std::min(c0.size(), c1.size()));
    VecDArr twists(n, 0.0);
    if (n < 2) return twists;

    auto tangent = [](const Vec3Arr& C, int i) -> Vec3 {
        int m = static_cast<int>(C.size());
        if (m < 2) return Vec3(0, 0, 0);
        Vec3 d;
        if (i == 0) d = C[1] - C[0];
        else if (i == m - 1) d = C[m - 1] - C[m - 2];
        else d = C[i + 1] - C[i - 1];
        double len = d.norm();
        return len < 1e-12 ? Vec3(0, 0, 0) : d / len;
    };

    for (int i = 0; i < n; ++i) {
        Vec3 r = c1[i] - c0[i];
        if (r.norm() < 1e-12) continue;
        Vec3 t0 = tangent(c0, i);
        Vec3 t1 = tangent(c1, i);
        Vec3 n0 = t0.cross(r);
        Vec3 n1 = t1.cross(r);
        double m0 = n0.norm(), m1 = n1.norm();
        if (m0 < 1e-12 || m1 < 1e-12) continue;
        double d = clamp(n0.dot(n1) / (m0 * m1), -1.0, 1.0);
        twists[i] = radToDeg(std::acos(d));
    }
    return twists;
}

static double computeRuledTwist(const Vec3Arr& c0, const Vec3Arr& c1) {
    VecDArr tw = computeRuledTwists(c0, c1);
    double mx = 0.0;
    for (double t : tw) mx = std::max(mx, t);
    return mx;
}

RuledCellFit fitCellRuled(const SurfaceWrapper& surf,
                          double u0, double u1, double v0, double v1,
                          ParamDir dir,
                          int nUSamples, int nVSamples,
                          int nRibs, double lambda)
{
    RuledCellFit r;
    r.fitDir = dir;
    r.maxError = 0.0;
    r.rmsError = 0.0;
    r.twist = 0.0;

    bool wrap = surf.isWrapU();

    int nCurveSamples = (dir == ParamDir::V) ? nUSamples : std::max(nUSamples, nVSamples);
    int nAcrossSamples = (dir == ParamDir::V) ? nVSamples : std::min(nUSamples, nVSamples);

    if (dir == ParamDir::V) {
        r.curveC0 = surf.extractIsoCurveV(v0);
        r.curveC1 = surf.extractIsoCurveV(v1);
        if (wrap && u0 > u1) {
            r.curveC0Samples = sampleCurveRangeWrap(r.curveC0, u0, u1, nCurveSamples);
            r.curveC1Samples = sampleCurveRangeWrap(r.curveC1, u0, u1, nCurveSamples);
        } else {
            r.curveC0Samples = sampleCurveRange(r.curveC0, u0, u1, nCurveSamples);
            r.curveC1Samples = sampleCurveRange(r.curveC1, u0, u1, nCurveSamples);
        }
    } else {
        r.curveC0 = surf.extractIsoCurveU(u0);
        r.curveC1 = surf.extractIsoCurveU(u1);
        r.curveC0Samples = sampleCurveRange(r.curveC0, v0, v1, nCurveSamples);
        r.curveC1Samples = sampleCurveRange(r.curveC1, v0, v1, nCurveSamples);
    }

    if (r.curveC0Samples.size() < 2 || r.curveC1Samples.size() < 2) return r;

    optimizeDirectrices(surf, u0, u1, v0, v1, dir,
        r.curveC0Samples, r.curveC1Samples, nRibs, lambda);

    r.twistPerRuling = computeRuledTwists(r.curveC0Samples, r.curveC1Samples);
    r.twist = 0.0;
    for (double t : r.twistPerRuling) r.twist = std::max(r.twist, t);

    r.ruledMeshVerts = generateRuledMesh(
        r.curveC0Samples, r.curveC1Samples,
        nCurveSamples, nAcrossSamples, r.ruledMeshFaces);

    auto [mx, rms] = computeError(surf, r, u0, u1, v0, v1, nCurveSamples, nAcrossSamples);
    r.maxError = mx; r.rmsError = rms;

    return r;
}

} // namespace simple
