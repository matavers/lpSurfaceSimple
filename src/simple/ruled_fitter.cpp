#include "simple/ruled_fitter.hpp"

#include <fstream>
#include <iomanip>
#include <limits>

#include <Geom_BSplineCurve.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAbs_Shape.hxx>
#include <NCollection_Array1.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <gp_Pnt.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_HArray1OfReal.hxx>

namespace simple {

// 计算两条准线采样点的共同弦长参数（归一化到 [0,1]），保证 C0/C1 共享同一参数化。
static std::vector<double> chordLengthParams(const Vec3Arr& c0, const Vec3Arr& c1) {
    int n = static_cast<int>(std::min(c0.size(), c1.size()));
    std::vector<double> t(n, 0.0);
    if (n == 0) return t;
    for (int i = 1; i < n; ++i) {
        double d = 0.5 * ((c0[i] - c0[i - 1]).norm() + (c1[i] - c1[i - 1]).norm());
        t[i] = t[i - 1] + std::max(d, 1e-9);
    }
    double L = t.back();
    if (L > 0.0)
        for (double& x : t) x /= L;
    return t;
}

// 用给定参数把采样点拟合成 3 次 B 样条曲线（用于导出节点向量与控制点）。
static Handle(Geom_BSplineCurve) fitBSpline(const Vec3Arr& pts,
                                            const std::vector<double>& params) {
    int n = static_cast<int>(pts.size());
    if (n < 2) return Handle(Geom_BSplineCurve)();

    NCollection_Array1<gp_Pnt> points(1, n);
    NCollection_Array1<double> pars(1, n);
    for (int i = 0; i < n; ++i) {
        points.SetValue(i + 1, gp_Pnt(pts[i].x(), pts[i].y(), pts[i].z()));
        pars.SetValue(i + 1, params[i]);
    }
    GeomAPI_PointsToBSpline approx(points, pars, 3, 3, GeomAbs_C2, 1e-6);
    if (!approx.IsDone()) return Handle(Geom_BSplineCurve)();
    Handle(Geom_BSplineCurve) curve = approx.Curve();
    return curve;
}

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

bool exportCurveBSplineTXT(const std::string& path,
                           const Handle(Geom_BSplineCurve)& curveC0,
                           const Handle(Geom_BSplineCurve)& curveC1,
                           int nSamples)
{
    std::ofstream out(path);
    if (!out || curveC0.IsNull() || curveC1.IsNull()) return false;

    auto writeCurve = [&](const char* name, const Handle(Geom_BSplineCurve)& c) {
        out << "[" << name << "]\n";
        out << "degree = " << c->Degree() << "\n";
        out << "rational = " << (c->IsRational() ? "yes" : "no") << "\n";
        out << "nbPoles = " << c->NbPoles() << "\n";
        out << "poles (x y z w):\n";
        for (int i = 1; i <= c->NbPoles(); ++i) {
            gp_Pnt p = c->Pole(i);
            out << "  " << std::fixed << std::setprecision(6)
                << p.X() << " " << p.Y() << " " << p.Z() << " " << c->Weight(i) << "\n";
        }
        out << "nbKnots = " << c->NbKnots() << "\n";
        out << "knots:";
        for (int i = 1; i <= c->NbKnots(); ++i)
            out << " " << std::setprecision(8) << c->Knot(i);
        out << "\nmultiplicities:";
        for (int i = 1; i <= c->NbKnots(); ++i)
            out << " " << c->Multiplicity(i);
        out << "\n";
    };

    out << "nSamples = " << nSamples << "\n\n";
    writeCurve("C0", curveC0);
    out << "\n";
    writeCurve("C1", curveC1);
    out << "\n[mapping]\n";
    out << "description = identity, C0(u) and C1(u) share the same parameterization (i-th point/control pair forms a ruling)\n";
    return true;
}

bool exportDirectrixBSplineTXT(const std::string& path,
                               const Vec3Arr& c0Samples,
                               const Vec3Arr& c1Samples)
{
    if (c0Samples.size() < 2 || c1Samples.size() < 2) return false;

    std::vector<double> params = chordLengthParams(c0Samples, c1Samples);
    Handle(Geom_BSplineCurve) c0 = fitBSpline(c0Samples, params);
    Handle(Geom_BSplineCurve) c1 = fitBSpline(c1Samples, params);
    if (c0.IsNull() || c1.IsNull()) return false;

    return exportCurveBSplineTXT(path, c0, c1, static_cast<int>(c0Samples.size()));
}

bool exportErrorsCSV(const std::string& path,
                     const std::vector<RuledResult>& results)
{
    std::ofstream out(path);
    if (!out) return false;
    out << "surface,version,segment,splitDir,directrixDir,maxError,rmsError\n";
    out << std::fixed << std::setprecision(6);
    for (const auto& res : results) {
        for (const auto& seg : res.segments) {
            out << res.name << "," << res.version << ","
                << seg.segmentIndex << ","
                << (int)res.splitDir << ","
                << (int)seg.directrixDir << ","
                << seg.maxError << "," << seg.rmsError << "\n";
        }
    }
    return true;
}

bool exportMetaJSON(const std::string& path,
                    const std::string& file1, const std::string& file2,
                    const std::vector<RuledResult>& results)
{
    auto jsonSafe = [](std::string s) {
        for (char& c : s) { if (c == '\\') c = '/'; if (c == '"') c = '\''; }
        return s;
    };
    std::ofstream out(path);
    if (!out) return false;
    out << "{\n";
    out << "  \"files\": [\"" << jsonSafe(file1) << "\", \"" << jsonSafe(file2) << "\"],\n";
    out << "  \"surfaces\": [\n";
    for (size_t i = 0; i < results.size(); ++i) {
        out << "    {\n";
        out << "      \"name\": \"" << results[i].name << "\",\n";
        out << "      \"version\": " << results[i].version << ",\n";
        out << "      \"splitDir\": " << (int)results[i].splitDir << ",\n";
        out << "      \"numSegments\": " << results[i].segments.size() << ",\n";
        out << "      \"segments\": [\n";
        for (size_t j = 0; j < results[i].segments.size(); ++j) {
            out << "        { \"index\": " << results[i].segments[j].segmentIndex
                << ", \"directrixDir\": " << (int)results[i].segments[j].directrixDir
                << ", \"maxError\": " << results[i].segments[j].maxError
                << ", \"rmsError\": " << results[i].segments[j].rmsError << " }";
            if (j + 1 < results[i].segments.size()) out << ",";
            out << "\n";
        }
        out << "      ]\n";
        out << "    }";
        if (i + 1 < results.size()) out << ",";
        out << "\n";
    }
    out << "  ]\n}\n";
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

    r.ruledMeshVerts = generateRuledMesh(
        r.curveC0Samples, r.curveC1Samples,
        nCurveSamples, nAcrossSamples, r.ruledMeshFaces);

    auto [mx, rms] = computeError(surf, r, u0, u1, v0, v1, nCurveSamples, nAcrossSamples);
    r.maxError = mx; r.rmsError = rms;

    return r;
}

RuledCellFit fitCellRuledAuto(const SurfaceWrapper& surf,
                              double u0, double u1, double v0, double v1,
                              int nUSamples, int nVSamples,
                              int nRibs, double lambda)
{
    RuledCellFit fv = fitCellRuled(surf, u0, u1, v0, v1, ParamDir::V,
                                   nUSamples, nVSamples, nRibs, lambda);
    RuledCellFit fu = fitCellRuled(surf, u0, u1, v0, v1, ParamDir::U,
                                   nUSamples, nVSamples, nRibs, lambda);
    return (fv.maxError <= fu.maxError) ? fv : fu;
}

RuledResult fitRuledSegments(const SurfaceWrapper& surf,
                              int numSegments,
                              ParamDir splitDir,
                              const std::vector<ParamDir>& directrixDirs,
                              int nUSamples,
                              int nVSamples,
                              int version,
                              const std::string& name)
{
    RuledResult result;
    result.name = name;
    result.version = version;
    result.splitDir = splitDir;

    auto [uFull0, uFull1] = surf.paramDomainU();
    auto [vFull0, vFull1] = surf.paramDomainV();

    for (int seg = 0; seg < numSegments; ++seg) {
        ParamDir ddir = (seg < (int)directrixDirs.size()) ? directrixDirs[seg] : ParamDir::V;

        double uSeg0, uSeg1, vSeg0, vSeg1;
        if (splitDir == ParamDir::V) {
            uSeg0 = uFull0; uSeg1 = uFull1;
            vSeg0 = vFull0 + (vFull1 - vFull0) * seg / numSegments;
            vSeg1 = vFull0 + (vFull1 - vFull0) * (seg + 1) / numSegments;
        } else {
            uSeg0 = uFull0 + (uFull1 - uFull0) * seg / numSegments;
            uSeg1 = uFull0 + (uFull1 - uFull0) * (seg + 1) / numSegments;
            vSeg0 = vFull0; vSeg1 = vFull1;
        }

        RuledCellFit cell = fitCellRuled(surf, uSeg0, uSeg1, vSeg0, vSeg1, ddir,
                                         nUSamples, nVSamples, 20, 1.0);

        RuledSegment rseg;
        rseg.segmentIndex = seg;
        rseg.directrixDir = ddir;
        rseg.curveC0 = cell.curveC0;
        rseg.curveC1 = cell.curveC1;
        rseg.curveC0Samples = cell.curveC0Samples;
        rseg.curveC1Samples = cell.curveC1Samples;
        rseg.ruledMeshVerts = cell.ruledMeshVerts;
        rseg.ruledMeshFaces = cell.ruledMeshFaces;
        rseg.maxError = cell.maxError;
        rseg.rmsError = cell.rmsError;

        result.segments.push_back(rseg);
    }

    return result;
}

} // namespace simple
