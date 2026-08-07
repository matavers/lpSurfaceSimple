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
    std::ofstream out(path);
    if (!out) return false;
    out << "{\n";
    out << "  \"files\": [\"" << file1 << "\", \"" << file2 << "\"],\n";
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

bool exportCurveParamsTXT(const std::string& path,
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
    out << "description = identity, C0(u) and C1(u) share same U-parameterization from surface iso-curves\n";
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

std::pair<double, double> computeError(const SurfaceWrapper& surf,
                                        const RuledSegment& seg,
                                        double uSeg0, double uSeg1,
                                        double vSeg0, double vSeg1,
                                        int nAlong, int nAcross)
{
    double totalDist = 0.0;
    double maxDist = 0.0;
    int count = 0;

    bool wrap = surf.isWrapU();
    if (seg.directrixDir == ParamDir::V) {
        for (int i = 0; i < nAlong; ++i) {
            for (int j = 0; j < nAcross; ++j) {
                double u = evalU(uSeg0, uSeg1, i / (nAlong - 1.0), wrap);
                double v = vSeg0 + (vSeg1 - vSeg0) * j / (nAcross - 1.0);
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
    } else {
        for (int i = 0; i < nAlong; ++i) {
            for (int j = 0; j < nAcross; ++j) {
                double v = vSeg0 + (vSeg1 - vSeg0) * i / (nAlong - 1.0);
                double u = evalU(uSeg0, uSeg1, j / (nAcross - 1.0), wrap);
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
    bool wrap = surf.isWrapU();

    for (int seg = 0; seg < numSegments; ++seg) {
        RuledSegment rseg;
        rseg.segmentIndex = seg;
        ParamDir ddir = (seg < (int)directrixDirs.size()) ? directrixDirs[seg] : ParamDir::V;
        rseg.directrixDir = ddir;

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

        int nCurveSamples = (ddir == ParamDir::V) ? nUSamples : std::max(nUSamples, nVSamples);
        int nAcrossSamples = (ddir == ParamDir::V) ? nVSamples : std::min(nUSamples, nVSamples);
        int nRibs = 20;

        if (ddir == ParamDir::V) {
            rseg.curveC0 = surf.extractIsoCurveV(vSeg0);
            rseg.curveC1 = surf.extractIsoCurveV(vSeg1);
            if (wrap && uSeg0 > uSeg1) {
                rseg.curveC0Samples = sampleCurveRangeWrap(rseg.curveC0, uSeg0, uSeg1, nCurveSamples);
                rseg.curveC1Samples = sampleCurveRangeWrap(rseg.curveC1, uSeg0, uSeg1, nCurveSamples);
            } else {
                rseg.curveC0Samples = sampleCurveRange(rseg.curveC0, uSeg0, uSeg1, nCurveSamples);
                rseg.curveC1Samples = sampleCurveRange(rseg.curveC1, uSeg0, uSeg1, nCurveSamples);
            }
        } else {
            rseg.curveC0 = surf.extractIsoCurveU(uSeg0);
            rseg.curveC1 = surf.extractIsoCurveU(uSeg1);
            rseg.curveC0Samples = sampleCurveRange(rseg.curveC0, vSeg0, vSeg1, nCurveSamples);
            rseg.curveC1Samples = sampleCurveRange(rseg.curveC1, vSeg0, vSeg1, nCurveSamples);
        }

        optimizeDirectrices(surf, uSeg0, uSeg1, vSeg0, vSeg1, ddir,
            rseg.curveC0Samples, rseg.curveC1Samples, nRibs, 1.0);

        rseg.ruledMeshVerts = generateRuledMesh(
            rseg.curveC0Samples, rseg.curveC1Samples,
            nCurveSamples, nAcrossSamples, rseg.ruledMeshFaces);

        auto [mx, rms] = computeError(surf, rseg,
            uSeg0, uSeg1, vSeg0, vSeg1, nCurveSamples, nAcrossSamples);
        rseg.maxError = mx; rseg.rmsError = rms;

        result.segments.push_back(rseg);
    }

    return result;
}

} // namespace simple
