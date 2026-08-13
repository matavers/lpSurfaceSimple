#include "simple/planar_fitter.hpp"

#include <fstream>
#include <functional>
#include <Eigen/Eigenvalues>

namespace simple {

static Vec3Arr samplePoints(const SurfaceWrapper& surf,
                             double u0, double u1, double v0, double v1,
                             int nU, int nV)
{
    Vec3Arr pts;
    pts.reserve(nU * nV);
    for (int i = 0; i < nU; ++i) {
        for (int j = 0; j < nV; ++j) {
            double u = u0 + (u1 - u0) * i / std::max(1, nU - 1);
            double v = v0 + (v1 - v0) * j / std::max(1, nV - 1);
            pts.push_back(surf.evaluate(u, v));
        }
    }
    return pts;
}

PlanarResult fitPlanarSegments(const SurfaceWrapper& surf,
                                int numSegments,
                                ParamDir splitDir,
                                int nSamplesU, int nSamplesV,
                                int version,
                                const std::string& name)
{
    PlanarResult result;
    result.name = name;
    result.version = version;
    result.splitDir = splitDir;

    auto [uFull0, uFull1] = surf.paramDomainU();
    auto [vFull0, vFull1] = surf.paramDomainV();

    for (int seg = 0; seg < numSegments; ++seg) {
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

        Vec3Arr pts = samplePoints(surf, uSeg0, uSeg1, vSeg0, vSeg1, nSamplesU, nSamplesV);
        if (pts.size() < 3) continue;

        Vec3 c(0, 0, 0);
        for (auto& p : pts) c += p;
        c /= static_cast<double>(pts.size());

        Mat3 cov = Mat3::Zero();
        for (auto& p : pts) {
            Vec3 d = p - c;
            cov += d * d.transpose();
        }
        cov /= static_cast<double>(pts.size());

        Eigen::SelfAdjointEigenSolver<Mat3> eig(cov);
        Vec3 normal = eig.eigenvectors().col(0);
        normal.normalize();

        double maxErr = 0.0, sumSq = 0.0;
        for (auto& p : pts) {
            double dist = std::abs(normal.dot(p - c));
            sumSq += dist * dist;
            if (dist > maxErr) maxErr = dist;
        }
        double rmsErr = std::sqrt(sumSq / pts.size());

        Vec3Arr meshVerts;
        FaceArr meshFaces;
        int gr, gc;
        if (splitDir == ParamDir::V) {
            gr = 24; gc = 3;
        } else {
            gr = 3; gc = 24;
        }
        for (int i = 0; i <= gr; ++i) {
            double t = (double)i / gr;
            double u = uSeg0 + (uSeg1 - uSeg0) * t;
            for (int j = 0; j <= gc; ++j) {
                double s = (double)j / gc;
                double v = vSeg0 + (vSeg1 - vSeg0) * s;
                Vec3 pt = surf.evaluate(u, v);
                double d = normal.dot(pt - c);
                meshVerts.push_back(pt - d * normal);
            }
        }
        for (int i = 0; i < gr; ++i) {
            for (int j = 0; j < gc; ++j) {
                int a = i * (gc + 1) + j;
                meshFaces.push_back({a, a + 1, a + gc + 1});
                meshFaces.push_back({a + 1, a + gc + 2, a + gc + 1});
            }
        }

        PlanarSegment rseg;
        rseg.segmentIndex = seg;
        rseg.centroid = c;
        rseg.normal = normal;
        rseg.meshVerts = meshVerts;
        rseg.meshFaces = meshFaces;
        rseg.maxError = maxErr;
        rseg.rmsError = rmsErr;

        result.segments.push_back(rseg);
    }

    return result;
}

static PlanarSegment fitSinglePlane(const SurfaceWrapper& surf,
                                     double u0, double u1, double v0, double v1,
                                     int nSamplesU, int nSamplesV, int index,
                                     ParamDir splitDir)
{
    Vec3Arr pts = samplePoints(surf, u0, u1, v0, v1, nSamplesU, nSamplesV);
    PlanarSegment seg;
    seg.segmentIndex = index;
    seg.centroid = Vec3(0, 0, 0);
    seg.normal = Vec3(0, 0, 1);
    seg.maxError = 0.0;
    seg.rmsError = 0.0;

    if (pts.size() < 3) return seg;

    Vec3 c(0, 0, 0);
    for (auto& p : pts) c += p;
    c /= static_cast<double>(pts.size());

    Mat3 cov = Mat3::Zero();
    for (auto& p : pts) {
        Vec3 d = p - c;
        cov += d * d.transpose();
    }
    cov /= static_cast<double>(pts.size());

    Eigen::SelfAdjointEigenSolver<Mat3> eig(cov);
    Vec3 n = eig.eigenvectors().col(0);
    n.normalize();

    double maxErr = 0.0, sumSq = 0.0;
    for (auto& p : pts) {
        double dist = std::abs(n.dot(p - c));
        sumSq += dist * dist;
        if (dist > maxErr) maxErr = dist;
    }
    double rmsErr = std::sqrt(sumSq / pts.size());

    int gr, gc;
    if (splitDir == ParamDir::V) {
        gr = 24; gc = 3;
    } else {
        gr = 3; gc = 24;
    }

    Vec3Arr mv;
    FaceArr mf;
    for (int i = 0; i <= gr; ++i) {
        double t = (double)i / gr;
        double u = u0 + (u1 - u0) * t;
        for (int j = 0; j <= gc; ++j) {
            double s = (double)j / gc;
            double v = v0 + (v1 - v0) * s;
            Vec3 pt = surf.evaluate(u, v);
            double d = n.dot(pt - c);
            mv.push_back(pt - d * n);
        }
    }
    for (int i = 0; i < gr; ++i) {
        for (int j = 0; j < gc; ++j) {
            int a = i * (gc + 1) + j;
            mf.push_back({a, a + 1, a + gc + 1});
            mf.push_back({a + 1, a + gc + 2, a + gc + 1});
        }
    }

    seg.centroid = c;
    seg.normal = n;
    seg.meshVerts = mv;
    seg.meshFaces = mf;
    seg.maxError = maxErr;
    seg.rmsError = rmsErr;
    return seg;
}

PlanarResult fitPlanarSegmentsAdaptive(const SurfaceWrapper& surf,
                                       const std::vector<double>& targetDensity,
                                       int initSegments,
                                       ParamDir splitDir,
                                       int nSamplesU, int nSamplesV,
                                       int maxDepth,
                                       const std::string& name)
{
    PlanarResult result;
    result.name = name;
    result.splitDir = splitDir;

    auto [uFull0, uFull1] = surf.paramDomainU();
    auto [vFull0, vFull1] = surf.paramDomainV();
    double totalRange = (splitDir == ParamDir::V) ? (vFull1 - vFull0) : (uFull1 - uFull0);
    double minDomainLen = totalRange * 0.005;
    const int maxTotalSegments = 60;

    int globalIdx = 0;
    std::function<void(double, double, double, double, double, double, int)> subdivide;
    subdivide = [&](double u0, double u1, double v0, double v1,
                     double domLen, double target, int depth)
    {
        if (domLen < minDomainLen || globalIdx >= maxTotalSegments) {
            PlanarSegment seg = fitSinglePlane(surf, u0, u1, v0, v1,
                                                nSamplesU, nSamplesV, globalIdx, splitDir);
            seg.segmentIndex = globalIdx++;
            result.segments.push_back(seg);
            return;
        }

        PlanarSegment seg = fitSinglePlane(surf, u0, u1, v0, v1,
                                            nSamplesU, nSamplesV, globalIdx, splitDir);

        if (seg.maxError <= target || depth <= 0) {
            seg.segmentIndex = globalIdx++;
            result.segments.push_back(seg);
            return;
        }

        if (splitDir == ParamDir::V) {
            double vm = (v0 + v1) * 0.5;
            double half = domLen * 0.5;

            auto segA = fitSinglePlane(surf, u0, u1, v0, vm, nSamplesU, nSamplesV, 0, splitDir);
            auto segB = fitSinglePlane(surf, u0, u1, vm, v1, nSamplesU, nSamplesV, 0, splitDir);
            if (std::max(segA.maxError, segB.maxError) > seg.maxError * 0.95) {
                seg.segmentIndex = globalIdx++;
                result.segments.push_back(seg);
                return;
            }
            subdivide(u0, u1, v0, vm, half, target, depth - 1);
            subdivide(u0, u1, vm, v1, half, target, depth - 1);
        } else {
            double um = (u0 + u1) * 0.5;
            double half = domLen * 0.5;

            auto segA = fitSinglePlane(surf, u0, um, v0, v1, nSamplesU, nSamplesV, 0, splitDir);
            auto segB = fitSinglePlane(surf, um, u1, v0, v1, nSamplesU, nSamplesV, 0, splitDir);
            if (std::max(segA.maxError, segB.maxError) > seg.maxError * 0.95) {
                seg.segmentIndex = globalIdx++;
                result.segments.push_back(seg);
                return;
            }
            subdivide(u0, um, v0, v1, half, target, depth - 1);
            subdivide(um, u1, v0, v1, half, target, depth - 1);
        }
    };

    for (int seg = 0; seg < initSegments; ++seg) {
        double u0, u1, v0, v1, domLen;
        if (splitDir == ParamDir::V) {
            u0 = uFull0; u1 = uFull1;
            v0 = vFull0 + (vFull1 - vFull0) * seg / initSegments;
            v1 = vFull0 + (vFull1 - vFull0) * (seg + 1) / initSegments;
            domLen = (vFull1 - vFull0) / initSegments;
        } else {
            u0 = uFull0 + (uFull1 - uFull0) * seg / initSegments;
            u1 = uFull0 + (uFull1 - uFull0) * (seg + 1) / initSegments;
            v0 = vFull0; v1 = vFull1;
            domLen = (uFull1 - uFull0) / initSegments;
        }
        double target = (seg < (int)targetDensity.size())
                            ? targetDensity[seg]
                            : targetDensity.back();
        subdivide(u0, u1, v0, v1, domLen, target, maxDepth);
    }

    return result;
}

bool exportOBJ(const std::string& path,
               const Vec3Arr& verts, const FaceArr& faces);

bool exportPlanarOBJs(const std::string& outDir,
                       const std::vector<PlanarResult>& results)
{
    for (size_t bi = 0; bi < results.size(); ++bi) {
        std::string prefix = (bi == 0) ? "blade1" : "blade2";
        for (const auto& seg : results[bi].segments) {
            std::string path = outDir + "/" + prefix
                             + "_plane" + std::to_string(seg.segmentIndex) + ".obj";
            exportOBJ(path, seg.meshVerts, seg.meshFaces);
        }
    }
    return true;
}

} // namespace simple
