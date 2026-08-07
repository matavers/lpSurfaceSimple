#include "simple/planar_fitter.hpp"

#include <fstream>
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

        Vec3 corners[4] = {
            surf.evaluate(uSeg0, vSeg0),
            surf.evaluate(uSeg1, vSeg0),
            surf.evaluate(uSeg1, vSeg1),
            surf.evaluate(uSeg0, vSeg1)
        };
        Vec3Arr meshVerts;
        for (int i = 0; i < 4; ++i) {
            double d = normal.dot(corners[i] - c);
            meshVerts.push_back(corners[i] - d * normal);
        }

        PlanarSegment rseg;
        rseg.segmentIndex = seg;
        rseg.centroid = c;
        rseg.normal = normal;
        rseg.meshVerts = meshVerts;
        rseg.meshFaces = {{0, 1, 2}, {0, 2, 3}};
        rseg.maxError = maxErr;
        rseg.rmsError = rmsErr;

        result.segments.push_back(rseg);
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
