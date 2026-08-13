#include "simple/planar_fitter.hpp"

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

PlaneCellFit fitCellPlane(const SurfaceWrapper& surf,
                          double u0, double u1, double v0, double v1,
                          int nSamplesU, int nSamplesV)
{
    PlaneCellFit cell;
    cell.centroid = Vec3(0, 0, 0);
    cell.normal = Vec3(0, 0, 1);
    cell.maxError = 0.0;
    cell.rmsError = 0.0;

    Vec3Arr pts = samplePoints(surf, u0, u1, v0, v1, nSamplesU, nSamplesV);
    if (pts.size() < 3) return cell;

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

    int gr = 8, gc = 8;
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

    cell.centroid = c;
    cell.normal = n;
    cell.meshVerts = mv;
    cell.meshFaces = mf;
    cell.maxError = maxErr;
    cell.rmsError = rmsErr;
    return cell;
}

} // namespace simple
