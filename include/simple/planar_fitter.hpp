#pragma once

#include "common.hpp"
#include "surface_wrapper.hpp"

namespace simple {

// 单个格子（参数域矩形）的平面拟合结果（PCA）。
struct PlaneCellFit {
    Vec3 centroid;
    Vec3 normal;
    Vec3Arr meshVerts;
    FaceArr meshFaces;
    double maxError;
    double rmsError;
};

PlaneCellFit fitCellPlane(const SurfaceWrapper& surf,
                          double u0, double u1, double v0, double v1,
                          int nSamplesU, int nSamplesV);

} // namespace simple
