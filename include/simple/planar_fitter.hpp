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

// 单格平面拟合（PCA），供网格拟合 fitGridPlanar 使用。
PlaneCellFit fitCellPlane(const SurfaceWrapper& surf,
                          double u0, double u1, double v0, double v1,
                          int nSamplesU, int nSamplesV);

struct PlanarSegment {
    int segmentIndex;
    Vec3 centroid;
    Vec3 normal;
    Vec3Arr meshVerts;
    FaceArr meshFaces;
    double maxError;
    double rmsError;
};

struct PlanarResult {
    std::string name;
    int version;
    ParamDir splitDir;
    std::vector<PlanarSegment> segments;
};

PlanarResult fitPlanarSegments(const SurfaceWrapper& surf,
                                int numSegments,
                                ParamDir splitDir,
                                int nSamplesU, int nSamplesV,
                                int version,
                                const std::string& name);

PlanarResult fitPlanarSegmentsAdaptive(const SurfaceWrapper& surf,
                                       const std::vector<double>& targetDensity,
                                       int initSegments,
                                       ParamDir splitDir,
                                       int nSamplesU, int nSamplesV,
                                       int maxDepth,
                                       const std::string& name);

bool exportPlanarOBJs(const std::string& outDir,
                       const std::vector<PlanarResult>& results);

} // namespace simple
