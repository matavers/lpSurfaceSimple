#pragma once

#include "common.hpp"
#include "surface_wrapper.hpp"

namespace simple {

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
