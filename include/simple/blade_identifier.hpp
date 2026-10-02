#pragma once

#include "common.hpp"

#include <TopoDS_Face.hxx>

#include <string>
#include <vector>

namespace simple {

struct BladeIdentifyResult {
    int pressureFaceIndex = -1;
    int suctionFaceIndex = -1;
    bool success = false;
    std::string pressureLabel;
    std::string suctionLabel;
    std::string message;
};

BladeIdentifyResult identifyBladeSurfaces(
    const std::vector<TopoDS_Face>& faces,
    int normalSampleGrid = 9,
    double minParamArea = 0.01,
    double minDiag3D = 20.0,
    int curvatureSampleGrid = 9);

} // namespace simple
