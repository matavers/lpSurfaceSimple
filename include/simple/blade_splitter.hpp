#pragma once

#include "common.hpp"

#include <TopoDS_Face.hxx>

#include <string>
#include <vector>

namespace simple {

struct SplitRegion {
    double uStart = 0, uEnd = 0;
    double avgCurv = 0;
    std::string label;
};

struct BladeSplitResult {
    std::vector<SplitRegion> regions;
    bool success = false;
    std::string message;
};

BladeSplitResult splitBladeFaceBySection(
    const TopoDS_Face& face,
    int numSections = 7,
    int samplesPerSection = 200,
    int smoothingWindow = 7,
    double peakThresholdRatio = 0.15,
    double clusterGapRatio = 0.12);

} // namespace simple
