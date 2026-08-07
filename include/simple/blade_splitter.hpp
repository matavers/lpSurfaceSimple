#pragma once

#include "common.hpp"

#include <TopoDS_Face.hxx>

#include <string>

namespace simple {

struct BladeSplitResult {
    double uPeak1 = 0, uPeak2 = 0;
    double curvPeak1 = 0, curvPeak2 = 0;

    double uEdgeStart = 0, uEdgeEnd = 0;
    double uPressStart = 0, uPressEnd = 0;
    double uSuctStart = 0, uSuctEnd = 0;

    double curvPress = 0, curvSuct = 0;

    bool success = false;
    std::string message;
};

BladeSplitResult splitBladeFaceBySection(
    const TopoDS_Face& face,
    int numSections = 7,
    int samplesPerSection = 200,
    int smoothingWindow = 7);

} // namespace simple
