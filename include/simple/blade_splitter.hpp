#pragma once

#include "common.hpp"

#include <TopoDS_Face.hxx>

#include <string>

namespace simple {

struct BladeSplitResult {
    double uPressStart = 0, uPressEnd = 0;
    double uSuctStart = 0, uSuctEnd = 0;
    double uLE = 0, uTE = 0;
    double leCurv = 0, teCurv = 0;
    bool success = false;
    std::string message;
};

BladeSplitResult splitBladeFaceBySection(
    const TopoDS_Face& face,
    int numSections = 7,
    int samplesPerSection = 200,
    int smoothingWindow = 7);

} // namespace simple
