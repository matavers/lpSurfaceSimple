#pragma once

#include "common.hpp"

#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>

namespace simple {

struct StepLoadResult {
    std::string filename;
    std::vector<TopoDS_Face> faces;
    bool loaded;
    std::string errorMsg;
};

StepLoadResult loadStepFile(const std::string& filepath);

TopoDS_Face findLargestFace(const TopoDS_Shape& shape);

void generateFaceMesh(const TopoDS_Face& face, double deflection,
                      Vec3Arr& vertices, std::vector<std::array<int, 3>>& faces);

} // namespace simple
