#pragma once

#include "common.hpp"
#include "surface_wrapper.hpp"

#include <array>

#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>

namespace simple {

struct RuledSegment {
    int segmentIndex;
    ParamDir directrixDir;
    Handle(Geom_BSplineCurve) curveC0;
    Handle(Geom_BSplineCurve) curveC1;
    Vec3Arr curveC0Samples;
    Vec3Arr curveC1Samples;
    Vec3Arr ruledMeshVerts;
    FaceArr ruledMeshFaces;
    double maxError;
    double rmsError;
};

struct RuledResult {
    std::string name;
    int version;
    ParamDir splitDir;
    std::vector<RuledSegment> segments;
    Vec3Arr origMeshVerts;
    FaceArr origMeshFaces;
};

RuledResult fitRuledSegments(const SurfaceWrapper& surf,
                              int numSegments,
                              ParamDir splitDir,
                              const std::vector<ParamDir>& directrixDirs,
                              int nUSamples,
                              int nVSamples,
                              int version,
                              const std::string& name);

Vec3Arr generateRuledMesh(const Vec3Arr& c0Samples,
                           const Vec3Arr& c1Samples,
                           int nU, int nV,
                           FaceArr& faces);

std::pair<double, double> computeError(const SurfaceWrapper& surf,
                                        const RuledSegment& seg,
                                        double uSeg0, double uSeg1,
                                        double vSeg0, double vSeg1,
                                        int nAlong, int nAcross);

void optimizeDirectrices(
    const SurfaceWrapper& surf,
    double uSeg0, double uSeg1,
    double vSeg0, double vSeg1,
    ParamDir directrixDir,
    Vec3Arr& c0Samples,
    Vec3Arr& c1Samples,
    int nRibs, double lambda);

bool exportOBJ(const std::string& path,
               const Vec3Arr& verts, const FaceArr& faces);

bool exportErrorsCSV(const std::string& path,
                     const std::vector<RuledResult>& results);

bool exportMetaJSON(const std::string& path,
                    const std::string& file1, const std::string& file2,
                    const std::vector<RuledResult>& results);

bool exportCurveParamsTXT(const std::string& path,
                          const Handle(Geom_BSplineCurve)& curveC0,
                          const Handle(Geom_BSplineCurve)& curveC1,
                          int nSamples);

} // namespace simple
