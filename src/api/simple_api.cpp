#include "ruledSurfaceFitting.h"

#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <filesystem>

#include "simple/common.hpp"
#include "simple/step_reader.hpp"
#include "simple/surface_wrapper.hpp"
#include "simple/ruled_fitter.hpp"
#include "simple/planar_fitter.hpp"
#include "simple/grid_fitter.hpp"
#include "simple/blade_identifier.hpp"
#include "simple/blade_splitter.hpp"

#include <BRepAdaptor_Surface.hxx>
#include <BRep_Tool.hxx>
#include <GeomConvert.hxx>
#include <Geom_Surface.hxx>
#include <TopExp_Explorer.hxx>

using namespace simple;

namespace {

TopoDS_Face pickFace(const std::vector<TopoDS_Face>& faces, int forceIdx) {
    if (faces.empty()) return TopoDS_Face();
    if (forceIdx >= 0 && forceIdx < (int)faces.size()) return faces[forceIdx];
    if (faces.size() == 1) return faces[0];
    int bestIdx = 0; double bestDiag = 0.0;
    for (size_t i = 0; i < faces.size(); ++i) {
        BRepAdaptor_Surface a(faces[i], true);
        Vec3 bmin(1e30,1e30,1e30), bmax(-1e30,-1e30,-1e30);
        for (int ui = 0; ui <= 4; ++ui)
            for (int vi = 0; vi <= 4; ++vi) {
                double u = a.FirstUParameter() + (a.LastUParameter()-a.FirstUParameter())*ui/4.0;
                double v = a.FirstVParameter() + (a.LastVParameter()-a.FirstVParameter())*vi/4.0;
                gp_Pnt p = a.Value(u, v);
                bmin=Vec3(std::min(bmin.x(),p.X()),std::min(bmin.y(),p.Y()),std::min(bmin.z(),p.Z()));
                bmax=Vec3(std::max(bmax.x(),p.X()),std::max(bmax.y(),p.Y()),std::max(bmax.z(),p.Z()));
            }
        double d = (bmax-bmin).norm();
        if (d > bestDiag) { bestDiag = d; bestIdx = (int)i; }
    }
    return faces[bestIdx];
}

Handle(Geom_BSplineSurface) faceToSurf(const TopoDS_Face& f) {
    BRepAdaptor_Surface a(f, true);
    if (a.GetType() == GeomAbs_BSplineSurface) return a.BSpline();
    return GeomConvert::SurfaceToBSplineSurface(BRep_Tool::Surface(f));
}

GridConfig makeGridConfig(int nUSamples, int nVSamples, int nRibs, double lambda,
                          int nSplitU, int nSplitV, double tolerance, int maxDepth) {
    GridConfig c;
    c.nUSamples = nUSamples;
    c.nVSamples = nVSamples;
    c.nRibs = nRibs;
    c.lambda = lambda;
    c.nSplitU = nSplitU;
    c.nSplitV = nSplitV;
    c.tolerance = tolerance;
    c.maxDepth = maxDepth;
    return c;
}

void fillResult(RuledFittingResult* res, const std::vector<GridResult>& grs,
                const std::string& mode) {
    res->errorCode = RULED_OK;
    res->numSurfaces = (int)grs.size();
    std::ostringstream summary;
    summary << "{\"ok\":true,\"mode\":\"" << mode << "\",\"surfaces\":[";
    for (size_t si = 0; si < grs.size(); ++si) {
        const auto& g = grs[si];
        if (si) summary << ",";
        std::strncpy(res->surfaces[si].name, g.name.c_str(), 63);
        res->surfaces[si].name[63] = 0;
        res->surfaces[si].numCells = (int)std::min<size_t>(g.cells.size(), 512);
        double worst = 0.0;
        for (size_t k = 0; k < g.cells.size() && k < 512; ++k) {
            const auto& c = g.cells[k];
            res->surfaces[si].cells[k].index = c.row * g.nCols + c.col;
            res->surfaces[si].cells[k].row = c.row;
            res->surfaces[si].cells[k].col = c.col;
            res->surfaces[si].cells[k].fitDir = (c.fitDir == ParamDir::U) ? RULED_DIR_U : RULED_DIR_V;
            res->surfaces[si].cells[k].maxError = c.maxError;
            res->surfaces[si].cells[k].rmsError = c.rmsError;
            worst = std::max(worst, c.maxError);
        }
        summary << "{\"name\":\"" << g.name << "\",\"numCells\":" << g.cells.size()
                << ",\"maxErr\":" << worst << "}";
    }
    summary << "]}";
    std::string s = summary.str();
    std::strncpy(res->metaJson, s.c_str(), sizeof(res->metaJson) - 1);
}

void writeMetaFile(const std::string& path, const std::vector<GridResult>& grs,
                   const std::string& mode, const std::string& f1, const std::string& f2) {
    std::ofstream o(path);
    if (o) o << buildGridMetaJson(grs, mode, f1, f2) << "\n";
}

void exportMeshAndGrid(const std::string& od, const TopoDS_Face& face,
                       const BRepAdaptor_Surface& ad, const SurfaceWrapper& sw,
                       const std::string& prefix,
                       const GridResult& gr, bool planar) {
    Vec3Arr mv; std::vector<std::array<int,3>> mf;
    double defl = std::max(0.2, (ad.LastUParameter() - ad.FirstUParameter()
                               + ad.LastVParameter() - ad.FirstVParameter()) * 0.01);
    generateFaceMesh(face, defl, mv, mf);
    exportOBJ(od + "/" + prefix + "_mesh.obj", mv, mf);
    exportGridOBJs(od, prefix, gr, planar);
    exportGridLinesVTK(od + "/" + prefix + "_grid.vtk", sw, gr);
}

} // anon

extern "C" {

RULED_API RuledFittingResult* ruled_surface_fitting(const RuledConfig* cfg) {
    RuledFittingResult* res = new RuledFittingResult(); std::memset(res, 0, sizeof(RuledFittingResult));
    if (!cfg || !cfg->stepFile1 || !cfg->stepFile2) {
        res->errorCode = RULED_ERR_INVALID_PARAMS; return res;
    }
    RuledConfig c = *cfg;
    if (c.nUSamples <= 0) c.nUSamples = 50;
    if (c.nVSamples <= 0) c.nVSamples = 10;
    if (c.nRibs <= 0) c.nRibs = 20;
    if (c.lambda <= 0.0) c.lambda = 1.0;
    if (c.nSplitU <= 0) c.nSplitU = 2;
    if (c.nSplitV <= 0) c.nSplitV = 2;
    if (c.tolerance <= 0.0) c.tolerance = 0.1;
    if (c.maxDepth <= 0) c.maxDepth = 20;

    std::filesystem::create_directories(std::filesystem::path(c.outputDir ? c.outputDir : "."));

    StepLoadResult r1 = loadStepFile(c.stepFile1), r2 = loadStepFile(c.stepFile2);
    if (!r1.loaded || !r2.loaded) { res->errorCode = RULED_ERR_STEP_READ_FAILED; return res; }

    TopoDS_Face f1 = pickFace(r1.faces, c.faceIdx[0]), f2 = pickFace(r2.faces, c.faceIdx[1]);
    if (f1.IsNull() || f2.IsNull()) { res->errorCode = RULED_ERR_NO_VALID_FACE; return res; }

    Handle(Geom_BSplineSurface) sf1 = faceToSurf(f1), sf2 = faceToSurf(f2);
    BRepAdaptor_Surface a1(f1, true), a2(f2, true);

    SurfaceWrapper sw1(sf1, a1.FirstUParameter(), a1.LastUParameter(), a1.FirstVParameter(), a1.LastVParameter());
    SurfaceWrapper sw2(sf2, a2.FirstUParameter(), a2.LastUParameter(), a2.FirstVParameter(), a2.LastVParameter());

    GridConfig g = makeGridConfig(c.nUSamples, c.nVSamples, c.nRibs, c.lambda,
                                  c.nSplitU, c.nSplitV, c.tolerance, c.maxDepth);
    GridResult gr1 = fitGridRuled(sw1, g, "Blade-1");
    GridResult gr2 = fitGridRuled(sw2, g, "Blade-2");

    std::string od = c.outputDir ? c.outputDir : ".";
    exportMeshAndGrid(od, f1, a1, sw1, "blade1", gr1, false);
    exportMeshAndGrid(od, f2, a2, sw2, "blade2", gr2, false);
    writeMetaFile(od + "/meta.json", {gr1, gr2}, "ruled", c.stepFile1, c.stepFile2);

    fillResult(res, {gr1, gr2}, "ruled");
    return res;
}

RULED_API RuledFittingResult* plane_surface_fitting(const PlanarConfig* cfg) {
    RuledFittingResult* res = new RuledFittingResult(); std::memset(res, 0, sizeof(RuledFittingResult));
    if (!cfg || !cfg->stepFile1 || !cfg->stepFile2) {
        res->errorCode = RULED_ERR_INVALID_PARAMS; return res;
    }
    PlanarConfig c = *cfg;
    if (c.nUSamples <= 0) c.nUSamples = 50;
    if (c.nVSamples <= 0) c.nVSamples = 10;
    if (c.nSplitU <= 0) c.nSplitU = 2;
    if (c.nSplitV <= 0) c.nSplitV = 2;
    if (c.tolerance <= 0.0) c.tolerance = 0.1;
    if (c.maxDepth <= 0) c.maxDepth = 20;

    std::filesystem::create_directories(std::filesystem::path(c.outputDir ? c.outputDir : "."));

    StepLoadResult r1 = loadStepFile(c.stepFile1), r2 = loadStepFile(c.stepFile2);
    if (!r1.loaded || !r2.loaded) { res->errorCode = RULED_ERR_STEP_READ_FAILED; return res; }

    TopoDS_Face f1 = pickFace(r1.faces, c.faceIdx[0]), f2 = pickFace(r2.faces, c.faceIdx[1]);
    if (f1.IsNull() || f2.IsNull()) { res->errorCode = RULED_ERR_NO_VALID_FACE; return res; }

    Handle(Geom_BSplineSurface) sf1 = faceToSurf(f1), sf2 = faceToSurf(f2);
    BRepAdaptor_Surface a1(f1, true), a2(f2, true);

    SurfaceWrapper sw1(sf1, a1.FirstUParameter(), a1.LastUParameter(), a1.FirstVParameter(), a1.LastVParameter());
    SurfaceWrapper sw2(sf2, a2.FirstUParameter(), a2.LastUParameter(), a2.FirstVParameter(), a2.LastVParameter());

    GridConfig g = makeGridConfig(c.nUSamples, c.nVSamples, 20, 1.0,
                                  c.nSplitU, c.nSplitV, c.tolerance, c.maxDepth);
    GridResult gr1 = fitGridPlanar(sw1, g, "Blade-1");
    GridResult gr2 = fitGridPlanar(sw2, g, "Blade-2");

    std::string od = c.outputDir ? c.outputDir : ".";
    exportMeshAndGrid(od, f1, a1, sw1, "blade1", gr1, true);
    exportMeshAndGrid(od, f2, a2, sw2, "blade2", gr2, true);
    writeMetaFile(od + "/meta.json", {gr1, gr2}, "planar", c.stepFile1, c.stepFile2);

    fillResult(res, {gr1, gr2}, "planar");
    return res;
}

namespace {

bool prepareBladeSurface(const char* filepath, bool wantPressure,
                          TopoDS_Face& face, double& vStart, double& vEnd,
                          int& faceIdx, std::string& errMsg)
{
    StepLoadResult r = loadStepFile(filepath);
    if (!r.loaded) { errMsg = r.errorMsg; return false; }

    BladeIdentifyResult ident = identifyBladeSurfaces(r.faces);
    if (!ident.success) { errMsg = ident.message; return false; }

    faceIdx = ident.pressureFaceIndex;
    if (faceIdx < 0 || faceIdx >= (int)r.faces.size()) {
        errMsg = "Invalid face index from auto-identify"; return false;
    }
    face = r.faces[faceIdx];

    BladeSplitResult split = splitBladeFaceBySection(face);
    if (!split.success) { errMsg = split.message; return false; }

    vStart = -1; vEnd = -1;
    for (auto& reg : split.regions) {
        if ((wantPressure && reg.label == "pressure") ||
            (!wantPressure && reg.label == "suction")) {
            vStart = reg.uStart; vEnd = reg.uEnd; break;
        }
    }
    if (vStart < 0) {
        errMsg = wantPressure ? "Pressure region not found in split"
                              : "Suction region not found in split";
        return false;
    }
    return true;
}

RuledFittingResult* runGridOnSide(const char* filepath, bool wantPressure,
                                  bool planar, int nUSamples, int nVSamples,
                                  int nSplitU, int nSplitV,
                                  double tolerance, int maxDepth)
{
    RuledFittingResult* res = new RuledFittingResult();
    std::memset(res, 0, sizeof(RuledFittingResult));

    if (!filepath) { res->errorCode = RULED_ERR_INVALID_PARAMS; return res; }

    TopoDS_Face face;
    double v0, v1;
    int fi;
    std::string err;
    if (!prepareBladeSurface(filepath, wantPressure, face, v0, v1, fi, err)) {
        res->errorCode = RULED_ERR_NO_VALID_FACE;
        std::strncpy(res->errorMsg, err.c_str(), 255);
        return res;
    }

    Handle(Geom_BSplineSurface) sf = faceToSurf(face);
    BRepAdaptor_Surface ad(face, true);
    SurfaceWrapper sw(sf, ad.FirstUParameter(), ad.LastUParameter(), v0, v1, false);

    const char* sideName = wantPressure ? "Pressure" : "Suction";
    GridConfig g = makeGridConfig(nUSamples, nVSamples, 20, 1.0,
                                  nSplitU, nSplitV, tolerance, maxDepth);
    GridResult gr = planar ? fitGridPlanar(sw, g, sideName)
                           : fitGridRuled(sw, g, sideName);

    std::string od = ".";
    // outputDir resolution handled by caller via config; here default to cwd-relative
    std::string pfx = wantPressure ? "pressure" : "suction";
    exportMeshAndGrid(od, face, ad, sw, pfx, gr, planar);
    writeMetaFile(od + "/meta.json", {gr}, planar ? "planar" : "ruled", filepath, "");

    fillResult(res, {gr}, planar ? "planar" : "ruled");
    return res;
}

} // anon

RULED_API RuledFittingResult* pressure_ruled_fitting(const RuledConfig* cfg) {
    if (!cfg || !cfg->outputDir) { RuledFittingResult* r = new RuledFittingResult(); std::memset(r,0,sizeof(*r)); r->errorCode = RULED_ERR_INVALID_PARAMS; return r; }
    RuledConfig c = *cfg;
    if (c.nUSamples <= 0) c.nUSamples = 50;
    if (c.nVSamples <= 0) c.nVSamples = 10;
    if (c.nSplitU <= 0) c.nSplitU = 2;
    if (c.nSplitV <= 0) c.nSplitV = 2;
    if (c.tolerance <= 0.0) c.tolerance = 0.1;
    if (c.maxDepth <= 0) c.maxDepth = 20;
    std::filesystem::create_directories(std::filesystem::path(c.outputDir));
    return runGridOnSide(c.stepFile1, true, false, c.nUSamples, c.nVSamples,
                         c.nSplitU, c.nSplitV, c.tolerance, c.maxDepth);
}
RULED_API RuledFittingResult* pressure_plane_fitting(const PlanarConfig* cfg) {
    if (!cfg || !cfg->outputDir) { RuledFittingResult* r = new RuledFittingResult(); std::memset(r,0,sizeof(*r)); r->errorCode = RULED_ERR_INVALID_PARAMS; return r; }
    PlanarConfig c = *cfg;
    if (c.nUSamples <= 0) c.nUSamples = 50;
    if (c.nVSamples <= 0) c.nVSamples = 10;
    if (c.nSplitU <= 0) c.nSplitU = 2;
    if (c.nSplitV <= 0) c.nSplitV = 2;
    if (c.tolerance <= 0.0) c.tolerance = 0.1;
    if (c.maxDepth <= 0) c.maxDepth = 20;
    std::filesystem::create_directories(std::filesystem::path(c.outputDir));
    return runGridOnSide(c.stepFile1, true, true, c.nUSamples, c.nVSamples,
                         c.nSplitU, c.nSplitV, c.tolerance, c.maxDepth);
}
RULED_API RuledFittingResult* suction_ruled_fitting(const RuledConfig* cfg) {
    if (!cfg || !cfg->outputDir) { RuledFittingResult* r = new RuledFittingResult(); std::memset(r,0,sizeof(*r)); r->errorCode = RULED_ERR_INVALID_PARAMS; return r; }
    RuledConfig c = *cfg;
    if (c.nUSamples <= 0) c.nUSamples = 50;
    if (c.nVSamples <= 0) c.nVSamples = 10;
    if (c.nSplitU <= 0) c.nSplitU = 2;
    if (c.nSplitV <= 0) c.nSplitV = 2;
    if (c.tolerance <= 0.0) c.tolerance = 0.1;
    if (c.maxDepth <= 0) c.maxDepth = 20;
    std::filesystem::create_directories(std::filesystem::path(c.outputDir));
    return runGridOnSide(c.stepFile1, false, false, c.nUSamples, c.nVSamples,
                         c.nSplitU, c.nSplitV, c.tolerance, c.maxDepth);
}
RULED_API RuledFittingResult* suction_plane_fitting(const PlanarConfig* cfg) {
    if (!cfg || !cfg->outputDir) { RuledFittingResult* r = new RuledFittingResult(); std::memset(r,0,sizeof(*r)); r->errorCode = RULED_ERR_INVALID_PARAMS; return r; }
    PlanarConfig c = *cfg;
    if (c.nUSamples <= 0) c.nUSamples = 50;
    if (c.nVSamples <= 0) c.nVSamples = 10;
    if (c.nSplitU <= 0) c.nSplitU = 2;
    if (c.nSplitV <= 0) c.nSplitV = 2;
    if (c.tolerance <= 0.0) c.tolerance = 0.1;
    if (c.maxDepth <= 0) c.maxDepth = 20;
    std::filesystem::create_directories(std::filesystem::path(c.outputDir));
    return runGridOnSide(c.stepFile1, false, true, c.nUSamples, c.nVSamples,
                         c.nSplitU, c.nSplitV, c.tolerance, c.maxDepth);
}

RULED_API void free_result(RuledFittingResult* result) { delete result; }

} // extern "C"
