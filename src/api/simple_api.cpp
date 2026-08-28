#include "ruledSurfaceFitting.h"

#include <cstring>
#include <cctype>
#include <string>
#include <sstream>
#include <filesystem>
#include <algorithm>

#include "simple/common.hpp"
#include "simple/step_reader.hpp"
#include "simple/surface_wrapper.hpp"
#include "simple/ruled_fitter.hpp"
#include "simple/grid_fitter.hpp"
#include "simple/blade_identifier.hpp"
#include "simple/blade_splitter.hpp"

#include <BRepAdaptor_Surface.hxx>
#include <BRep_Tool.hxx>
#include <GeomConvert.hxx>
#include <Geom_Surface.hxx>
#include <TopoDS_Face.hxx>

using namespace simple;
namespace fs = std::filesystem;

namespace {

struct Defaults {
    static constexpr int nUSamples = 50;
    static constexpr int nVSamples = 10;
    static constexpr int nRibs = 20;
    static constexpr double lambda = 1.0;
    static constexpr double tolerance = 0.1;
};

Handle(Geom_BSplineSurface) faceToSurf(const TopoDS_Face& f) {
    BRepAdaptor_Surface a(f, true);
    if (a.GetType() == GeomAbs_BSplineSurface) return a.BSpline();
    return GeomConvert::SurfaceToBSplineSurface(BRep_Tool::Surface(f));
}

// 取压力/吸力侧在弦向分割中的参数区间（沿用 splitBladeFaceBySection 的 uStart/uEnd 语义，即 V 区间）
bool sideVRange(const TopoDS_Face& face, bool wantPressure,
                double& v0, double& v1, std::string& err) {
    BladeSplitResult split = splitBladeFaceBySection(face);
    if (!split.success) { err = split.message; return false; }
    v0 = -1; v1 = -1;
    for (const auto& reg : split.regions) {
        if ((wantPressure && reg.label == "pressure") ||
            (!wantPressure && reg.label == "suction")) {
            v0 = reg.uStart; v1 = reg.uEnd; break;
        }
    }
    if (v0 < 0) {
        err = wantPressure ? "pressure region not found" : "suction region not found";
        return false;
    }
    return true;
}

struct SideOutcome {
    bool ok = false;
    int numPieces = 0;
    double maxError = 0.0;
    std::vector<double> pieceMaxErr;
    std::vector<double> pieceRmsErr;
};

// 拟合单个侧面：先 3 等分；若最大误差 < tolerance 直接输出，否则井字形网格细分。
SideOutcome fitAndExportSide(const SurfaceWrapper& sw, const std::string& prefix,
                             const std::string& outDir, double tolerance) {
    SideOutcome out;
    const int nU = Defaults::nUSamples, nV = Defaults::nVSamples, nRibs = Defaults::nRibs;
    const double lam = Defaults::lambda;

    std::vector<ParamDir> dd = { ParamDir::V, ParamDir::V, ParamDir::V };
    RuledResult rr = fitRuledSegments(sw, 3, ParamDir::V, dd, nU, nV, 0, prefix);

    double maxErr = 0.0;
    for (const auto& seg : rr.segments) maxErr = std::max(maxErr, seg.maxError);

    if (maxErr < tolerance) {
        for (size_t i = 0; i < rr.segments.size(); ++i) {
            const auto& seg = rr.segments[i];
            exportOBJ(outDir + "/" + prefix + "_seg" + std::to_string(i) + ".obj",
                      seg.ruledMeshVerts, seg.ruledMeshFaces);
            exportDirectrixTXT(outDir + "/" + prefix + "_seg" + std::to_string(i) + "_params.txt",
                               seg.curveC0Samples, seg.curveC1Samples);
            out.pieceMaxErr.push_back(seg.maxError);
            out.pieceRmsErr.push_back(seg.rmsError);
        }
        out.numPieces = (int)rr.segments.size();
        out.maxError = maxErr;
        out.ok = true;
        return out;
    }

    GridConfig gcfg;
    gcfg.nSplitU = 2;
    gcfg.nSplitV = 2;
    gcfg.tolerance = tolerance;
    gcfg.nUSamples = nU;
    gcfg.nVSamples = nV;
    gcfg.nRibs = nRibs;
    gcfg.lambda = lam;

    GridResult gr = fitGridRuled(sw, gcfg, prefix);
    for (const auto& cell : gr.cells) {
        int idx = cell.row * gr.nCols + cell.col;
        exportOBJ(outDir + "/" + prefix + "_seg" + std::to_string(idx) + ".obj",
                  cell.ruled.ruledMeshVerts, cell.ruled.ruledMeshFaces);
        exportDirectrixTXT(outDir + "/" + prefix + "_seg" + std::to_string(idx) + "_params.txt",
                           cell.ruled.curveC0Samples, cell.ruled.curveC1Samples);
        out.pieceMaxErr.push_back(cell.maxError);
        out.pieceRmsErr.push_back(cell.rmsError);
    }
    out.numPieces = (int)gr.cells.size();
    out.maxError = gr.maxError;
    out.ok = true;
    return out;
}

std::string jsonSafeStr(std::string s) {
    for (char& c : s) { if (c == '\\') c = '/'; if (c == '"') c = '\''; }
    return s;
}

std::string buildSummaryJson(const std::string& mode, const std::string& file,
                             const std::vector<std::pair<std::string, SideOutcome>>& sides,
                             double tolerance) {
    std::ostringstream o;
    o << "{\"mode\":\"" << mode << "\",\"tolerance\":" << tolerance
      << ",\"file\":\"" << jsonSafeStr(file) << "\",\"surfaces\":[";
    for (size_t i = 0; i < sides.size(); ++i) {
        if (i) o << ",";
        o << "{\"name\":\"" << sides[i].first << "\",\"numSegments\":" << sides[i].second.numPieces
          << ",\"maxError\":" << sides[i].second.maxError << "}";
    }
    o << "]}";
    return o.str();
}

bool hasSupportedExt(const fs::path& p) {
    std::string ext = p.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(),
                   [](unsigned char c) { return (char)std::tolower(c); });
    return ext == ".step" || ext == ".stp" || ext == ".igs" || ext == ".iges";
}

} // anon

extern "C" {

RULED_API RuledFittingResult* ruled_fitting(const RuledFitConfig* cfg) {
    RuledFittingResult* res = new RuledFittingResult();
    std::memset(res, 0, sizeof(*res));

    if (!cfg || !cfg->inputPath || !cfg->inputPath[0] || !cfg->outputDir || !cfg->outputDir[0]) {
        res->errorCode = RULED_ERR_INVALID_PARAMS;
        return res;
    }

    double tol = cfg->tolerance;
    if (tol <= 0.0) tol = Defaults::tolerance;

    std::string inPath = cfg->inputPath;
    std::string outDir = cfg->outputDir;
    std::error_code ec;
    fs::create_directories(outDir, ec);

    StepLoadResult r = loadStepFile(inPath);
    if (!r.loaded) {
        res->errorCode = RULED_ERR_STEP_READ_FAILED;
        std::strncpy(res->errorMsg, r.errorMsg.c_str(), sizeof(res->errorMsg) - 1);
        return res;
    }

    BladeIdentifyResult ident = identifyBladeSurfaces(r.faces);
    if (!ident.success) {
        res->errorCode = RULED_ERR_NO_VALID_FACE;
        std::strncpy(res->errorMsg, ident.message.c_str(), sizeof(res->errorMsg) - 1);
        return res;
    }

    struct SideSpec { int faceIdx; bool pressure; const char* prefix; const char* name; };
    SideSpec sides[2] = {
        { ident.pressureFaceIndex, true,  "pressure", "Pressure" },
        { ident.suctionFaceIndex,  false, "suction",  "Suction"  },
    };

    std::vector<std::pair<std::string, SideOutcome>> outcomes;
    for (int si = 0; si < 2; ++si) {
        int fi = sides[si].faceIdx;
        if (fi < 0 || fi >= (int)r.faces.size()) continue;

        TopoDS_Face face = r.faces[fi];
        double v0, v1;
        std::string err;
        if (!sideVRange(face, sides[si].pressure, v0, v1, err)) continue;

        Handle(Geom_BSplineSurface) sf = faceToSurf(face);
        BRepAdaptor_Surface ad(face, true);
        SurfaceWrapper sw(sf, ad.FirstUParameter(), ad.LastUParameter(), v0, v1, false);

        SideOutcome oc = fitAndExportSide(sw, sides[si].prefix, outDir, tol);
        if (!oc.ok) continue;

        RuledSurfaceResult& sr = res->surfaces[res->numSurfaces];
        std::strncpy(sr.name, sides[si].name, sizeof(sr.name) - 1);
        sr.maxError = oc.maxError;
        sr.numSegments = oc.numPieces;
        int cap = std::min(oc.numPieces, 64);
        for (int k = 0; k < cap; ++k) {
            sr.segments[k].segmentIndex = k;
            sr.segments[k].maxError = oc.pieceMaxErr[k];
            sr.segments[k].rmsError = oc.pieceRmsErr[k];
        }
        res->numSurfaces++;
        outcomes.push_back({ sides[si].name, oc });
    }

    if (res->numSurfaces == 0) {
        res->errorCode = RULED_ERR_NO_VALID_FACE;
        std::strncpy(res->errorMsg, "no pressure/suction surface fitted", sizeof(res->errorMsg) - 1);
        return res;
    }

    std::string meta = buildSummaryJson("ruled", inPath, outcomes, tol);
    std::strncpy(res->metaJson, meta.c_str(), sizeof(res->metaJson) - 1);

    res->errorCode = RULED_OK;
    return res;
}

RULED_API RuledFittingResult* ruled_fitting_simple(const char* inputDir, const char* outputDir) {
    RuledFittingResult* res = new RuledFittingResult();
    std::memset(res, 0, sizeof(*res));

    if (!inputDir || !inputDir[0] || !outputDir || !outputDir[0]) {
        res->errorCode = RULED_ERR_INVALID_PARAMS;
        return res;
    }

    std::string inDir = inputDir;
    std::string outDir = outputDir;
    std::error_code ec;
    fs::create_directories(outDir, ec);

    const int nU = Defaults::nUSamples, nV = Defaults::nVSamples;
    std::vector<ParamDir> dd = { ParamDir::V, ParamDir::V, ParamDir::V };

    int processed = 0;
    for (const auto& entry : fs::directory_iterator(inDir)) {
        if (!entry.is_regular_file() || !hasSupportedExt(entry.path())) continue;

        std::string file = entry.path().string();
        std::string stem = entry.path().stem().string();

        StepLoadResult r = loadStepFile(file);
        if (!r.loaded) continue;

        BladeIdentifyResult ident = identifyBladeSurfaces(r.faces);
        if (!ident.success) continue;

        struct SideSpec { int faceIdx; bool pressure; const char* prefix; };
        SideSpec sides[2] = {
            { ident.pressureFaceIndex, true,  "pressure" },
            { ident.suctionFaceIndex,  false, "suction"  },
        };

        for (int si = 0; si < 2; ++si) {
            int fi = sides[si].faceIdx;
            if (fi < 0 || fi >= (int)r.faces.size()) continue;

            TopoDS_Face face = r.faces[fi];
            double v0, v1;
            std::string err;
            if (!sideVRange(face, sides[si].pressure, v0, v1, err)) continue;

            Handle(Geom_BSplineSurface) sf = faceToSurf(face);
            BRepAdaptor_Surface ad(face, true);
            SurfaceWrapper sw(sf, ad.FirstUParameter(), ad.LastUParameter(), v0, v1, false);

            RuledResult rr = fitRuledSegments(sw, 3, ParamDir::V, dd, nU, nV, 0, sides[si].prefix);
            std::string pfx = stem + "_" + sides[si].prefix;
            for (size_t k = 0; k < rr.segments.size(); ++k) {
                const auto& seg = rr.segments[k];
                exportOBJ(outDir + "/" + pfx + "_seg" + std::to_string(k) + ".obj",
                          seg.ruledMeshVerts, seg.ruledMeshFaces);
                exportDirectrixTXT(outDir + "/" + pfx + "_seg" + std::to_string(k) + "_params.txt",
                                   seg.curveC0Samples, seg.curveC1Samples);
            }
        }
        ++processed;
    }

    if (processed == 0) {
        res->errorCode = RULED_ERR_FILE_NOT_FOUND;
        return res;
    }

    res->errorCode = RULED_OK;
    res->numSurfaces = 0;
    std::ostringstream meta;
    meta << "{\"mode\":\"ruled-simple\",\"filesProcessed\":" << processed << "}";
    std::strncpy(res->metaJson, meta.str().c_str(), sizeof(res->metaJson) - 1);
    return res;
}

RULED_API void free_result(RuledFittingResult* result) { delete result; }

} // extern "C"
