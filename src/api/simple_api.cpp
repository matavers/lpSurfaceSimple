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

ParamDir toDir(RuledDirection d) { return (d == RULED_DIR_U) ? ParamDir::U : ParamDir::V; }

void exportCurveJSON(const std::string& path, const Vec3Arr& c0, const Vec3Arr& c1) {
    std::ofstream o(path);
    if (!o) return;
    o << "{\"C0\":["; for (size_t i = 0; i < c0.size(); ++i) {
        if (i) o << ","; o << "[" << c0[i].x() << "," << c0[i].y() << "," << c0[i].z() << "]"; }
    o << "],\"C1\":["; for (size_t i = 0; i < c1.size(); ++i) {
        if (i) o << ","; o << "[" << c1[i].x() << "," << c1[i].y() << "," << c1[i].z() << "]"; }
    o << "]}\n";
}

void exportPlaneTXT(const std::string& path, const Vec3& centroid, const Vec3& normal,
                     const Vec3Arr& corners) {
    std::ofstream o(path);
    if (!o) return;
    o << std::fixed << std::setprecision(6);
    o << "centroid = " << centroid.x() << " " << centroid.y() << " " << centroid.z() << "\n";
    o << "normal = " << normal.x() << " " << normal.y() << " " << normal.z() << "\n";
    if (corners.size() >= 4) {
        for (size_t i = 0; i < 4; ++i)
            o << "corner" << i << " = " << corners[i].x() << " "
              << corners[i].y() << " " << corners[i].z() << "\n";
    }
}

void exportMeta(const std::string& path, const std::string& mode,
                const std::string& f1, const std::string& f2,
                const std::vector<double>& errors1, const std::vector<double>& errors2) {
    std::ofstream o(path);
    o << "{\"mode\":\"" << mode << "\",\"files\":[\"" << f1 << "\",\"" << f2 << "\"],\"surfaces\":[";
    for (int s = 0; s < 2; ++s) {
        if (s) o << ",";
        o << "{\"name\":\"" << ((s==0)?"Blade-1":"Blade-2") << "\",\"segments\":[";
        const auto& e = (s==0) ? errors1 : errors2;
        for (size_t j = 0; j < e.size(); j+=2) {
            if (j) o << ",";
            o << "{\"index\":" << (j/2) << ",\"maxErr\":" << e[j] << ",\"rmsErr\":" << e[j+1] << "}";
        }
        o << "]}";
    }
    o << "]}\n";
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
    for (int i = 0; i < 2; ++i) if (c.numDirectrixDirs[i] <= 0) c.numDirectrixDirs[i] = 3;
    std::filesystem::create_directories(std::filesystem::path(c.outputDir ? c.outputDir : "."));

    StepLoadResult r1 = loadStepFile(c.stepFile1), r2 = loadStepFile(c.stepFile2);
    if (!r1.loaded || !r2.loaded) { res->errorCode = RULED_ERR_STEP_READ_FAILED; return res; }

    TopoDS_Face f1 = pickFace(r1.faces, c.faceIdx[0]), f2 = pickFace(r2.faces, c.faceIdx[1]);
    if (f1.IsNull() || f2.IsNull()) { res->errorCode = RULED_ERR_NO_VALID_FACE; return res; }

    Handle(Geom_BSplineSurface) sf1 = faceToSurf(f1), sf2 = faceToSurf(f2);
    BRepAdaptor_Surface a1(f1, true), a2(f2, true);

    SurfaceWrapper sw1(sf1, a1.FirstUParameter(), a1.LastUParameter(), a1.FirstVParameter(), a1.LastVParameter());
    SurfaceWrapper sw2(sf2, a2.FirstUParameter(), a2.LastUParameter(), a2.FirstVParameter(), a2.LastVParameter());

    std::vector<ParamDir> dd1, dd2;
    for (int j = 0; j < c.numDirectrixDirs[0]; ++j) dd1.push_back(toDir((RuledDirection)c.directrixDirs[0][j]));
    for (int j = 0; j < c.numDirectrixDirs[1]; ++j) dd2.push_back(toDir((RuledDirection)c.directrixDirs[1][j]));

    int ns = c.numDirectrixDirs[0];
    auto rr1 = fitRuledSegments(sw1, ns, toDir(c.splitDir[0]), dd1, c.nUSamples, c.nVSamples, 0, "Blade-1");
    auto rr2 = fitRuledSegments(sw2, ns, toDir(c.splitDir[1]), dd2, c.nUSamples, c.nVSamples, 0, "Blade-2");

    std::string od = c.outputDir ? c.outputDir : ".";
    // mesh + seg OBJs
    for (int bi = 0; bi < 2; ++bi) {
        auto& f = (bi==0)?f1:f2; auto& a = (bi==0)?a1:a2; auto& rr = (bi==0)?rr1:rr2;
        std::string pfx = (bi==0)?"blade1":"blade2";
        Vec3Arr mv; std::vector<std::array<int,3>> mf;
        double defl = (a.FirstUParameter() + a.FirstVParameter()) * 0.005;
        if (defl < 0.1) defl = 0.1;
        generateFaceMesh(f, defl, mv, mf);
        exportOBJ(od + "/" + pfx + "_mesh.obj", mv, mf);
        for (auto& seg : rr.segments) {
            std::string sp = od + "/" + pfx + "_seg" + std::to_string(seg.segmentIndex) + ".obj";
            exportOBJ(sp, seg.ruledMeshVerts, seg.ruledMeshFaces);
            std::string cp = od + "/" + pfx + "_seg" + std::to_string(seg.segmentIndex) + "_params.txt";
            exportCurveParamsTXT(cp, seg.curveC0, seg.curveC1, c.nUSamples);
        }
    }
    std::vector<double> e1, e2;
    for (auto& s : rr1.segments) { e1.push_back(s.maxError); e1.push_back(s.rmsError); }
    for (auto& s : rr2.segments) { e2.push_back(s.maxError); e2.push_back(s.rmsError); }
    exportMeta(od + "/meta.json", "ruled", c.stepFile1, c.stepFile2, e1, e2);

    std::ostringstream o; o << "{\"ok\":true,\"mode\":\"ruled\",\"surfaces\":[{";
    o << "\"name\":\"Blade-1\",\"segments\":[";
    for (size_t j = 0; j < rr1.segments.size(); ++j) {
        if (j) o << ","; o << "{\"maxErr\":" << rr1.segments[j].maxError << ",\"rmsErr\":" << rr1.segments[j].rmsError << "}"; }
    o << "]},{"; o << "\"name\":\"Blade-2\",\"segments\":[";
    for (size_t j = 0; j < rr2.segments.size(); ++j) {
        if (j) o << ","; o << "{\"maxErr\":" << rr2.segments[j].maxError << ",\"rmsErr\":" << rr2.segments[j].rmsError << "}"; }
    o << "]}]}";
    std::string j = o.str(); std::strncpy(res->metaJson, j.c_str(), sizeof(res->metaJson)-1);

    res->errorCode = RULED_OK; res->numSurfaces = 2;
    for (int si = 0; si < 2; ++si) {
        auto& r = (si==0)?rr1:rr2;
        std::strncpy(res->surfaces[si].name, r.name.c_str(), 63);
        res->surfaces[si].numSegments = (int)r.segments.size();
        for (size_t k = 0; k < r.segments.size() && k < 10; ++k) {
            res->surfaces[si].segments[k].segmentIndex = (int)k;
            res->surfaces[si].segments[k].maxError = r.segments[k].maxError;
            res->surfaces[si].segments[k].rmsError = r.segments[k].rmsError;
        }
    }
    return res;
}

RULED_API RuledFittingResult* plane_surface_fitting(const PlanarConfig* cfg) {
    // ... existing implementation ... unchanged
    RuledFittingResult* res = new RuledFittingResult(); std::memset(res, 0, sizeof(RuledFittingResult));
    if (!cfg || !cfg->stepFile1 || !cfg->stepFile2) {
        res->errorCode = RULED_ERR_INVALID_PARAMS; return res;
    }
    PlanarConfig c = *cfg;
    if (c.nUSamples <= 0) c.nUSamples = 50;
    if (c.nVSamples <= 0) c.nVSamples = 10;
    std::filesystem::create_directories(std::filesystem::path(c.outputDir ? c.outputDir : "."));

    StepLoadResult r1 = loadStepFile(c.stepFile1), r2 = loadStepFile(c.stepFile2);
    if (!r1.loaded || !r2.loaded) { res->errorCode = RULED_ERR_STEP_READ_FAILED; return res; }

    TopoDS_Face f1 = pickFace(r1.faces, c.faceIdx[0]), f2 = pickFace(r2.faces, c.faceIdx[1]);
    if (f1.IsNull() || f2.IsNull()) { res->errorCode = RULED_ERR_NO_VALID_FACE; return res; }

    Handle(Geom_BSplineSurface) sf1 = faceToSurf(f1), sf2 = faceToSurf(f2);
    BRepAdaptor_Surface a1(f1, true), a2(f2, true);

    SurfaceWrapper sw1(sf1, a1.FirstUParameter(), a1.LastUParameter(), a1.FirstVParameter(), a1.LastVParameter());
    SurfaceWrapper sw2(sf2, a2.FirstUParameter(), a2.LastUParameter(), a2.FirstVParameter(), a2.LastVParameter());

    auto pr1 = fitPlanarSegments(sw1, 3, toDir(c.splitDir[0]), c.nUSamples, c.nVSamples, 0, "Blade-1");
    auto pr2 = fitPlanarSegments(sw2, 3, toDir(c.splitDir[1]), c.nUSamples, c.nVSamples, 0, "Blade-2");

    std::string od = c.outputDir ? c.outputDir : ".";
    for (int bi = 0; bi < 2; ++bi) {
        auto& f = (bi==0)?f1:f2; auto& a = (bi==0)?a1:a2; auto& pr = (bi==0)?pr1:pr2;
        std::string pfx = (bi==0)?"blade1":"blade2";
        Vec3Arr mv; std::vector<std::array<int,3>> mf;
        double defl = (a.FirstUParameter() + a.FirstVParameter()) * 0.005;
        if (defl < 0.1) defl = 0.1;
        generateFaceMesh(f, defl, mv, mf);
        exportOBJ(od + "/" + pfx + "_mesh.obj", mv, mf);
        for (auto& seg : pr.segments) {
            std::string sp = od + "/" + pfx + "_plane" + std::to_string(seg.segmentIndex) + ".obj";
            exportOBJ(sp, seg.meshVerts, seg.meshFaces);
            std::string dp = od + "/" + pfx + "_plane" + std::to_string(seg.segmentIndex) + "_desc.txt";
            exportPlaneTXT(dp, seg.centroid, seg.normal, seg.meshVerts);
        }
    }
    std::vector<double> e1, e2;
    for (auto& s : pr1.segments) { e1.push_back(s.maxError); e1.push_back(s.rmsError); }
    for (auto& s : pr2.segments) { e2.push_back(s.maxError); e2.push_back(s.rmsError); }
    exportMeta(od + "/meta.json", "planar", c.stepFile1, c.stepFile2, e1, e2);

    std::ostringstream o; o << "{\"ok\":true,\"mode\":\"planar\",\"surfaces\":[{";
    o << "\"name\":\"Blade-1\",\"segments\":[";
    for (size_t j = 0; j < pr1.segments.size(); ++j) {
        if (j) o << ","; o << "{\"maxErr\":" << pr1.segments[j].maxError << ",\"rmsErr\":" << pr1.segments[j].rmsError << "}"; }
    o << "]},{"; o << "\"name\":\"Blade-2\",\"segments\":[";
    for (size_t j = 0; j < pr2.segments.size(); ++j) {
        if (j) o << ","; o << "{\"maxErr\":" << pr2.segments[j].maxError << ",\"rmsErr\":" << pr2.segments[j].rmsError << "}"; }
    o << "]}]}";
    std::string j = o.str(); std::strncpy(res->metaJson, j.c_str(), sizeof(res->metaJson)-1);

    res->errorCode = RULED_OK; res->numSurfaces = 2;
    for (int si = 0; si < 2; ++si) {
        auto& p = (si==0)?pr1:pr2;
        std::strncpy(res->surfaces[si].name, p.name.c_str(), 63);
        res->surfaces[si].numSegments = (int)p.segments.size();
        for (size_t k = 0; k < p.segments.size() && k < 10; ++k) {
            res->surfaces[si].segments[k].segmentIndex = (int)k;
            res->surfaces[si].segments[k].maxError = p.segments[k].maxError;
            res->surfaces[si].segments[k].rmsError = p.segments[k].rmsError;
        }
    }
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

RuledFittingResult* runRuledOnSide(const char* filepath, bool wantPressure,
                                    const RuledConfig* cfg) {
    RuledFittingResult* res = new RuledFittingResult();
    std::memset(res, 0, sizeof(RuledFittingResult));

    if (!cfg || !cfg->outputDir) {
        res->errorCode = RULED_ERR_INVALID_PARAMS; return res;
    }

    RuledConfig c = *cfg;
    if (c.nUSamples <= 0) c.nUSamples = 50;
    if (c.nVSamples <= 0) c.nVSamples = 10;
    if (c.nRibs <= 0) c.nRibs = 20;
    if (c.lambda <= 0.0) c.lambda = 1.0;
    if (c.numDirectrixDirs[0] <= 0) c.numDirectrixDirs[0] = 3;
    if (c.numDirectrixDirs[1] <= 0) c.numDirectrixDirs[1] = 3;

    std::filesystem::create_directories(std::filesystem::path(c.outputDir));

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

    SurfaceWrapper sw(sf, ad.FirstUParameter(), ad.LastUParameter(), v0, v1,
                      false);

    std::vector<ParamDir> dd;
    for (int j = 0; j < c.numDirectrixDirs[0]; ++j)
        dd.push_back(toDir((RuledDirection)c.directrixDirs[0][j]));

    int ns = c.numDirectrixDirs[0];
    const char* sideName = wantPressure ? "Pressure" : "Suction";
    auto rr = fitRuledSegments(sw, ns, toDir(c.splitDir[0]), dd, c.nUSamples,
                                c.nVSamples, 0, sideName);

    std::string od = c.outputDir;
    std::string pfx = wantPressure ? "pressure" : "suction";

    Vec3Arr mv; std::vector<std::array<int,3>> mf;
    double defl = std::max(0.2, (ad.LastUParameter() - ad.FirstUParameter()
                               + ad.LastVParameter() - ad.FirstVParameter()) * 0.01);
    generateFaceMesh(face, defl, mv, mf);
    exportOBJ(od + "/" + pfx + "_mesh.obj", mv, mf);

    for (auto& seg : rr.segments) {
        std::string sp = od + "/" + pfx + "_seg" + std::to_string(seg.segmentIndex) + ".obj";
        exportOBJ(sp, seg.ruledMeshVerts, seg.ruledMeshFaces);
        std::string cp = od + "/" + pfx + "_seg" + std::to_string(seg.segmentIndex) + "_params.txt";
        exportCurveParamsTXT(cp, seg.curveC0, seg.curveC1, c.nUSamples);
    }

    std::vector<double> e1, e2;
    for (auto& s : rr.segments) { e1.push_back(s.maxError); e1.push_back(s.rmsError); }
    exportMeta(od + "/meta.json", "ruled", filepath, "", e1, e2);

    std::ostringstream o;
    o << "{\"ok\":true,\"mode\":\"ruled\",\"side\":\"" << sideName
      << "\",\"surfaces\":[{\"name\":\"" << sideName << "\",\"segments\":[";
    for (size_t j = 0; j < rr.segments.size(); ++j) {
        if (j) o << ",";
        o << "{\"maxErr\":" << rr.segments[j].maxError << ",\"rmsErr\":" << rr.segments[j].rmsError << "}";
    }
    o << "]}]}";
    std::string j = o.str();
    std::strncpy(res->metaJson, j.c_str(), sizeof(res->metaJson)-1);

    res->errorCode = RULED_OK; res->numSurfaces = 1;
    std::strncpy(res->surfaces[0].name, sideName, 63);
    res->surfaces[0].numSegments = (int)rr.segments.size();
    for (size_t k = 0; k < rr.segments.size() && k < 10; ++k) {
        res->surfaces[0].segments[k].segmentIndex = (int)k;
        res->surfaces[0].segments[k].maxError = rr.segments[k].maxError;
        res->surfaces[0].segments[k].rmsError = rr.segments[k].rmsError;
    }
    return res;
}

RuledFittingResult* runPlanarOnSide(const char* filepath, bool wantPressure,
                                     const PlanarConfig* cfg) {
    RuledFittingResult* res = new RuledFittingResult();
    std::memset(res, 0, sizeof(RuledFittingResult));

    if (!cfg || !cfg->outputDir) {
        res->errorCode = RULED_ERR_INVALID_PARAMS; return res;
    }

    PlanarConfig c = *cfg;
    if (c.nUSamples <= 0) c.nUSamples = 50;
    if (c.nVSamples <= 0) c.nVSamples = 10;
    std::filesystem::create_directories(std::filesystem::path(c.outputDir));

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
    auto pr = fitPlanarSegments(sw, 3, toDir(c.splitDir[0]), c.nUSamples, c.nVSamples, 0, sideName);

    std::string od = c.outputDir;
    std::string pfx = wantPressure ? "pressure" : "suction";

    Vec3Arr mv; std::vector<std::array<int,3>> mf;
    double defl = std::max(0.2, (ad.LastUParameter() - ad.FirstUParameter()
                               + ad.LastVParameter() - ad.FirstVParameter()) * 0.01);
    generateFaceMesh(face, defl, mv, mf);
    exportOBJ(od + "/" + pfx + "_mesh.obj", mv, mf);

    for (auto& seg : pr.segments) {
        std::string sp = od + "/" + pfx + "_plane" + std::to_string(seg.segmentIndex) + ".obj";
        exportOBJ(sp, seg.meshVerts, seg.meshFaces);
        std::string dp = od + "/" + pfx + "_plane" + std::to_string(seg.segmentIndex) + "_desc.txt";
        exportPlaneTXT(dp, seg.centroid, seg.normal, seg.meshVerts);
    }

    std::vector<double> e1, e2;
    for (auto& s : pr.segments) { e1.push_back(s.maxError); e1.push_back(s.rmsError); }
    exportMeta(od + "/meta.json", "planar", filepath, "", e1, e2);

    std::ostringstream o;
    o << "{\"ok\":true,\"mode\":\"planar\",\"side\":\"" << sideName
      << "\",\"surfaces\":[{\"name\":\"" << sideName << "\",\"segments\":[";
    for (size_t j = 0; j < pr.segments.size(); ++j) {
        if (j) o << ",";
        o << "{\"maxErr\":" << pr.segments[j].maxError << ",\"rmsErr\":" << pr.segments[j].rmsError << "}";
    }
    o << "]}]}";
    std::string j = o.str();
    std::strncpy(res->metaJson, j.c_str(), sizeof(res->metaJson)-1);

    res->errorCode = RULED_OK; res->numSurfaces = 1;
    std::strncpy(res->surfaces[0].name, sideName, 63);
    res->surfaces[0].numSegments = (int)pr.segments.size();
    for (size_t k = 0; k < pr.segments.size() && k < 10; ++k) {
        res->surfaces[0].segments[k].segmentIndex = (int)k;
        res->surfaces[0].segments[k].maxError = pr.segments[k].maxError;
        res->surfaces[0].segments[k].rmsError = pr.segments[k].rmsError;
    }
    return res;
}

} // anon

RULED_API RuledFittingResult* pressure_ruled_fitting(const RuledConfig* cfg) {
    return runRuledOnSide(cfg ? cfg->stepFile1 : nullptr, true, cfg);
}
RULED_API RuledFittingResult* pressure_plane_fitting(const PlanarConfig* cfg) {
    return runPlanarOnSide(cfg ? cfg->stepFile1 : nullptr, true, cfg);
}
RULED_API RuledFittingResult* suction_ruled_fitting(const RuledConfig* cfg) {
    return runRuledOnSide(cfg ? cfg->stepFile1 : nullptr, false, cfg);
}
RULED_API RuledFittingResult* suction_plane_fitting(const PlanarConfig* cfg) {
    return runPlanarOnSide(cfg ? cfg->stepFile1 : nullptr, false, cfg);
}

RULED_API void free_result(RuledFittingResult* result) { delete result; }

} // extern "C"
