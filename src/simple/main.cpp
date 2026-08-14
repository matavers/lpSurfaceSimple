#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>
#include <filesystem>

#include "simple/common.hpp"
#include "simple/step_reader.hpp"
#include "simple/surface_wrapper.hpp"
#include "simple/ruled_fitter.hpp"
#include "simple/planar_fitter.hpp"
#include "simple/grid_fitter.hpp"
#include "simple/blend_fitter.hpp"
#include "simple/blade_identifier.hpp"
#include "simple/blade_splitter.hpp"

#include <IGESControl_Writer.hxx>
#include <BRep_Builder.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRep_Tool.hxx>
#include <GeomConvert.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>

using namespace simple;
namespace fs = std::filesystem;

static std::string jsonSafe(const std::string& s) {
    std::string r = s;
    for (char& c : r) {
        if (c == '\\') c = '/';
        if (c == '"') c = '\'';
    }
    return r;
}

TopoDS_Face pickLargestFace(
    const std::vector<TopoDS_Face>& faces,
    const std::string& label)
{
    if (faces.empty()) return TopoDS_Face();
    if (faces.size() == 1) return faces[0];

    int bestIdx = 0;
    double bestDiag = 0.0;

    std::cout << "  " << label << ": " << faces.size() << " faces, picking largest:" << std::endl;

    for (size_t i = 0; i < faces.size(); ++i) {
        BRepAdaptor_Surface adapt(faces[i], true);
        double umin = adapt.FirstUParameter(), umax = adapt.LastUParameter();
        double vmin = adapt.FirstVParameter(), vmax = adapt.LastVParameter();
        int nWires = 0;
        for (TopExp_Explorer we(faces[i], TopAbs_WIRE); we.More(); we.Next()) ++nWires;

        Vec3 bbMin(1e30, 1e30, 1e30), bbMax(-1e30, -1e30, -1e30);
        for (int ui = 0; ui <= 4; ++ui) {
            for (int vi = 0; vi <= 4; ++vi) {
                double u = umin + (umax - umin) * ui / 4.0;
                double v = vmin + (vmax - vmin) * vi / 4.0;
                gp_Pnt p = adapt.Value(u, v);
                bbMin = Vec3(std::min(bbMin.x(), p.X()), std::min(bbMin.y(), p.Y()), std::min(bbMin.z(), p.Z()));
                bbMax = Vec3(std::max(bbMax.x(), p.X()), std::max(bbMax.y(), p.Y()), std::max(bbMax.z(), p.Z()));
            }
        }
        double diag = (bbMax - bbMin).norm();
        std::cout << "    face[" << i << "] type=" << adapt.GetType()
                  << " wires=" << nWires
                  << " U=[" << umin << "," << umax << "] V=[" << vmin << "," << vmax << "]"
                  << " bbox=" << std::fixed << std::setprecision(3) << diag;

        if (diag > bestDiag) {
            bestDiag = diag;
            bestIdx = static_cast<int>(i);
            std::cout << "  <-- selected";
        }
        std::cout << std::endl;
    }
    return faces[bestIdx];
}

void printUsage() {
    std::cout << "Usage: simple.exe <step_file1> <step_file2> [options]\n"
              << "  step_file1, step_file2 : Path to STEP blade model files\n"
              << "  Options:\n"
              << "    --mode <ruled|planar|info|face-obj|extract-faces|auto-identify|split-blade>  Algorithm mode (default: ruled)\n"
              << "    --outdir <DIR>          Output directory (default: ./output)\n"
              << "    --face-idx1 <N>         Face index for file 1 (default: auto-pick largest)\n"
              << "    --face-idx2 <N>         Face index for file 2 (default: auto-pick largest)\n"
              << "    --nusamples <N>         U-direction samples (default: 50)\n"
              << "    --nvsamples <N>         V-direction samples (default: 10)\n"
              << "    --nsplit-u <N>          U-direction (vertical) splits -> columns (default: 2)\n"
              << "    --nsplit-v <N>          V-direction (horizontal) splits -> rows (default: 2)\n"
              << "    --tolerance <T>         Tolerance for adaptive refinement (default: 0.1)\n"
              << "    --max-depth <N>         Max refinement steps (default: 20)\n"
              << "    --blend                 Enable trim+blend post-pass (ruled mode only)\n"
              << "    --fmax <F>              Max per-cell trim fraction for blend (default: 0.3)\n"
              << "    --help                  Show this help\n"
              << std::endl;
}

int main(int argc, char* argv[]) {
    std::string stepFile1, stepFile2, outDir = "output";
    std::string modeStr = "ruled";
    int nUSamples = 50, nVSamples = 10;
    int nSplitU = 2, nSplitV = 2, maxDepth = 20;
    double tolerance = 0.1;
    bool doBlend = false;
    double fMax = 0.3;
    int faceIdx1 = -1, faceIdx2 = -1;
    std::string faceOutPath;
    double uRange1Min = -1, uRange1Max = -1;
    double uRange2Min = -1, uRange2Max = -1;
    double vRange1Min = -1, vRange1Max = -1;
    double vRange2Min = -1, vRange2Max = -1;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--help" || arg == "-h") { printUsage(); return 0; }
        else if (arg == "--mode" && i + 1 < argc) { modeStr = argv[++i]; }
        else if (arg == "--outdir" && i + 1 < argc) { outDir = argv[++i]; }
        else if (arg == "--face-idx1" && i + 1 < argc) { faceIdx1 = std::stoi(argv[++i]); }
        else if (arg == "--face-idx2" && i + 1 < argc) { faceIdx2 = std::stoi(argv[++i]); }
        else if (arg == "--face-idx" && i + 1 < argc) { faceIdx1 = std::stoi(argv[++i]); }
        else if (arg == "--face-out" && i + 1 < argc) { faceOutPath = argv[++i]; }
        else if (arg == "--extract-out" && i + 1 < argc) { faceOutPath = argv[++i]; }
        else if (arg == "--u-range1" && i + 1 < argc) {
            std::string s(argv[++i]);
            auto comma = s.find(',');
            if (comma != std::string::npos) {
                uRange1Min = std::stod(s.substr(0, comma));
                uRange1Max = std::stod(s.substr(comma + 1));
            }
        }
        else if (arg == "--u-range2" && i + 1 < argc) {
            std::string s(argv[++i]);
            auto comma = s.find(',');
            if (comma != std::string::npos) {
                uRange2Min = std::stod(s.substr(0, comma));
                uRange2Max = std::stod(s.substr(comma + 1));
            }
        }
        else if (arg == "--v-range1" && i + 1 < argc) {
            std::string s(argv[++i]);
            auto comma = s.find(',');
            if (comma != std::string::npos) {
                vRange1Min = std::stod(s.substr(0, comma));
                vRange1Max = std::stod(s.substr(comma + 1));
            }
        }
        else if (arg == "--v-range2" && i + 1 < argc) {
            std::string s(argv[++i]);
            auto comma = s.find(',');
            if (comma != std::string::npos) {
                vRange2Min = std::stod(s.substr(0, comma));
                vRange2Max = std::stod(s.substr(comma + 1));
            }
        }
        else if (arg == "--nusamples" && i + 1 < argc) { nUSamples = std::stoi(argv[++i]); }
        else if (arg == "--nvsamples" && i + 1 < argc) { nVSamples = std::stoi(argv[++i]); }
        else if (arg == "--nsplit-u" && i + 1 < argc) { nSplitU = std::stoi(argv[++i]); }
        else if (arg == "--nsplit-v" && i + 1 < argc) { nSplitV = std::stoi(argv[++i]); }
        else if (arg == "--tolerance" && i + 1 < argc) { tolerance = std::stod(argv[++i]); }
        else if (arg == "--max-depth" && i + 1 < argc) { maxDepth = std::stoi(argv[++i]); }
        else if (arg == "--blend") { doBlend = true; }
        else if (arg == "--fmax" && i + 1 < argc) { fMax = std::stod(argv[++i]); }
        else if (stepFile1.empty()) { stepFile1 = arg; }
        else if (stepFile2.empty()) { stepFile2 = arg; }
    }

    bool isPlanar = (modeStr == "planar");

    if (stepFile1.empty() || stepFile2.empty()) {
        // --info mode
        if (!stepFile1.empty() && modeStr == "info") {
            auto r = simple::loadStepFile(stepFile1);
            if (!r.loaded) {
                std::cerr << "[Error] " << r.errorMsg << std::endl;
                return 1;
            }
            std::cout << "{\"file\":\"" << r.filename << "\",\"numFaces\":" << r.faces.size() << ",\"faces\":[";
            for (size_t i = 0; i < r.faces.size(); ++i) {
                if (i > 0) std::cout << ",";
                BRepAdaptor_Surface adapt(r.faces[i], true);
                Bnd_Box bbox;
                BRepBndLib::Add(r.faces[i], bbox);
                double x1,y1,z1,x2,y2,z2;
                bbox.Get(x1,y1,z1,x2,y2,z2);
                double diag = std::sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
                double area = (adapt.LastUParameter()-adapt.FirstUParameter())
                            * (adapt.LastVParameter()-adapt.FirstVParameter());
                int nWires = 0;
                for (TopExp_Explorer we(r.faces[i], TopAbs_WIRE); we.More(); we.Next()) ++nWires;
                std::cout << "{\"index\":" << i
                          << ",\"type\":" << (int)adapt.GetType()
                          << ",\"uMin\":" << adapt.FirstUParameter()
                          << ",\"uMax\":" << adapt.LastUParameter()
                          << ",\"vMin\":" << adapt.FirstVParameter()
                          << ",\"vMax\":" << adapt.LastVParameter()
                          << ",\"diag\":" << diag
                          << ",\"area\":" << area
                          << ",\"nWires\":" << nWires
                          << "}";
            }
            std::cout << "]}" << std::endl;
            return 0;
        }
        // --mode face-obj: export single face mesh as OBJ
        if (!stepFile1.empty() && modeStr == "face-obj" && !faceOutPath.empty()) {
            auto r = simple::loadStepFile(stepFile1);
            if (!r.loaded) {
                std::cerr << "[Error] " << r.errorMsg << std::endl;
                return 1;
            }
            int fi = (faceIdx1 >= 0 && faceIdx1 < (int)r.faces.size()) ? faceIdx1 : 0;
            BRepAdaptor_Surface ad(r.faces[fi], true);
            Bnd_Box bbox;
            BRepBndLib::Add(r.faces[fi], bbox);
            double x1, y1, z1, x2, y2, z2;
            bbox.Get(x1, y1, z1, x2, y2, z2);
            double diag3d = std::sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
            double defl = std::max(0.02, diag3d * 0.0005);
            Vec3Arr mv; std::vector<std::array<int,3>> mf;
            generateFaceMesh(r.faces[fi], defl, mv, mf);
            exportOBJ(faceOutPath, mv, mf);
            std::cout << "wrote " << faceOutPath << std::endl;
            return 0;
        }
        // --mode extract-faces: export selected faces as clean IGES
        if (!stepFile1.empty() && modeStr == "extract-faces" && !faceOutPath.empty()) {
            auto r = simple::loadStepFile(stepFile1);
            if (!r.loaded) { std::cerr << "[Error] " << r.errorMsg << std::endl; return 1; }
            int fi1 = (faceIdx1 >= 0 && faceIdx1 < (int)r.faces.size()) ? faceIdx1 : 0;
            int fi2 = (faceIdx2 >= 0 && faceIdx2 < (int)r.faces.size()) ? faceIdx2 : 1;
            if (fi1 == fi2 && r.faces.size() > 1) fi2 = (fi1 + 1) % r.faces.size();

            TopoDS_Compound comp;
            BRep_Builder b; b.MakeCompound(comp);
            b.Add(comp, r.faces[fi1]);
            b.Add(comp, r.faces[fi2]);

            IGESControl_Writer writer;
            writer.AddShape(comp);
            writer.Write(faceOutPath.c_str());
            std::cout << "wrote " << faceOutPath << std::endl;
            return 0;
        }
        // --mode auto-identify: automatically detect pressure/suction faces
        if (!stepFile1.empty() && modeStr == "auto-identify") {
            auto r = simple::loadStepFile(stepFile1);
            if (!r.loaded) { std::cerr << "[Error] " << r.errorMsg << std::endl; return 1; }
            auto ident = simple::identifyBladeSurfaces(r.faces, 9, 0.1, 20.0, 5);
            std::string escapedMsg = ident.message;
            for (auto& c : escapedMsg) {
                if (c == '\n') c = ' ';
                if (c == '\r') c = ' ';
                if (c == '"' || c == '\\') c = '\'';
            }
            double pDiag = 0, sDiag = 0;
            if (ident.pressureFaceIndex >= 0 && ident.pressureFaceIndex < (int)r.faces.size()) {
                Bnd_Box bb; BRepBndLib::Add(r.faces[ident.pressureFaceIndex], bb);
                double x1,y1,z1,x2,y2,z2; bb.Get(x1,y1,z1,x2,y2,z2);
                pDiag = std::sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
            }
            if (ident.suctionFaceIndex >= 0 && ident.suctionFaceIndex < (int)r.faces.size()) {
                Bnd_Box bb; BRepBndLib::Add(r.faces[ident.suctionFaceIndex], bb);
                double x1,y1,z1,x2,y2,z2; bb.Get(x1,y1,z1,x2,y2,z2);
                sDiag = std::sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
            }
            std::cout << "{\"file\":\"" << r.filename
                      << "\",\"numFaces\":" << r.faces.size()
                      << ",\"success\":" << (ident.success ? "true" : "false")
                      << ",\"pressureIndex\":" << ident.pressureFaceIndex
                      << ",\"suctionIndex\":" << ident.suctionFaceIndex
                      << ",\"pressureDiag\":" << pDiag
                      << ",\"suctionDiag\":" << sDiag
                      << ",\"message\":\"" << escapedMsg << "\"}" << std::endl;
            return ident.success ? 0 : 1;
        }
        // --mode split-blade: detect curvature peaks and split into dynamic regions
        if (!stepFile1.empty() && modeStr == "split-blade") {
            auto r = simple::loadStepFile(stepFile1);
            if (!r.loaded) { std::cerr << "[Error] " << r.errorMsg << std::endl; return 1; }
            int fi = (faceIdx1 >= 0 && faceIdx1 < (int)r.faces.size()) ? faceIdx1 : 0;
            auto sp = simple::splitBladeFaceBySection(r.faces[fi]);
            std::string escapedMsg = sp.message;
            for (auto& c : escapedMsg) {
                if (c == '\n') c = ' ';
                if (c == '\r') c = ' ';
                if (c == '"' || c == '\\') c = '\'';
            }
            std::cout << "{\"faceIndex\":" << fi
                      << ",\"success\":" << (sp.success ? "true" : "false")
                      << ",\"dir\":\"" << sp.splitDir << "\""
                      << ",\"regions\":[";
            for (size_t i = 0; i < sp.regions.size(); ++i) {
                if (i > 0) std::cout << ",";
                std::cout << "{\"uStart\":" << sp.regions[i].uStart
                          << ",\"uEnd\":" << sp.regions[i].uEnd
                          << ",\"avgCurv\":" << sp.regions[i].avgCurv
                          << ",\"label\":\"" << sp.regions[i].label << "\"}";
            }
            std::cout << "],\"message\":\"" << escapedMsg << "\"}" << std::endl;
            return sp.success ? 0 : 1;
        }
        // --mode face-obj-range: export partial face mesh (U-range or V-range restriction)
        if (!stepFile1.empty() && modeStr == "face-obj-range" && !faceOutPath.empty()) {
            auto r = simple::loadStepFile(stepFile1);
            if (!r.loaded) { std::cerr << "[Error] " << r.errorMsg << std::endl; return 1; }
            int fi = (faceIdx1 >= 0 && faceIdx1 < (int)r.faces.size()) ? faceIdx1 : 0;
            BRepAdaptor_Surface ad(r.faces[fi], true);
            if (ad.GetType() != GeomAbs_BSplineSurface) { std::cerr << "[Error] Face not BSpline" << std::endl; return 1; }
            Handle(Geom_BSplineSurface) s = ad.BSpline();
            if (s.IsNull()) s = Handle(Geom_BSplineSurface)::DownCast(BRep_Tool::Surface(r.faces[fi]));

            double us = ad.FirstUParameter(), ue = ad.LastUParameter();
            double vs = ad.FirstVParameter(), ve = ad.LastVParameter();
            bool restrictU = (uRange1Min >= 0);
            bool restrictV = (vRange1Min >= 0);

            if (restrictU) { us = uRange1Min; ue = uRange1Max; if (us > ue) std::swap(us, ue); }
            if (restrictV) { vs = vRange1Min; ve = vRange1Max; if (vs > ve) std::swap(vs, ve); }

            int ru = restrictU ? 30 : 60;
            int rv = restrictV ? 30 : 15;
            Vec3Arr verts; FaceArr faces;
            for (int i = 0; i <= ru; ++i) {
                double u = us + (ue - us) * i / ru;
                for (int j = 0; j <= rv; ++j) {
                    double v = vs + (ve - vs) * j / rv;
                    gp_Pnt p = s->Value(u, v);
                    verts.push_back(Vec3(p.X(), p.Y(), p.Z()));
                }
            }
            for (int i = 0; i < ru; ++i) {
                for (int j = 0; j < rv; ++j) {
                    int a = i * (rv + 1) + j;
                    faces.push_back({a, a + 1, a + rv + 1});
                    faces.push_back({a + 1, a + rv + 2, a + rv + 1});
                }
            }
            exportOBJ(faceOutPath, verts, faces);
            std::cout << "wrote " << faceOutPath << std::endl;
            return 0;
        }
        // --mode iso-curve: extract iso-parametric curve as OBJ polyline
        if (!stepFile1.empty() && modeStr == "iso-curve" && !faceOutPath.empty()) {
            auto r = simple::loadStepFile(stepFile1);
            if (!r.loaded) { std::cerr << "[Error] " << r.errorMsg << std::endl; return 1; }
            int fi = (faceIdx1 >= 0 && faceIdx1 < (int)r.faces.size()) ? faceIdx1 : 0;
            BRepAdaptor_Surface ad(r.faces[fi], true);
            Handle(Geom_BSplineSurface) s = ad.BSpline();
            if (s.IsNull()) s = Handle(Geom_BSplineSurface)::DownCast(BRep_Tool::Surface(r.faces[fi]));
            if (s.IsNull()) { std::cerr << "[Error] No BSpline" << std::endl; return 1; }
            double paramF = (uRange1Min >= 0) ? uRange1Min : 0.5;
            bool useVIso = (vRange1Min >= 0);
            double vParam = useVIso ? vRange1Min : 0.5;
            Handle(Geom_Curve) curve = useVIso ? s->VIso(vParam) : s->UIso(paramF);
            if (curve.IsNull()) { std::cerr << "[Error] Null iso-curve" << std::endl; return 1; }
            GeomAdaptor_Curve gac(curve);
            double t0 = gac.FirstParameter(), t1 = gac.LastParameter();
            int nP = 200;
            std::ofstream out(faceOutPath);
            out << std::fixed << std::setprecision(6);
            for (int i = 0; i < nP; ++i) {
                double t = t0 + (t1 - t0) * i / (nP - 1);
                gp_Pnt p; gac.D0(t, p);
                out << "v " << p.X() << " " << p.Y() << " " << p.Z() << "\n";
            }
            for (int i = 0; i < nP - 1; ++i)
                out << "l " << (i+1) << " " << (i+2) << "\n";
            std::cout << "wrote " << faceOutPath << std::endl;
            return 0;
        }
        // allow single-file: if it has 2 faces, use them directly
        if (!stepFile1.empty() && stepFile2.empty()) {
            stepFile2 = stepFile1;
        } else if (stepFile1.empty() || stepFile2.empty()) {
            std::cerr << "[Error] At least one file required.\n";
            printUsage();
            return 1;
        }
    }

    fs::create_directories(outDir);

    std::cout << "[Step 1] Loading STEP files..." << std::endl;

    auto result1 = simple::loadStepFile(stepFile1);
    auto result2 = simple::loadStepFile(stepFile2);

    if (!result1.loaded) {
        std::cerr << "[Error] " << result1.errorMsg << std::endl;
        return 1;
    }
    if (!result2.loaded) {
        std::cerr << "[Error] " << result2.errorMsg << std::endl;
        return 1;
    }

    std::cout << "  File 1: " << result1.filename
              << " - " << result1.faces.size() << " face(s)" << std::endl;
    std::cout << "  File 2: " << result2.filename
              << " - " << result2.faces.size() << " face(s)" << std::endl;

    TopoDS_Face face1, face2;
    if (faceIdx1 >= 0 && faceIdx1 < (int)result1.faces.size()) {
        face1 = result1.faces[faceIdx1];
        std::cout << "  Using face[" << faceIdx1 << "] from file 1" << std::endl;
    } else {
        face1 = pickLargestFace(result1.faces, "Blade-raw1");
    }
    if (faceIdx2 >= 0 && faceIdx2 < (int)result2.faces.size()) {
        face2 = result2.faces[faceIdx2];
        std::cout << "  Using face[" << faceIdx2 << "] from file 2" << std::endl;
    } else {
        face2 = pickLargestFace(result2.faces, "Blade-raw2");
    }

    if (face1.IsNull() || face2.IsNull()) {
        std::cerr << "[Error] No valid face found.\n";
        return 1;
    }

    std::cout << "[Step 2] Creating surface wrappers..." << std::endl;

    BRepAdaptor_Surface adapt1(face1, true);
    BRepAdaptor_Surface adapt2(face2, true);

    Handle(Geom_BSplineSurface) surf1;
    Handle(Geom_BSplineSurface) surf2;

    if (adapt1.GetType() == GeomAbs_BSplineSurface) {
        surf1 = adapt1.BSpline();
    } else {
        Handle(Geom_Surface) gs = BRep_Tool::Surface(face1);
        surf1 = GeomConvert::SurfaceToBSplineSurface(gs);
    }
    if (adapt2.GetType() == GeomAbs_BSplineSurface) {
        surf2 = adapt2.BSpline();
    } else {
        Handle(Geom_Surface) gs = BRep_Tool::Surface(face2);
        surf2 = GeomConvert::SurfaceToBSplineSurface(gs);
    }

    double u1Min = adapt1.FirstUParameter(), u1Max = adapt1.LastUParameter();
    double v1Min = adapt1.FirstVParameter(), v1Max = adapt1.LastVParameter();
    double u2Min = adapt2.FirstUParameter(), u2Max = adapt2.LastUParameter();
    double v2Min = adapt2.FirstVParameter(), v2Max = adapt2.LastVParameter();

    bool wrap1 = (uRange1Min >= 0 && uRange1Max < uRange1Min);
    bool wrap2 = (uRange2Min >= 0 && uRange2Max < uRange2Min);
    double u1m = (uRange1Min >= 0) ? uRange1Min : u1Min;
    double u1M = (uRange1Min >= 0) ? uRange1Max : u1Max;
    double u2m = (uRange2Min >= 0) ? uRange2Min : u2Min;
    double u2M = (uRange2Min >= 0) ? uRange2Max : u2Max;
    double v1m = (vRange1Min >= 0) ? vRange1Min : v1Min;
    double v1M = (vRange1Min >= 0) ? vRange1Max : v1Max;
    double v2m = (vRange2Min >= 0) ? vRange2Min : v2Min;
    double v2M = (vRange2Min >= 0) ? vRange2Max : v2Max;

    simple::SurfaceWrapper sw1(surf1, u1m, u1M, v1m, v1M, wrap1);
    simple::SurfaceWrapper sw2(surf2, u2m, u2M, v2m, v2M, wrap2);

    auto [u1min, u1max] = sw1.paramDomainU();
    auto [v1min, v1max] = sw1.paramDomainV();
    auto [u2min, u2max] = sw2.paramDomainU();
    auto [v2min, v2max] = sw2.paramDomainV();

    std::cout << "  Surface 1: U=[" << u1min << "," << u1max << "] V=["
              << v1min << "," << v1max << "] deg=("
              << surf1->UDegree() << "," << surf1->VDegree() << ")" << std::endl;
    std::cout << "  Surface 2: U=[" << u2min << "," << u2max << "] V=["
              << v2min << "," << v2max << "] deg=("
              << surf2->UDegree() << "," << surf2->VDegree() << ")" << std::endl;

    std::cout << "[Step 3] Fitting " << (isPlanar ? "planar" : "ruled")
              << " grid with adaptive refinement..." << std::endl;
    std::cout << "  Splits: U=" << nSplitU << " V=" << nSplitV
              << "  U-samples: " << nUSamples
              << "  V-samples: " << nVSamples
              << "  tolerance: " << tolerance
              << "  max-depth: " << maxDepth << std::endl;

    GridConfig gcfg;
    gcfg.nSplitU = nSplitU;
    gcfg.nSplitV = nSplitV;
    gcfg.tolerance = tolerance;
    gcfg.maxDepth = maxDepth;
    gcfg.nUSamples = nUSamples;
    gcfg.nVSamples = nVSamples;

    GridResult gr1, gr2;
    if (isPlanar) {
        gr1 = fitGridPlanar(sw1, gcfg, "Blade-1");
        gr2 = fitGridPlanar(sw2, gcfg, "Blade-2");
    } else {
        gr1 = fitGridRuled(sw1, gcfg, "Blade-1");
        gr2 = fitGridRuled(sw2, gcfg, "Blade-2");
    }

    for (const auto& gr : {gr1, gr2}) {
        std::cout << "  " << gr.name << ": " << gr.nRows << "x" << gr.nCols
                  << " = " << gr.cells.size() << " cells" << std::endl;
        for (const auto& c : gr.cells) {
            std::cout << "    cell[" << c.row << "," << c.col << "]"
                      << " dir=" << (c.fitDir == ParamDir::U ? "U" : "V")
                      << "  maxError=" << std::fixed << std::setprecision(5) << c.maxError
                      << "  rmsError=" << c.rmsError << std::endl;
        }
    }

    BlendResult br1, br2;
    if (doBlend && !isPlanar) {
        BlendConfig bcfg;
        bcfg.fMax = fMax;
        std::cout << "[Step 3b] Trim+blend post-pass (fMax=" << fMax << ")..." << std::endl;
        std::vector<double> f1 = optimizeTrim(sw1, gr1, bcfg);
        std::vector<double> f2 = optimizeTrim(sw2, gr2, bcfg);
        br1 = buildBlend(sw1, gr1, f1, bcfg);
        br2 = buildBlend(sw2, gr2, f2, bcfg);
        auto printF = [](const GridResult& gr, const std::vector<double>& f) {
            std::cout << "    trim f:";
            for (size_t i = 0; i < gr.cells.size(); ++i)
                std::cout << " [" << gr.cells[i].row << "," << gr.cells[i].col << "]="
                          << std::fixed << std::setprecision(3) << f[i];
            std::cout << std::endl;
        };
        printF(gr1, f1);
        printF(gr2, f2);
        for (const auto& br : {br1, br2}) {
            std::cout << "  blend: totalMaxError=" << std::fixed << std::setprecision(5)
                      << br.totalMaxError << "  totalRmsError=" << br.totalRmsError
                      << "  strips=" << br.strips.size() << std::endl;
        }
    }

    std::cout << "[Step 4] Generating trimmed face meshes..." << std::endl;
    Vec3Arr mv1, mv2; FaceArr mf1, mf2; Vec2Arr dummyUV;
    sw1.generateMesh(40, 20, mv1, mf1, dummyUV);
    if (!(faceIdx1 >= 0 && faceIdx2 >= 0 && faceIdx1 == faceIdx2)) {
        sw2.generateMesh(40, 20, mv2, mf2, dummyUV);
    } else {
        mv2 = mv1; mf2 = mf1;
    }

    std::cout << "[Step 5] Exporting data to: " << outDir << std::endl;
    exportOBJ(outDir + "/blade1_mesh.obj", mv1, mf1);
    std::cout << "  wrote " << outDir << "/blade1_mesh.obj" << std::endl;
    exportOBJ(outDir + "/blade2_mesh.obj", mv2, mf2);
    std::cout << "  wrote " << outDir << "/blade2_mesh.obj" << std::endl;

    exportGridOBJs(outDir, "blade1", gr1, isPlanar);
    exportGridOBJs(outDir, "blade2", gr2, isPlanar);
    exportGridLinesVTK(outDir + "/blade1_grid.vtk", sw1, gr1);
    exportGridLinesVTK(outDir + "/blade2_grid.vtk", sw2, gr2);

    if (doBlend && !isPlanar) {
        exportBlendOBJs(outDir, "blade1", br1);
        exportBlendOBJs(outDir, "blade2", br2);
        std::cout << "  wrote blend OBJs (trimmed cells + strips)" << std::endl;
    }

    {
        std::ofstream metaOut(outDir + "/meta.json");
        metaOut << buildGridMetaJson({gr1, gr2},
                                     isPlanar ? "planar" : "ruled",
                                     stepFile1, stepFile2) << "\n";
    }
    std::cout << "  wrote " << outDir << "/meta.json" << std::endl;

    std::cout << "\n[Done] All files exported to " << outDir << std::endl;
    return 0;
}
