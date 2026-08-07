#include "simple/step_reader.hpp"

#include <STEPControl_Reader.hxx>
#include <IGESControl_Reader.hxx>
#include <TopExp_Explorer.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <Poly_Triangulation.hxx>
#include <Poly_Triangle.hxx>
#include <ShapeFix_Shape.hxx>
#include <TopLoc_Location.hxx>

#include <algorithm>
#include <cstring>

namespace simple {

StepLoadResult loadStepFile(const std::string& filepath) {
    StepLoadResult result;
    result.filename = filepath;
    result.loaded = false;

    IFSelect_ReturnStatus status = IFSelect_RetVoid;
    TopoDS_Shape shape;

    std::string ext;
    auto dotPos = filepath.rfind('.');
    if (dotPos != std::string::npos)
        ext = filepath.substr(dotPos);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    if (ext == ".igs" || ext == ".iges") {
        IGESControl_Reader reader;
        status = reader.ReadFile(filepath.c_str());
        if (status == IFSelect_RetDone) {
            reader.TransferRoots();
            shape = reader.OneShape();
        }
    } else {
        STEPControl_Reader reader;
        status = reader.ReadFile(filepath.c_str());
        if (status == IFSelect_RetDone) {
            reader.TransferRoots();
            shape = reader.OneShape();
        }
    }

    if (status != IFSelect_RetDone) {
        result.errorMsg = "Failed to read file: " + filepath;
        return result;
    }
    if (shape.IsNull()) {
        result.errorMsg = "Null shape from file: " + filepath;
        return result;
    }

    Handle(ShapeFix_Shape) fixer = new ShapeFix_Shape(shape);
    fixer->Perform();
    shape = fixer->Shape();

    for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next()) {
        TopoDS_Face face = TopoDS::Face(exp.Current());
        if (face.IsNull()) continue;

        BRepAdaptor_Surface adapt(face, true);
        GeomAbs_SurfaceType st = adapt.GetType();
        if (st == GeomAbs_BSplineSurface || st == GeomAbs_Plane ||
            st == GeomAbs_Cylinder || st == GeomAbs_Cone ||
            st == GeomAbs_Sphere || st == GeomAbs_Torus ||
            st == GeomAbs_BezierSurface || st == GeomAbs_SurfaceOfRevolution ||
            st == GeomAbs_SurfaceOfExtrusion || st == GeomAbs_OffsetSurface ||
            st == GeomAbs_OtherSurface) {
            result.faces.push_back(face);
        }
    }

    if (result.faces.empty()) {
        result.errorMsg = "No faces found in: " + filepath;
        return result;
    }

    result.loaded = true;
    return result;
}

TopoDS_Face findLargestFace(const TopoDS_Shape& shape) {
    TopoDS_Face best;
    double bestArea = -1.0;
    for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next()) {
        TopoDS_Face f = TopoDS::Face(exp.Current());
        if (f.IsNull()) continue;
        BRepAdaptor_Surface adapt(f, true);
        double u1 = adapt.FirstUParameter(), u2 = adapt.LastUParameter();
        double v1 = adapt.FirstVParameter(), v2 = adapt.LastVParameter();
        double area = (u2 - u1) * (v2 - v1);
        if (area > bestArea) {
            bestArea = area;
            best = f;
        }
    }
    return best;
}

void generateFaceMesh(const TopoDS_Face& face, double deflection,
                      Vec3Arr& vertices, std::vector<std::array<int, 3>>& facesOut)
{
    TopoDS_Face f = face;
    f.Location(TopLoc_Location());

    BRepMesh_IncrementalMesh mesher(f, deflection);
    mesher.Perform();

    TopLoc_Location loc;
    Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(f, loc);

    if (tri.IsNull()) return;

    int nVerts = tri->NbNodes();
    vertices.resize(nVerts);
    for (int i = 1; i <= nVerts; ++i) {
        gp_Pnt p = tri->Node(i).Transformed(loc.Transformation());
        vertices[i - 1] = Vec3(p.X(), p.Y(), p.Z());
    }

    int nTris = tri->NbTriangles();
    facesOut.resize(nTris);
    for (int i = 1; i <= nTris; ++i) {
        const Poly_Triangle& t = tri->Triangle(i);
        int n1, n2, n3;
        t.Get(n1, n2, n3);
        facesOut[i - 1] = {n1 - 1, n2 - 1, n3 - 1};
    }
}

} // namespace simple
