#include "simple/blade_splitter.hpp"

#include <BRepAdaptor_Surface.hxx>
#include <BRep_Tool.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Curve.hxx>
#include <GeomAdaptor_Curve.hxx>

#include <algorithm>
#include <cmath>
#include <sstream>

namespace simple {

BladeSplitResult splitBladeFaceBySection(
    const TopoDS_Face& face, int nSec, int nPts, int hw, double, double)
{
    BladeSplitResult result;
    std::ostringstream log;

    BRepAdaptor_Surface adapt(face, true);
    Handle(Geom_BSplineSurface) surf = adapt.BSpline();
    if (surf.IsNull()) surf = Handle(Geom_BSplineSurface)::DownCast(BRep_Tool::Surface(face));
    if (surf.IsNull()) { result.message = "No BSpline"; return result; }

    double uMin = adapt.FirstUParameter(), uMax = adapt.LastUParameter();
    double vMin = adapt.FirstVParameter(), vMax = adapt.LastVParameter();
    double vRng = vMax - vMin, uRng = uMax - uMin;
    if (vRng < 1e-6) { result.message = "V too small"; return result; }

    nSec = std::max(3, nSec); nPts = std::max(50, nPts);
    double vm = vRng * 0.08;
    std::vector<double> vH;
    for (int i = 0; i < nSec; ++i) vH.push_back(vMin + vm + (vRng - 2*vm) * i / (nSec-1));

    std::vector<double> cv(nPts, 0.0);
    std::vector<int> ct(nPts, 0);
    for (double vh : vH) {
        Handle(Geom_Curve) c = surf->VIso(vh);
        if (c.IsNull()) continue;
        GeomAdaptor_Curve gac(c);
        double t0 = gac.FirstParameter(), t1 = gac.LastParameter();
        for (int i = 0; i < nPts; ++i) {
            double t = t0 + (t1 - t0) * i / (nPts - 1);
            gp_Pnt p; gp_Vec d1, d2;
            gac.D2(t, p, d1, d2);
            double m = d1.Magnitude();
            cv[i] += (m > 1e-10) ? d1.Crossed(d2).Magnitude() / (m*m*m) : 0;
            ct[i]++;
        }
    }
    for (int i = 0; i < nPts; ++i) if (ct[i] > 0) cv[i] /= ct[i];

    std::vector<double> sm(nPts);
    int h = hw / 2;
    for (int i = 0; i < nPts; ++i) {
        double s = 0; int c = 0;
        for (int j = std::max(0,i-h); j < std::min(nPts,i+h+1); ++j) { s += cv[j]; ++c; }
        sm[i] = c > 0 ? s / c : 0;
    }

    double mx = *std::max_element(sm.begin(), sm.end());
    double mn = *std::min_element(sm.begin(), sm.end());

    std::vector<double> bp = {uMin};
    for (int i = 2; i < nPts - 2; ++i) {
        if (sm[i] > sm[i-1] && sm[i] > sm[i+1] && sm[i] > mn + (mx-mn)*0.08) {
            double u = uMin + uRng * i / (nPts - 1);
            if (bp.empty() || u - bp.back() > uRng * 0.015) bp.push_back(u);
        }
    }
    bp.push_back(uMax);

    auto rc = [&](double a, double b) {
        int is = std::max(0, (int)((a-uMin)/uRng*(nPts-1)));
        int ie = std::min(nPts-1, (int)((b-uMin)/uRng*(nPts-1)));
        double s = 0; for (int i=is;i<=ie;++i)s+=sm[i];
        return s/std::max(1,ie-is+1);
    };

    std::vector<SplitRegion> regs;
    for (size_t i = 0; i+1 < bp.size(); ++i) {
        SplitRegion r; r.uStart = bp[i]; r.uEnd = bp[i+1]; r.avgCurv = rc(r.uStart, r.uEnd);
        regs.push_back(r);
    }

    double rMax = 0, rMin = 1e30;
    for (auto& r : regs) { if(r.avgCurv>rMax)rMax=r.avgCurv; if(r.avgCurv<rMin)rMin=r.avgCurv; }
    double rMed = (rMax + rMin) * 0.5;

    for (auto& r : regs) {
        double L = r.uEnd - r.uStart;
        if (L < uRng * 0.10) r.label = "edge";
        else if (r.avgCurv >= rMed) r.label = "suction";
        else r.label = "pressure";
    }

    std::vector<SplitRegion> mr;
    for (auto& r : regs) {
        if (!mr.empty() && mr.back().label == r.label) {
            mr.back().uEnd = r.uEnd;
            mr.back().avgCurv = (mr.back().avgCurv + r.avgCurv) * 0.5;
        } else mr.push_back(r);
    }

    for (auto& r : mr) if (r.label.empty()) r.label = "pressure";

    result.regions = mr;
    result.success = true;
    result.splitDir = 'U';
    log << "  maxC=" << mx << "  rawReg=" << regs.size() << " merged=" << mr.size();
    for (auto& r : mr)
        log << "\n  [" << r.uStart << "," << r.uEnd << "] L="
            << (r.uEnd-r.uStart) << " C=" << r.avgCurv << " " << r.label;
    result.message = log.str();
    return result;
}

} // namespace simple
