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
    double uRng = uMax - uMin, vRng = vMax - vMin;
    if (uRng < 1e-6) { result.message = "U too small"; return result; }

    nSec = std::max(3, nSec); nPts = std::max(50, nPts);
    double um = uRng * 0.08;
    std::vector<double> uH;
    for (int i = 0; i < nSec; ++i)
        uH.push_back(uMin + um + (uRng - 2*um) * i / (nSec-1));

    std::vector<double> cv(nPts, 0.0);
    std::vector<int> ct(nPts, 0);
    for (double uh : uH) {
        Handle(Geom_Curve) c = surf->UIso(uh);
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
    int h = 1;
    for (int i = 0; i < nPts; ++i) {
        double s = 0; int c = 0;
        for (int j = std::max(0,i-h); j < std::min(nPts,i+h+1); ++j) { s += cv[j]; ++c; }
        sm[i] = c > 0 ? s / c : 0;
    }

    double mx = *std::max_element(sm.begin(), sm.end());
    double mn = *std::min_element(sm.begin(), sm.end());
    log << " maxC=" << mx;

    std::vector<double> bp = {vMin};
    for (int i = 2; i < nPts - 2; ++i) {
        if (sm[i] > sm[i-1] && sm[i] > sm[i+1] &&
            sm[i] > sm[i-2] && sm[i] > sm[i+2] && sm[i] > mn + (mx-mn)*0.05) {
            double v = vMin + vRng * i / (nPts - 1);
            if (bp.empty() || v - bp.back() > vRng * 0.005) bp.push_back(v);
        }
    }
    bp.push_back(vMax);

    // 前缘/后缘分别检测：边界处曲率相对自身下降即切出边缘条带（不再共用全局阈值，
    // 避免一边叶缘曲率远小于另一边时被漏掉）。
    auto cutEdge = [&](bool fromStart) {
        double edgeC = fromStart ? sm[0] : sm[nPts - 1];
        if (edgeC < 1e-12) return -1.0;
        double th = edgeC * 0.25;
        if (fromStart) {
            for (int i = 1; i < nPts; ++i)
                if (sm[i] < th) return vMin + vRng * i / (nPts - 1.0);
        } else {
            for (int i = nPts - 2; i >= 0; --i)
                if (sm[i] < th) return vMin + vRng * i / (nPts - 1.0);
        }
        return -1.0;
    };
    double e0 = cutEdge(true);
    if (e0 > vMin && e0 < vMax) bp.insert(bp.begin() + 1, e0);
    double e1 = cutEdge(false);
    if (e1 > vMin && e1 < vMax) bp.insert(bp.end() - 1, e1);

    auto rc = [&](double a, double b) {
        int is = std::max(0, (int)((a-vMin)/vRng*(nPts-1)));
        int ie = std::min(nPts-1, (int)((b-vMin)/vRng*(nPts-1)));
        double s = 0;
        for (int i=is;i<=ie;++i) s+=sm[i];
        return s / std::max(1, ie-is+1);
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
        if (L < vRng * 0.08) r.label = "edge";
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

    if (mr.size() == 3 && mr[1].label == "pressure" && mr[1].uEnd - mr[1].uStart > vRng * 0.6) {
        SplitRegion side = mr[1];
        mr.erase(mr.begin() + 1);
        double mid = (side.uStart + side.uEnd) * 0.5;
        SplitRegion r1, r2;
        r1.uStart = side.uStart; r1.uEnd = mid;
        r2.uStart = mid; r2.uEnd = side.uEnd;
        int is = std::max(0, (int)((side.uStart-vMin)/vRng*(nPts-1)));
        int im = std::max(0, (int)((mid-vMin)/vRng*(nPts-1)));
        int ie = std::min(nPts-1, (int)((side.uEnd-vMin)/vRng*(nPts-1)));
        double s1=0,s2=0;int c1=0,c2=0;
        for(int i=is;i<im;i++){s1+=sm[i];c1++;}
        for(int i=im;i<=ie;i++){s2+=sm[i];c2++;}
        r1.avgCurv=c1>0?s1/c1:0; r2.avgCurv=c2>0?s2/c2:0;
        r1.label=(r1.avgCurv>=r2.avgCurv)?"suction":"pressure";
        r2.label=(r1.label=="suction")?"pressure":"suction";
        mr.insert(mr.begin() + 1, r1);
        mr.insert(mr.begin() + 2, r2);
    }

    result.regions = mr;
    result.success = true;
    result.splitDir = 'V';
    for (auto& r : mr)
        log << "\n  V[" << r.uStart << "," << r.uEnd << "] L="
            << (r.uEnd-r.uStart) << " C=" << r.avgCurv << " " << r.label;
    result.message = log.str();
    return result;
}

} // namespace simple
