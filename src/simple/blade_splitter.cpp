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

    // —— 叶缘识别：封闭面环绕，找曲率峰（前缘/后缘）与谷（叶盆/叶背）——
    // 把最深谷旋转到环绕点，保证两个叶缘峰都在内部、峰谷交替、不跨环绕。
    int sh = 0;
    for (int i = 1; i < nPts; ++i) if (sm[i] < sm[sh]) sh = i;
    auto s = [&](int i) { return sm[(i + sh) % nPts]; };
    auto vOf = [&](int i) { return vMin + vRng * (double)i / (nPts - 1.0); };

    std::vector<int> peaks, valleys;
    for (int i = 1; i < nPts - 1; ++i) {
        if (s(i) >= s(i - 1) && s(i) >= s(i + 1) && s(i) > mn + (mx - mn) * 0.05)
            peaks.push_back(i);
        if (s(i) <= s(i - 1) && s(i) <= s(i + 1))
            valleys.push_back(i);
    }
    if (peaks.size() > 2) {
        std::sort(peaks.begin(), peaks.end(), [&](int a, int b) { return s(a) > s(b); });
        peaks.resize(2);
    }
    if (valleys.size() > 2) {
        std::sort(valleys.begin(), valleys.end(), [&](int a, int b) { return s(a) < s(b); });
        valleys.resize(2);
    }

    // 合并峰谷，排序，相邻中点作为切分点（含环绕）
    std::vector<int> ext;
    for (int p : peaks) ext.push_back(p);
    for (int v : valleys) ext.push_back(v);
    std::sort(ext.begin(), ext.end());

    std::vector<double> bp;
    if (ext.size() >= 2) {
        for (size_t k = 0; k < ext.size(); ++k) {
            int a = ext[k], b = ext[(k + 1) % ext.size()];
            int d = (b - a + nPts) % nPts;
            int mid = (a + d / 2) % nPts;
            bp.push_back(vOf(mid));
        }
        std::sort(bp.begin(), bp.end());
        std::vector<double> bpu;
        for (double c : bp)
            if (bpu.empty() || c - bpu.back() > vRng * 0.005) bpu.push_back(c);
        bp = bpu;
    }

    auto rcWrap = [&](double a, double b) {
        int ia = (int)((a - vMin) / vRng * (nPts - 1.0));
        int ib = (int)((b - vMin) / vRng * (nPts - 1.0));
        double s = 0; int c = 0;
        for (int i = ia; ; ++i) {
            s += sm[i % nPts]; ++c;
            if (i % nPts == ib % nPts) break;
        }
        return c > 0 ? s / c : 0.0;
    };

    // 构建区域（环绕）：相邻切分点之间为一段
    std::vector<SplitRegion> regs;
    int nCuts = (int)bp.size();
    if (nCuts >= 2) {
        for (int k = 0; k < nCuts; ++k) {
            double a = bp[k];
            double b = bp[(k + 1) % nCuts];
            double end = (k == nCuts - 1) ? b + vRng : b;  // 最后一段跨环绕
            SplitRegion r;
            r.uStart = a;
            r.uEnd = end;
            r.avgCurv = rcWrap(a, end);
            regs.push_back(r);
        }
    } else {
        SplitRegion r; r.uStart = vMin; r.uEnd = vMax; r.avgCurv = rcWrap(vMin, vMax);
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
