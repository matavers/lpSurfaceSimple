#include "simple/blade_splitter.hpp"

#include <BRepAdaptor_Surface.hxx>
#include <BRep_Tool.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Curve.hxx>
#include <GeomAdaptor_Curve.hxx>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <numeric>
#include <set>

namespace simple {

static double curveCurv(const GeomAdaptor_Curve& gac, double t) {
    gp_Pnt p; gp_Vec d1, d2;
    gac.D2(t, p, d1, d2);
    double m = d1.Magnitude();
    return m > 1e-10 ? d1.Crossed(d2).Magnitude() / (m * m * m) : 0.0;
}

struct CProfile {
    std::vector<double> sm;
    double dMin, dMax, dRange;
};

static CProfile buildProfile(Handle(Geom_BSplineSurface)& surf,
    double uMin, double uMax, double vMin, double vMax,
    int nSec, int nPts, int hw, bool alongU)
{
    CProfile cp;
    cp.dMin = alongU ? uMin : vMin;
    cp.dMax = alongU ? uMax : vMax;
    cp.dRange = cp.dMax - cp.dMin;

    double v0 = alongU ? vMin : uMin, v1 = alongU ? vMax : uMax;
    double vr = v1 - v0, vm = vr * 0.08;
    nSec = std::max(3, nSec); nPts = std::max(50, nPts);

    cp.sm.resize(nPts, 0.0);
    std::vector<int> cnt(nPts, 0);

    for (int si = 0; si < nSec; ++si) {
        double vh = v0 + vm + (vr - 2 * vm) * si / (nSec - 1);
        Handle(Geom_Curve) crv = alongU ? surf->VIso(vh) : surf->UIso(vh);
        if (crv.IsNull()) continue;
        GeomAdaptor_Curve gac(crv);
        double t0 = gac.FirstParameter(), t1 = gac.LastParameter();
        for (int i = 0; i < nPts; ++i) {
            double t = t0 + (t1 - t0) * i / (nPts - 1);
            cp.sm[i] += curveCurv(gac, t);
            cnt[i]++;
        }
    }
    for (int i = 0; i < nPts; ++i)
        if (cnt[i] > 0) cp.sm[i] /= cnt[i];

    std::vector<double> tmp = cp.sm;
    hw /= 2;
    for (int i = 0; i < nPts; ++i) {
        double s = 0; int c = 0;
        for (int j = std::max(0, i - hw); j < std::min(nPts, i + hw + 1); ++j)
            { s += tmp[j]; ++c; }
        cp.sm[i] = c > 0 ? s / c : 0;
    }
    return cp;
}

BladeSplitResult splitBladeFaceBySection(
    const TopoDS_Face& face, int nSec, int nPts, int hw, double, double)
{
    BladeSplitResult result;
    std::ostringstream log;

    BRepAdaptor_Surface adapt(face, true);
    Handle(Geom_BSplineSurface) surf = adapt.BSpline();
    if (surf.IsNull()) surf = Handle(Geom_BSplineSurface)::DownCast(BRep_Tool::Surface(face));
    if (surf.IsNull()) { result.message = "No BSpline surface"; return result; }

    double uMin = adapt.FirstUParameter(), uMax = adapt.LastUParameter();
    double vMin = adapt.FirstVParameter(), vMax = adapt.LastVParameter();

    auto cpU = buildProfile(surf, uMin, uMax, vMin, vMax, nSec, nPts, hw, true);
    auto cpV = buildProfile(surf, uMin, uMax, vMin, vMax, nSec, nPts, hw, false);

    auto featCount = [](const CProfile& cp) {
        double mx = *std::max_element(cp.sm.begin(), cp.sm.end());
        double mn = *std::min_element(cp.sm.begin(), cp.sm.end());
        double th = (mx - mn) * 0.1;
        int c = 0, n = (int)cp.sm.size();
        for (int i = 2; i < n - 2; ++i)
            if (cp.sm[i] > cp.sm[i-1] && cp.sm[i] > cp.sm[i+1] && cp.sm[i] > mn + th) ++c;
        return c;
    };

    int fU = featCount(cpU), fV = featCount(cpV);
    bool alongU = fU >= fV;
    auto& cp = alongU ? cpU : cpV;

    log << "  Direction: " << (alongU ? "U" : "V")
        << " (peaks U=" << fU << " V=" << fV << ")";
    log << " range [" << cp.dMin << "," << cp.dMax << "]";

    double mxC = *std::max_element(cp.sm.begin(), cp.sm.end());
    double mnC = *std::min_element(cp.sm.begin(), cp.sm.end());
    double rng = mxC - mnC;
    if (rng < 1e-12) { result.message = log.str(); return result; }

    int N = (int)cp.sm.size();
    double th = mnC + rng * 0.2;

    std::vector<std::pair<double, double>> peaksWithCurv;
    for (int i = 2; i < N - 2; ++i) {
        if (cp.sm[i] > cp.sm[i-1] && cp.sm[i] > cp.sm[i+1] &&
            cp.sm[i] > cp.sm[i-2] && cp.sm[i] > cp.sm[i+2] && cp.sm[i] > th) {
            double u = cp.dMin + cp.dRange * i / (N - 1);
            if (peaksWithCurv.empty() || (u - peaksWithCurv.back().first) > cp.dRange * 0.04)
                peaksWithCurv.push_back({u, cp.sm[i]});
        }
    }

    log << "\n  " << peaksWithCurv.size() << " raw peaks";

    if (peaksWithCurv.size() > 3) {
        std::sort(peaksWithCurv.begin(), peaksWithCurv.end(),
            [](auto& a, auto& b) { return a.second > b.second; });
        double maxC = peaksWithCurv[0].second;
        double keepTh = maxC * 0.65;
        std::vector<std::pair<double, double>> strong;
        for (auto& p : peaksWithCurv)
            if (p.second >= keepTh) strong.push_back(p);
        if (strong.size() >= 2) peaksWithCurv = strong;
        else peaksWithCurv.resize(std::min(2, (int)peaksWithCurv.size()));
    }
    std::sort(peaksWithCurv.begin(), peaksWithCurv.end(),
        [](auto& a, auto& b) { return a.first < b.first; });

    log << "\n  " << peaksWithCurv.size() << " significant peaks";
    for (auto& p : peaksWithCurv) log << " u=" << p.first;

    std::vector<double> bp = {cp.dMin};
    for (auto& p : peaksWithCurv) bp.push_back(p.first);
    bp.push_back(cp.dMax);

    auto regCurv = [&](double a, double b) {
        int is = std::max(0, (int)((a - cp.dMin) / cp.dRange * (N - 1)));
        int ie = std::min(N - 1, (int)((b - cp.dMin) / cp.dRange * (N - 1)));
        double s = 0; int c = 0;
        for (int i = is; i <= ie; ++i) { s += cp.sm[i]; ++c; }
        return c > 0 ? s / c : 0;
    };

    std::vector<SplitRegion> regs;
    for (size_t i = 0; i + 1 < bp.size(); ++i) {
        SplitRegion r; r.uStart = bp[i]; r.uEnd = bp[i+1];
        r.avgCurv = regCurv(r.uStart, r.uEnd);
        regs.push_back(r);
    }

    int n = (int)regs.size();
    double regMax = regs[0].avgCurv, regMin = regs[0].avgCurv;
    for (auto& r : regs) {
        if (r.avgCurv > regMax) regMax = r.avgCurv;
        if (r.avgCurv < regMin) regMin = r.avgCurv;
    }
    double regMed = (regMax + regMin) * 0.5;
    for (auto& r : regs) {
        double len = r.uEnd - r.uStart;
        if (len < cp.dRange * 0.15 && r.avgCurv > regMed)
            r.label = "edge";
        else if (r.avgCurv >= regMed)
            r.label = "suction";
        else
            r.label = "pressure";
    }

    result.regions = regs;
    result.success = true;
    for (auto& r : regs)
        log << "\n  [" << r.uStart << "," << r.uEnd << "] len="
            << (r.uEnd - r.uStart) << " C=" << r.avgCurv << " " << r.label;
    result.message = log.str();
    return result;
}

} // namespace simple
