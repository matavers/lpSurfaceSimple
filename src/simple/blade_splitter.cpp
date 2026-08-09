#include "simple/blade_splitter.hpp"

#include <BRepAdaptor_Surface.hxx>
#include <BRep_Tool.hxx>
#include <BRepLProp_SLProps.hxx>
#include <Geom_BSplineSurface.hxx>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <numeric>

namespace simple {

struct MCurvProfile {
    std::vector<double> sm;
    double dMin, dMax, dRange;
};

static MCurvProfile buildMCProfile(
    BRepAdaptor_Surface& adapt,
    double uMin, double uMax, double vMin, double vMax,
    int nSec, int nPts, int hw, bool alongU)
{
    MCurvProfile cp;
    cp.dMin = alongU ? uMin : vMin;
    cp.dMax = alongU ? uMax : vMax;
    cp.dRange = cp.dMax - cp.dMin;

    double crossMin = alongU ? vMin : uMin;
    double crossMax = alongU ? vMax : uMax;
    double crossRng = crossMax - crossMin;
    double crossMargin = crossRng * 0.08;
    nSec = std::max(3, nSec); nPts = std::max(50, nPts);

    cp.sm.resize(nPts, 0.0);
    std::vector<int> cnt(nPts, 0);

    for (int si = 0; si < nSec; ++si) {
        double cross = crossMin + crossMargin + (crossRng - 2.0 * crossMargin) * si / (nSec - 1);
        for (int i = 0; i < nPts; ++i) {
            double t = cp.dMin + cp.dRange * i / (nPts - 1);
            double u = alongU ? t : cross;
            double v = alongU ? cross : t;
            BRepLProp_SLProps props(adapt, u, v, 2, 1e-4);
            if (props.IsCurvatureDefined()) {
                cp.sm[i] += std::abs(props.MeanCurvature());
                cnt[i]++;
            }
        }
    }
    for (int i = 0; i < nPts; ++i)
        if (cnt[i] > 0) cp.sm[i] /= cnt[i];

    std::vector<double> tmp = cp.sm;
    hw = 2;
    for (int i = 0; i < nPts; ++i) {
        double s = 0; int c = 0;
        for (int j = std::max(0, i - hw); j < std::min(nPts, i + hw + 1); ++j)
            { s += tmp[j]; ++c; }
        cp.sm[i] = c > 0 ? s / c : 0;
    }
    return cp;
}

BladeSplitResult splitBladeFaceBySection(
    const TopoDS_Face& face, int nSec, int nPts, int hw,
    double peakThreshRatio, double)
{
    BladeSplitResult result;
    std::ostringstream log;

    BRepAdaptor_Surface adapt(face, true);
    Handle(Geom_BSplineSurface) surf = adapt.BSpline();
    if (surf.IsNull()) surf = Handle(Geom_BSplineSurface)::DownCast(BRep_Tool::Surface(face));
    if (!surf.IsNull()) {} // unused
    if (adapt.GetType() != GeomAbs_BSplineSurface) {
        result.message = "Not a BSpline surface"; return result;
    }

    double uMin = adapt.FirstUParameter(), uMax = adapt.LastUParameter();
    double vMin = adapt.FirstVParameter(), vMax = adapt.LastVParameter();

    auto cpU = buildMCProfile(adapt, uMin, uMax, vMin, vMax, nSec, nPts, hw, true);
    auto cpV = buildMCProfile(adapt, uMin, uMax, vMin, vMax, nSec, nPts, hw, false);

    auto countPeaks = [](const MCurvProfile& cp) {
        double mx = *std::max_element(cp.sm.begin(), cp.sm.end());
        double mn = *std::min_element(cp.sm.begin(), cp.sm.end());
        double th = mn + (mx - mn) * 0.1;
        int c = 0, n = (int)cp.sm.size();
        for (int i = 1; i < n - 1; ++i)
            if (cp.sm[i] > cp.sm[i-1] && cp.sm[i] > cp.sm[i+1] && cp.sm[i] > th) ++c;
        return c;
    };

    int pU = countPeaks(cpU), pV = countPeaks(cpV);

    double rngU = *std::max_element(cpU.sm.begin(), cpU.sm.end())
                - *std::min_element(cpU.sm.begin(), cpU.sm.end());
    double rngV = *std::max_element(cpV.sm.begin(), cpV.sm.end())
                - *std::min_element(cpV.sm.begin(), cpV.sm.end());

    bool alongU = rngU >= rngV && pU >= 2;
    if (pV >= 2 && rngV * 2 > rngU) alongU = false;

    auto& cp = alongU ? cpU : cpV;

    log << "  Direction:" << (alongU ? "U" : "V")
        << " rngU=" << rngU << " rngV=" << rngV
        << " peaks U=" << pU << " V=" << pV;
    log << " range [" << cp.dMin << "," << cp.dMax << "]";

    double mxC = *std::max_element(cp.sm.begin(), cp.sm.end());
    double mnC = *std::min_element(cp.sm.begin(), cp.sm.end());
    double rngC = mxC - mnC;
    if (rngC < 1e-10) { result.message = log.str(); return result; }

    int N = (int)cp.sm.size();
    double th = mnC + rngC * 0.05;

    std::vector<std::pair<double, double>> peaks;
    for (int i = 1; i < N - 1; ++i) {
        if (cp.sm[i] > cp.sm[i-1] && cp.sm[i] > cp.sm[i+1] && cp.sm[i] > th) {
            double u = cp.dMin + cp.dRange * i / (N - 1);
            if (peaks.empty() || (u - peaks.back().first) > cp.dRange * 0.02)
                peaks.push_back({u, cp.sm[i]});
        }
    }

    bool boundaryLow0 = cp.sm[0] > cp.sm[1] && cp.sm[0] > cp.sm[2] && cp.sm[0] > th;
    bool boundaryHigh1 = cp.sm[N-1] > cp.sm[N-2] && cp.sm[N-1] > cp.sm[N-3] && cp.sm[N-1] > th;

    double edgeCurvTh = mxC * 0.3;
    int edgeStartU = -1, edgeEndU = -1;

    if (boundaryLow0 && cp.sm[0] > edgeCurvTh) {
        int cut = 0;
        for (int i = 1; i < N; ++i)
            if (cp.sm[i] < edgeCurvTh) { cut = i; break; }
        if (cut > 0) {
            double eu = cp.dMin + cp.dRange * cut / (N - 1);
            peaks.insert(peaks.begin(), {eu, cp.sm[0]});
        }
    }
    if (boundaryHigh1 && cp.sm[N-1] > edgeCurvTh) {
        int cut = N - 1;
        for (int i = N - 2; i >= 0; --i)
            if (cp.sm[i] < edgeCurvTh) { cut = i; break; }
        if (cut < N - 1) {
            double eu = cp.dMin + cp.dRange * cut / (N - 1);
            peaks.push_back({eu, cp.sm[N-1]});
        }
    }

    log << "\n  " << peaks.size() << " raw peaks at th=" << th;
    for (auto& p : peaks) log << " u=" << p.first;

    if (peaks.size() > 3) {
        std::sort(peaks.begin(), peaks.end(),
            [](auto& a, auto& b) { return a.second > b.second; });
        double mxP = peaks[0].second;
        double kTh = mxP * 0.5;
        std::vector<std::pair<double, double>> strong;
        for (auto& p : peaks) if (p.second >= kTh) strong.push_back(p);
        if (strong.size() >= 2) peaks = strong;
    }

    std::sort(peaks.begin(), peaks.end(),
        [](auto& a, auto& b) { return a.first < b.first; });

    log << "\n  " << peaks.size() << " significant peaks";

    std::vector<double> bp = {cp.dMin};
    for (auto& p : peaks) bp.push_back(p.first);
    bp.push_back(cp.dMax);

    auto regC = [&](double a, double b) {
        int is = std::max(0, (int)((a - cp.dMin) / cp.dRange * (N - 1)));
        int ie = std::min(N - 1, (int)((b - cp.dMin) / cp.dRange * (N - 1)));
        double s = 0; int c = 0;
        for (int i = is; i <= ie; ++i) { s += cp.sm[i]; ++c; }
        return c > 0 ? s / c : 0;
    };

    std::vector<SplitRegion> regs;
    for (size_t i = 0; i + 1 < bp.size(); ++i) {
        SplitRegion r; r.uStart = bp[i]; r.uEnd = bp[i+1];
        r.avgCurv = regC(r.uStart, r.uEnd);
        regs.push_back(r);
    }

    double regMax = 0, regMin = 1e30;
    bool hasNonEdge = false;
    for (auto& r : regs) {
        double len = r.uEnd - r.uStart;
        if (len < cp.dRange * 1e-6) continue;
        if (r.avgCurv > mxC * 0.3) continue;
        if (r.avgCurv > regMax) regMax = r.avgCurv;
        if (r.avgCurv < regMin) regMin = r.avgCurv;
        hasNonEdge = true;
    }
    if (!hasNonEdge) { regMax = regs[0].avgCurv; regMin = regs[0].avgCurv; }
    double regMed = (regMax + regMin) * 0.5;
    double edgeTh = (regMax + regMed) * 0.5;

    for (auto& r : regs) {
        double len = r.uEnd - r.uStart;
        if (len < cp.dRange * 1e-6) continue;
        if (r.avgCurv > mxC * 0.3)
            r.label = "edge";
        else if (r.avgCurv >= regMed)
            r.label = "suction";
        else
            r.label = "pressure";
    }

    std::vector<SplitRegion> filtered;
    for (auto& r : regs)
        if (r.uEnd - r.uStart > cp.dRange * 1e-6)
            filtered.push_back(r);
    regs = filtered;

    result.regions = regs;
    result.success = true;
    result.splitDir = alongU ? 'U' : 'V';
    for (auto& r : regs)
        log << "\n  [" << r.uStart << "," << r.uEnd << "] len="
            << (r.uEnd - r.uStart) << " C=" << r.avgCurv << " " << r.label;
    result.message = log.str();
    return result;
}

} // namespace simple
