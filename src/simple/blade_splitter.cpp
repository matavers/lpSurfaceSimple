#include "simple/blade_splitter.hpp"

#include <BRepAdaptor_Surface.hxx>
#include <BRep_Tool.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Curve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <BRepLProp_SLProps.hxx>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <numeric>

namespace simple {

BladeSplitResult splitBladeFaceBySection(
    const TopoDS_Face& face, int numSections, int samplesPerSection,
    int smoothingWindow, double peakThresholdRatio, double clusterGapRatio)
{
    BladeSplitResult result;
    std::ostringstream log;

    BRepAdaptor_Surface adapt(face, true);
    Handle(Geom_BSplineSurface) surf = adapt.BSpline();
    if (surf.IsNull()) surf = Handle(Geom_BSplineSurface)::DownCast(BRep_Tool::Surface(face));
    if (surf.IsNull()) { result.message = "No BSpline surface"; return result; }

    double uMin = adapt.FirstUParameter(), uMax = adapt.LastUParameter();
    double vMin = adapt.FirstVParameter(), vMax = adapt.LastVParameter();
    double uRange = uMax - uMin, vRange = vMax - vMin;

    if (vRange < 1e-6) { result.message = "V range too small"; return result; }

    int nSec = std::max(3, numSections);
    double vMargin = vRange * 0.08;
    std::vector<double> vHeights;
    for (int i = 0; i < nSec; ++i)
        vHeights.push_back(vMin + vMargin + (vRange - 2.0 * vMargin) * i / (nSec - 1));

    struct Peak { double u; double curv; };
    std::vector<Peak> allPeaks;
    int ns = std::max(50, samplesPerSection);

    for (double vh : vHeights) {
        Handle(Geom_Curve) curve = surf->VIso(vh);
        if (curve.IsNull()) continue;
        GeomAdaptor_Curve gac(curve);

        std::vector<double> curvProfile(ns);
        for (int i = 0; i < ns; ++i) {
            double t = gac.FirstParameter() + (gac.LastParameter() - gac.FirstParameter()) * i / (ns - 1);
            gp_Pnt p; gp_Vec d1, d2;
            gac.D2(t, p, d1, d2);
            double m = d1.Magnitude();
            curvProfile[i] = (m > 1e-10) ? d1.Crossed(d2).Magnitude() / (m * m * m) : 0;
        }

        int hw = smoothingWindow / 2;
        std::vector<double> sm(ns);
        for (int i = 0; i < ns; ++i) {
            double s = 0; int c = 0;
            for (int j = std::max(0, i - hw); j < std::min(ns, i + hw + 1); ++j) { s += curvProfile[j]; ++c; }
            sm[i] = c > 0 ? s / c : 0;
        }

        double maxC = *std::max_element(sm.begin(), sm.end());
        if (maxC < 1e-12) continue;

        for (int i = 1; i < ns - 1; ++i)
            if (sm[i] > sm[i - 1] && sm[i] > sm[i + 1] && sm[i] > maxC * peakThresholdRatio)
                allPeaks.push_back({gac.FirstParameter() + (gac.LastParameter() - gac.FirstParameter()) * i / (ns - 1), sm[i]});
    }

    log << "  Found " << allPeaks.size() << " peaks across " << vHeights.size() << " sections";
    if (allPeaks.empty()) { result.message = log.str(); return result; }

    std::vector<double> uVals;
    for (auto& p : allPeaks) uVals.push_back(p.u);
    std::sort(uVals.begin(), uVals.end());

    double gap = uRange * clusterGapRatio;
    std::vector<std::vector<double>> clusters;
    std::vector<double> cur; cur.push_back(uVals[0]);
    for (size_t i = 1; i < uVals.size(); ++i) {
        if (uVals[i] - uVals[i - 1] < gap) cur.push_back(uVals[i]);
        else { clusters.push_back(cur); cur.clear(); cur.push_back(uVals[i]); }
    }
    clusters.push_back(cur);

    std::vector<double> means, curvs;
    for (auto& cl : clusters) {
        double m = 0; for (double v : cl) m += v; m /= cl.size(); means.push_back(m);
        double cs = 0; int cc = 0;
        for (auto& p : allPeaks) if (std::abs(p.u - m) < gap) { cs += p.curv; ++cc; }
        curvs.push_back(cc > 0 ? cs / cc : 0);
    }

    std::vector<int> idx(clusters.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return curvs[a] > curvs[b]; });

    double maxCC = curvs[idx[0]], minCC = curvs[idx.back()];
    double curvThreshold = (maxCC + minCC) * 0.6;

    std::vector<double> keeMeans, keeCurvs;
    for (int k = 0; k < (int)clusters.size(); ++k) {
        int ci = idx[k];
        if (curvs[ci] >= curvThreshold) {
            keeMeans.push_back(means[ci]);
            keeCurvs.push_back(curvs[ci]);
        }
    }
    if (keeMeans.empty()) { keeMeans.push_back(means[idx[0]]); keeCurvs.push_back(curvs[idx[0]]); }

    log << "\n  Keeping " << keeMeans.size() << "/" << clusters.size()
        << " clusters (curv >= " << curvThreshold << " mid=" << (maxCC + minCC) / 2 << ")";

    std::vector<int> order(keeMeans.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return keeMeans[a] < keeMeans[b]; });

    std::vector<double> bp = {uMin};
    for (int o : order) bp.push_back(keeMeans[o]);
    bp.push_back(uMax);

    auto regionCurv = [&](double us, double ue) -> double {
        double sum = 0; int cnt = 0;
        for (double vh : vHeights) {
            Handle(Geom_Curve) curve = surf->VIso(vh);
            if (curve.IsNull()) continue;
            GeomAdaptor_Curve gac(curve);
            int nr = std::max(10, (int)((ue - us) / uRange * 30));
            for (int i = 0; i < nr; ++i) {
                double u = us + (ue - us) * (i + 1) / (nr + 1.0);
                gp_Pnt p; gp_Vec d1, d2;
                gac.D2(u, p, d1, d2);
                double m = d1.Magnitude();
                if (m > 1e-10) sum += d1.Crossed(d2).Magnitude() / (m * m * m);
                ++cnt;
            }
        }
        return cnt > 0 ? sum / cnt : 0;
    };

    std::vector<SplitRegion> regions;
    for (size_t i = 0; i + 1 < bp.size(); ++i) {
        SplitRegion r;
        r.uStart = bp[i]; r.uEnd = bp[i + 1];
        r.avgCurv = regionCurv(r.uStart, r.uEnd);
        regions.push_back(r);
    }

    int n = (int)regions.size();
    std::vector<int> ri(n);
    std::iota(ri.begin(), ri.end(), 0);
    std::sort(ri.begin(), ri.end(), [&](int a, int b) { return regions[a].avgCurv > regions[b].avgCurv; });

    for (int i = 0; i < n; ++i) {
        int r = ri[i];
        double len = regions[r].uEnd - regions[r].uStart;
        if (i < 2 && len < uRange * 0.15)
            regions[r].label = "edge";
        else if (regions[r].avgCurv >= regions[ri[n / 2]].avgCurv)
            regions[r].label = "suction";
        else
            regions[r].label = "pressure";
    }

    for (auto& r : regions)
        if (r.label.empty()) r.label = "pressure";

    result.regions = regions;
    result.success = true;
    for (auto& r : regions)
        log << "\n  [" << r.uStart << "," << r.uEnd << "] len=" << (r.uEnd - r.uStart)
            << " C=" << r.avgCurv << " " << r.label;
    result.message = log.str();
    return result;
}

} // namespace simple
