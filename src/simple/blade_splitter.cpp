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

namespace simple {

BladeSplitResult splitBladeFaceBySection(
    const TopoDS_Face& face,
    int numSections,
    int samplesPerSection,
    int smoothingWindow,
    double peakThresholdRatio,
    double clusterGapRatio)
{
    BladeSplitResult result;
    std::ostringstream log;

    BRepAdaptor_Surface adapt(face, true);
    if (adapt.GetType() != GeomAbs_BSplineSurface) {
        result.message = "Face is not a BSpline surface";
        return result;
    }

    Handle(Geom_BSplineSurface) surf = adapt.BSpline();
    if (surf.IsNull())
        surf = Handle(Geom_BSplineSurface)::DownCast(BRep_Tool::Surface(face));
    if (surf.IsNull()) {
        result.message = "Cannot obtain BSpline surface";
        return result;
    }

    double uMin = adapt.FirstUParameter(), uMax = adapt.LastUParameter();
    double vMin = adapt.FirstVParameter(), vMax = adapt.LastVParameter();
    double vRange = vMax - vMin;
    double uRange = uMax - uMin;

    if (vRange < 1e-6) {
        result.message = "V range too small";
        return result;
    }

    double vMargin = vRange * 0.08;
    int nSec = std::max(3, numSections);
    std::vector<double> vHeights;
    for (int i = 0; i < nSec; ++i)
        vHeights.push_back(vMin + vMargin + (vRange - 2.0 * vMargin) * i / (nSec - 1));

    struct PeakInfo { double u; double curv; double v; };
    std::vector<PeakInfo> allPeaks;

    for (double vh : vHeights) {
        Handle(Geom_Curve) curve = surf->VIso(vh);
        if (curve.IsNull()) continue;

        GeomAdaptor_Curve gac(curve);
        int ns = std::max(50, samplesPerSection);

        std::vector<std::pair<double, double>> cp;
        cp.reserve(ns);
        for (int i = 0; i < ns; ++i) {
            double t = gac.FirstParameter() + (gac.LastParameter() - gac.FirstParameter()) * i / (ns - 1);
            gp_Pnt p; gp_Vec d1, d2;
            gac.D2(t, p, d1, d2);
            double curv = 0.0;
            double m = d1.Magnitude();
            if (m > 1e-10)
                curv = d1.Crossed(d2).Magnitude() / (m * m * m);
            cp.emplace_back(t, curv);
        }

        int hw = smoothingWindow / 2;
        std::vector<double> smoothed(ns, 0.0);
        for (int i = 0; i < ns; ++i) {
            double sum = 0.0; int cnt = 0;
            for (int j = std::max(0, i - hw); j < std::min(ns, i + hw + 1); ++j) {
                sum += cp[j].second; ++cnt;
            }
            smoothed[i] = (cnt > 0) ? sum / cnt : 0.0;
        }

        double maxCurv = *std::max_element(smoothed.begin(), smoothed.end());
        if (maxCurv < 1e-12) continue;
        double threshold = maxCurv * peakThresholdRatio;

        std::vector<PeakInfo> peaks;
        for (int i = 1; i < ns - 1; ++i) {
            if (smoothed[i] > smoothed[i - 1] && smoothed[i] > smoothed[i + 1]
                && smoothed[i] > threshold) {
                peaks.push_back({cp[i].first, smoothed[i], vh});
            }
        }
        for (auto& p : peaks) allPeaks.push_back(p);
    }

    log << "  Found " << allPeaks.size() << " curvature peaks across "
        << vHeights.size() << " sections";

    if (allPeaks.empty()) {
        result.message = log.str() + " -- no curvature peaks found";
        return result;
    }

    std::vector<double> uVals;
    for (auto& pk : allPeaks) uVals.push_back(pk.u);
    std::sort(uVals.begin(), uVals.end());

    double gapThreshold = uRange * clusterGapRatio;
    std::vector<std::vector<double>> clusters;
    std::vector<double> curCluster;
    curCluster.push_back(uVals[0]);
    for (size_t i = 1; i < uVals.size(); ++i) {
        if (uVals[i] - uVals[i - 1] < gapThreshold)
            curCluster.push_back(uVals[i]);
        else {
            clusters.push_back(curCluster);
            curCluster.clear();
            curCluster.push_back(uVals[i]);
        }
    }
    clusters.push_back(curCluster);

    log << "\n  Clustered into " << clusters.size() << " groups";

    std::vector<double> clusterMeans, clusterCurvs;
    for (auto& cl : clusters) {
        double mean = 0; for (double v : cl) mean += v; mean /= cl.size();
        clusterMeans.push_back(mean);
        double curvSum = 0; int curvCnt = 0;
        for (auto& pk : allPeaks)
            if (std::abs(pk.u - mean) < gapThreshold) { curvSum += pk.curv; ++curvCnt; }
        clusterCurvs.push_back(curvCnt > 0 ? curvSum / curvCnt : 0);
    }

    double maxClusterCurv = *std::max_element(clusterCurvs.begin(), clusterCurvs.end());
    double curvThreshold = maxClusterCurv * 0.25;

    std::vector<int> keptIdx;
    for (int k = 0; k < (int)clusters.size(); ++k)
        if (clusterCurvs[k] > curvThreshold) keptIdx.push_back(k);

    if (keptIdx.empty())
        for (int k = 0; k < (int)clusters.size(); ++k) keptIdx.push_back(k);

    std::sort(keptIdx.begin(), keptIdx.end(), [&](int a, int b) {
        return clusterMeans[a] < clusterMeans[b];
    });

    log << "\n  Kept " << keptIdx.size() << "/" << clusters.size()
        << " clusters above curv " << curvThreshold;

    std::vector<double> boundaries = {uMin};
    for (int k : keptIdx) boundaries.push_back(clusterMeans[k]);
    boundaries.push_back(uMax);
    std::sort(boundaries.begin(), boundaries.end());

    auto computeCurv = [&](double us, double ue) -> double {
        double sum = 0; int cnt = 0;
        for (double v : vHeights) {
            Handle(Geom_Curve) curve = surf->VIso(v);
            if (curve.IsNull()) continue;
            GeomAdaptor_Curve gac(curve);
            int ns2 = std::max(8, (int)((ue - us) / uRange * 30));
            for (int i = 0; i < ns2; ++i) {
                double u = us + (ue - us) * (i + 1) / (ns2 + 1.0);
                gp_Pnt p; gp_Vec d1, d2;
                gac.D2(u, p, d1, d2);
                double cv = 0, m = d1.Magnitude();
                if (m > 1e-10) cv = d1.Crossed(d2).Magnitude() / (m * m * m);
                sum += cv; ++cnt;
            }
        }
        return cnt > 0 ? sum / cnt : 0;
    };

    std::vector<SplitRegion> regions;
    for (size_t i = 0; i + 1 < boundaries.size(); ++i) {
        SplitRegion r;
        r.uStart = boundaries[i];
        r.uEnd = boundaries[i + 1];
        r.avgCurv = computeCurv(r.uStart, r.uEnd);
        r.label = "";
        regions.push_back(r);
    }

    std::vector<int> idx(regions.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) {
        return regions[a].avgCurv > regions[b].avgCurv;
    });

    double maxC = regions[idx[0]].avgCurv;
    double minC = regions[idx.back()].avgCurv;

    for (int i = 0; i < (int)idx.size(); ++i) {
        int ri = idx[i];
        double len = regions[ri].uEnd - regions[ri].uStart;
        if (i == 0 && regions[ri].avgCurv > maxC * 0.5)
            regions[ri].label = "edge";
        else if (regions[ri].avgCurv >= (maxC + minC) * 0.5)
            regions[ri].label = "suction";
        else
            regions[ri].label = "pressure";
    }

    for (auto& r : regions)
        if (r.label.empty())
            r.label = "pressure";

    result.regions = regions;
    result.success = true;

    for (auto& r : regions) {
        log << "\n  [" << r.uStart << "," << r.uEnd << "] len="
            << (r.uEnd - r.uStart) << " curv=" << r.avgCurv
            << " -> " << r.label;
    }
    result.message = log.str();
    return result;
}

} // namespace simple
