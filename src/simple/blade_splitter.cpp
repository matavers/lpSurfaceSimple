#include "simple/blade_splitter.hpp"

#include <BRepAdaptor_Surface.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>

#include <Geom_BSplineSurface.hxx>
#include <Geom_Curve.hxx>
#include <GeomAdaptor_Curve.hxx>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <numeric>
#include <map>

namespace simple {

BladeSplitResult splitBladeFaceBySection(
    const TopoDS_Face& face,
    int numSections,
    int samplesPerSection,
    int smoothingWindow)
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
    double uRange = uMax - uMin;
    double vRange = vMax - vMin;

    bool isPeriodic = surf->IsUPeriodic();
    gp_Pnt p0 = adapt.Value(uMin, vMin + vRange * 0.5);
    gp_Pnt p1 = adapt.Value(uMax, vMin + vRange * 0.5);
    double closureDist = p0.Distance(p1);

    log << "Surface U=[" << uMin << "," << uMax << "] V=[" << vMin << "," << vMax << "]"
        << " periodic=" << (isPeriodic ? "yes" : "no")
        << " closureDist=" << closureDist;

    if (vRange < 1e-6) {
        result.message = "V range too small for sections";
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
        double cMin = gac.FirstParameter(), cMax = gac.LastParameter();
        int ns = std::max(50, samplesPerSection);

        std::vector<std::pair<double, double>> cp;
        cp.reserve(ns);
        for (int i = 0; i < ns; ++i) {
            double t = cMin + (cMax - cMin) * i / (ns - 1);
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
        double threshold = maxCurv * 0.15;

        std::vector<PeakInfo> peaks;
        for (int i = 1; i < ns - 1; ++i) {
            if (smoothed[i] > smoothed[i - 1] && smoothed[i] > smoothed[i + 1]
                && smoothed[i] > threshold) {
                peaks.push_back({cp[i].first, smoothed[i], vh});
            }
        }
        std::sort(peaks.begin(), peaks.end(),
            [](const PeakInfo& a, const PeakInfo& b) { return a.curv > b.curv; });
        int take = std::min(3, (int)peaks.size());
        for (int k = 0; k < take; ++k) allPeaks.push_back(peaks[k]);
    }

    log << "\n  Found " << allPeaks.size() << " curvature peaks across "
        << vHeights.size() << " sections";

    if (allPeaks.empty()) {
        result.message = log.str() + " -- no curvature peaks found";
        return result;
    }

    std::vector<double> uVals;
    for (auto& pk : allPeaks) uVals.push_back(pk.u);
    std::sort(uVals.begin(), uVals.end());

    double gapThreshold = uRange * 0.15;
    std::vector<std::vector<double>> clusters;
    std::vector<double> curCluster;
    curCluster.push_back(uVals[0]);
    for (size_t i = 1; i < uVals.size(); ++i) {
        if (uVals[i] - uVals[i - 1] < gapThreshold)
            curCluster.push_back(uVals[i]);
        else { clusters.push_back(curCluster); curCluster.clear(); curCluster.push_back(uVals[i]); }
    }
    clusters.push_back(curCluster);

    log << "\n  Clustered into " << clusters.size() << " groups";

    if (clusters.size() < 2) {
        result.message = log.str() + " -- need at least 2 curvature peak clusters";
        return result;
    }

    std::vector<double> clusterMeans, clusterCurvs;
    for (auto& cl : clusters) {
        double mean = 0; for (double v : cl) mean += v; mean /= cl.size();
        clusterMeans.push_back(mean);
        double curvSum = 0; int curvCnt = 0;
        for (auto& pk : allPeaks)
            if (std::abs(pk.u - mean) < gapThreshold) { curvSum += pk.curv; ++curvCnt; }
        clusterCurvs.push_back(curvCnt > 0 ? curvSum / curvCnt : 0);
    }

    if (clusters.size() > 2) {
        std::vector<int> idx(clusters.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(),
            [&](int a, int b) { return clusterCurvs[a] > clusterCurvs[b]; });
        clusters = {clusters[idx[0]], clusters[idx[1]]};
        clusterMeans = {clusterMeans[idx[0]], clusterMeans[idx[1]]};
        clusterCurvs = {clusterCurvs[idx[0]], clusterCurvs[idx[1]]};
        log << "\n  Keeping top 2 clusters by curvature";
    }

    result.uPeak1 = clusterMeans[0];
    result.uPeak2 = clusterMeans[1];
    result.curvPeak1 = clusterCurvs[0];
    result.curvPeak2 = clusterCurvs[1];
    if (result.uPeak1 > result.uPeak2) {
        std::swap(result.uPeak1, result.uPeak2);
        std::swap(result.curvPeak1, result.curvPeak2);
    }

    double up1 = result.uPeak1, up2 = result.uPeak2;
    double len1 = up1;
    double len2 = up2 - up1;
    double len3 = uMax - up2;

    log << "\n  Peak1 u=" << up1 << " curv=" << result.curvPeak1;
    log << "\n  Peak2 u=" << up2 << " curv=" << result.curvPeak2;

    auto computeRegionCurv = [&](double us, double ue, const std::vector<double>& vh) -> double {
        double sum = 0; int cnt = 0;
        for (double v : vh) {
            Handle(Geom_Curve) curve = surf->VIso(v);
            if (curve.IsNull()) continue;
            GeomAdaptor_Curve gac(curve);
            int ns = std::max(10, (int)((ue - us) / uRange * 30));
            for (int i = 0; i < ns; ++i) {
                double u = us + (ue - us) * (i + 1) / (ns + 1.0);
                gp_Pnt p; gp_Vec d1, d2;
                gac.D2(u, p, d1, d2);
                double cv = 0, m = d1.Magnitude();
                if (m > 1e-10) cv = d1.Crossed(d2).Magnitude() / (m * m * m);
                sum += cv; ++cnt;
            }
        }
        return cnt > 0 ? sum / cnt : 0;
    };

    double curvReg1 = computeRegionCurv(0, up1, vHeights);
    double curvReg2 = computeRegionCurv(up1, up2, vHeights);
    double curvReg3 = computeRegionCurv(up2, uMax, vHeights);

    log << "\n  Region [0," << up1 << "] len=" << len1 << " avgCurv=" << curvReg1;
    log << "\n  Region [" << up1 << "," << up2 << "] len=" << len2 << " avgCurv=" << curvReg2;
    log << "\n  Region [" << up2 << "," << uMax << "] len=" << len3 << " avgCurv=" << curvReg3;

    int edgeIdx = 0;
    if (curvReg1 >= curvReg2 && curvReg1 >= curvReg3) edgeIdx = 0;
    else if (curvReg2 >= curvReg3) edgeIdx = 1;
    else edgeIdx = 2;

    int sideA = -1, sideB = -1;
    for (int k = 0; k < 3; ++k) {
        if (k == edgeIdx) continue;
        if (sideA < 0) sideA = k; else sideB = k;
    }

    double curvSegs[3] = {curvReg1, curvReg2, curvReg3};
    double lens[3] = {len1, len2, len3};

    if (curvSegs[sideA] < 0 && curvSegs[sideB] > 0) {}
    else if (curvSegs[sideB] < 0 && curvSegs[sideA] > 0) { std::swap(sideA, sideB); }

    double pressU0, pressU1, suctU0, suctU1;
    if (curvSegs[sideA] >= curvSegs[sideB]) {
        suctU0 = (sideA == 0) ? 0.0 : (sideA == 1 ? up1 : up2);
        suctU1 = (sideA == 0) ? up1 : (sideA == 1 ? up2 : uMax);
        pressU0 = (sideB == 0) ? 0.0 : (sideB == 1 ? up1 : up2);
        pressU1 = (sideB == 0) ? up1 : (sideB == 1 ? up2 : uMax);
    } else {
        pressU0 = (sideA == 0) ? 0.0 : (sideA == 1 ? up1 : up2);
        pressU1 = (sideA == 0) ? up1 : (sideA == 1 ? up2 : uMax);
        suctU0 = (sideB == 0) ? 0.0 : (sideB == 1 ? up1 : up2);
        suctU1 = (sideB == 0) ? up1 : (sideB == 1 ? up2 : uMax);
    }
    std::string pressLabel = "pressure", suctLabel = "suction";

    if (pressU0 > pressU1) {
        std::swap(pressU0, pressU1);
        if (pressU0 < up1 && pressU1 > up2) { pressLabel = "suction"; suctLabel = "pressure"; }
    }
    if (suctU0 > suctU1) {
        std::swap(suctU0, suctU1);
        if (suctU0 < up1 && suctU1 > up2) { suctLabel = "pressure"; pressLabel = "suction"; }
    }

    result.uPressStart = pressU0;
    result.uPressEnd = pressU1;
    result.uSuctStart = suctU0;
    result.uSuctEnd = suctU1;
    result.curvPress = curvSegs[sideA];
    result.curvSuct = curvSegs[sideB];

    result.uEdgeStart = (edgeIdx == 0) ? 0.0 : (edgeIdx == 1 ? up1 : up2);
    result.uEdgeEnd = (edgeIdx == 0) ? up1 : (edgeIdx == 1 ? up2 : uMax);

    result.success = true;
    log << "\n  Edge    U=[" << result.uEdgeStart << "," << result.uEdgeEnd << "]"
        << " (" << lens[edgeIdx] << ")";
    log << "\n  " << pressLabel << " U=[" << result.uPressStart
        << "," << result.uPressEnd << "] curv=" << result.curvPress;
    log << "\n  " << suctLabel << " U=[" << result.uSuctStart
        << "," << result.uSuctEnd << "] curv=" << result.curvSuct;
    result.message = log.str();

    return result;
}

} // namespace simple
