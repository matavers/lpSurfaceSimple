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
    if (surf.IsNull()) {
        surf = Handle(Geom_BSplineSurface)::DownCast(BRep_Tool::Surface(face));
    }
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
    for (int i = 0; i < nSec; ++i) {
        double v = vMin + vMargin + (vRange - 2.0 * vMargin) * i / (nSec - 1);
        vHeights.push_back(v);
    }

    struct PeakInfo {
        double u;
        double curv;
        double v;
    };
    std::vector<PeakInfo> allPeaks;

    for (double vh : vHeights) {
        Handle(Geom_Curve) curve = surf->VIso(vh);
        if (curve.IsNull()) continue;

        GeomAdaptor_Curve gac(curve);
        double cMin = gac.FirstParameter(), cMax = gac.LastParameter();
        int ns = std::max(50, samplesPerSection);

        std::vector<std::pair<double, double>> cp; // (u, curvature)
        cp.reserve(ns);

        for (int i = 0; i < ns; ++i) {
            double t = cMin + (cMax - cMin) * i / (ns - 1);
            gp_Pnt p;
            gp_Vec d1, d2;
            gac.D2(t, p, d1, d2);

            double curv = 0.0;
            double d1Mag = d1.Magnitude();
            if (d1Mag > 1e-10) {
                gp_Vec crossVec = d1.Crossed(d2);
                curv = crossVec.Magnitude() / (d1Mag * d1Mag * d1Mag);
            }
            cp.emplace_back(t, curv);
        }

        int hw = smoothingWindow / 2;
        std::vector<double> smoothed(ns, 0.0);
        for (int i = 0; i < ns; ++i) {
            double sum = 0.0;
            int cnt = 0;
            for (int j = std::max(0, i - hw); j < std::min(ns, i + hw + 1); ++j) {
                sum += cp[j].second;
                ++cnt;
            }
            smoothed[i] = (cnt > 0) ? sum / cnt : 0.0;
        }

        double maxCurv = *std::max_element(smoothed.begin(), smoothed.end());
        if (maxCurv < 1e-12) continue;
        double threshold = maxCurv * 0.15;

        std::vector<PeakInfo> peaks;
        for (int i = 1; i < ns - 1; ++i) {
            if (smoothed[i] > smoothed[i - 1] && smoothed[i] > smoothed[i + 1]
                && smoothed[i] > threshold)
            {
                PeakInfo pk;
                pk.u = cp[i].first;
                pk.curv = smoothed[i];
                pk.v = vh;
                peaks.push_back(pk);
            }
        }

        std::sort(peaks.begin(), peaks.end(),
            [](const PeakInfo& a, const PeakInfo& b) { return a.curv > b.curv; });

        int take = std::min(3, (int)peaks.size());
        for (int k = 0; k < take; ++k)
            allPeaks.push_back(peaks[k]);
    }

    log << "\n  Found " << allPeaks.size() << " curvature peaks across "
        << vHeights.size() << " sections";

    if (allPeaks.empty()) {
        result.message = log.str() + " — no curvature peaks found";
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
        if (uVals[i] - uVals[i - 1] < gapThreshold) {
            curCluster.push_back(uVals[i]);
        } else {
            clusters.push_back(curCluster);
            curCluster.clear();
            curCluster.push_back(uVals[i]);
        }
    }
    clusters.push_back(curCluster);

    log << "\n  Clustered into " << clusters.size() << " groups";

    if (clusters.size() < 2) {
        if (clusters.size() == 1) {
            double medU = clusters[0][clusters[0].size() / 2];
            double u1 = medU;
            double u2 = std::fmod(medU + uRange * 0.5 + uMin, uRange) + uMin;
            log << "\n  Only 1 cluster — heuristic split at u=" << u1 << " and u=" << u2;
            result.uLE = u1;
            result.uTE = u2;
        } else {
            result.message = log.str() + " — need at least 2 curvature peak clusters";
            return result;
        }
    } else {
        std::vector<double> clusterMeans;
        std::vector<double> clusterCurvs;
        for (auto& cl : clusters) {
            double mean = 0;
            for (double v : cl) mean += v;
            mean /= cl.size();
            clusterMeans.push_back(mean);

            double curvSum = 0;
            int curvCnt = 0;
            for (auto& pk : allPeaks) {
                if (std::abs(pk.u - mean) < gapThreshold) {
                    curvSum += pk.curv;
                    ++curvCnt;
                }
            }
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

        result.uLE = clusterMeans[0];
        result.uTE = clusterMeans[1];
        result.leCurv = clusterCurvs[0];
        result.teCurv = clusterCurvs[1];
    }

    double u1 = result.uLE;
    double u2 = result.uTE;
    if (u1 > u2) std::swap(u1, u2);

    double arcLenDirect = u2 - u1;
    double arcLenWrap = uRange - arcLenDirect;

    std::vector<double> curvDirect, curvWrap;
    for (double vh : vHeights) {
        Handle(Geom_Curve) curve = surf->VIso(vh);
        if (curve.IsNull()) continue;
        GeomAdaptor_Curve gac(curve);

        for (int i = 0; i < 30; ++i) {
            double t = u1 + arcLenDirect * (i + 1) / 31.0;
            gp_Pnt p; gp_Vec d1, d2;
            gac.D2(t, p, d1, d2);
            double cv = 0;
            double m = d1.Magnitude();
            if (m > 1e-10) cv = d1.Crossed(d2).Magnitude() / (m * m * m);
            curvDirect.push_back(cv);
        }

        int nw = std::max(1, (int)(arcLenWrap / uRange * 30));
        for (int i = 0; i < nw; ++i) {
            double t = u2 + arcLenWrap * (i + 1) / (nw + 1.0);
            if (t > uMax) t -= uRange;
            gp_Pnt p; gp_Vec d1, d2;
            gac.D2(t, p, d1, d2);
            double cv = 0;
            double m = d1.Magnitude();
            if (m > 1e-10) cv = d1.Crossed(d2).Magnitude() / (m * m * m);
            curvWrap.push_back(cv);
        }
    }

    double avgDirect = curvDirect.empty() ? 0 :
        std::accumulate(curvDirect.begin(), curvDirect.end(), 0.0) / curvDirect.size();
    double avgWrap = curvWrap.empty() ? 0 :
        std::accumulate(curvWrap.begin(), curvWrap.end(), 0.0) / curvWrap.size();

    if (avgWrap > avgDirect) {
        result.uSuctStart = u2;
        result.uSuctEnd = u1;
        result.uPressStart = u1;
        result.uPressEnd = u2;
    } else {
        result.uPressStart = u2;
        result.uPressEnd = u1;
        result.uSuctStart = u1;
        result.uSuctEnd = u2;
    }

    result.success = true;
    log << "\n  LE u=" << result.uLE << " curv=" << result.leCurv;
    log << "\n  TE u=" << result.uTE << " curv=" << result.teCurv;
    log << "\n  Direct arc U=[" << u1 << "," << u2 << "] len=" << arcLenDirect
        << " avgCurv=" << avgDirect;
    log << "\n  Wrap   arc U=[" << u2 << "," << u1 << "] len=" << arcLenWrap
        << " avgCurv=" << avgWrap;
    log << "\n  Pressure U=[" << result.uPressStart << "," << result.uPressEnd << "]";
    log << "\n  Suction  U=[" << result.uSuctStart << "," << result.uSuctEnd << "]";
    result.message = log.str();

    return result;
}

} // namespace simple
