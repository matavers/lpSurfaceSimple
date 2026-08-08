#include "simple/blade_identifier.hpp"

#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <BRepLProp_SLProps.hxx>
#include <BRep_Tool.hxx>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>

namespace simple {

struct CandidateInfo {
    int originalIndex;
    double diag3D;
    double paramArea;
    gp_Vec avgNormal;
    double avgMeanCurv;
};

static double sampleMeanCurvature(BRepAdaptor_Surface& adapt,
                                   double u, double v,
                                   double tol = 1e-6)
{
    BRepLProp_SLProps props(adapt, u, v, 2, tol);
    if (props.IsCurvatureDefined())
        return props.MeanCurvature();
    return 0.0;
}

BladeIdentifyResult identifyBladeSurfaces(
    const std::vector<TopoDS_Face>& faces,
    int normalSampleGrid,
    double minParamArea,
    double minDiag3D,
    int curvatureSampleGrid)
{
    BladeIdentifyResult result;
    std::ostringstream log;

    if (faces.empty()) {
        result.message = "No faces provided";
        return result;
    }

    struct RawInfo {
        int index;
        double diag;
        double pArea;
        TopoDS_Face face;
    };
    std::vector<RawInfo> raw;

    for (size_t i = 0; i < faces.size(); ++i) {
        BRepAdaptor_Surface adapt(faces[i], true);
        if (adapt.GetType() != GeomAbs_BSplineSurface)
            continue;

        double uRange = adapt.LastUParameter() - adapt.FirstUParameter();
        double vRange = adapt.LastVParameter() - adapt.FirstVParameter();
        double pArea = uRange * vRange;
        if (pArea < minParamArea)
            continue;

        Bnd_Box bbox;
        BRepBndLib::Add(faces[i], bbox);
        double x1, y1, z1, x2, y2, z2;
        bbox.Get(x1, y1, z1, x2, y2, z2);
        double diag = std::sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
        if (diag < minDiag3D)
            continue;

        raw.push_back({static_cast<int>(i), diag, pArea, faces[i]});
    }

    if (raw.size() < 2) {
        log << "Filtered " << raw.size() << " NURBS candidate(s) from "
            << faces.size() << " face(s) — need at least 2";
        result.message = log.str();
        return result;
    }

    std::vector<double> diags;
    for (auto& r : raw) diags.push_back(r.diag);
    std::sort(diags.begin(), diags.end());
    double medDiag = diags[diags.size() / 2];
    double maxDiag3D = medDiag * 5.0;

    log << "Filtered " << raw.size() << " NURBS face(s) (minArea=" << minParamArea
        << " minDiag=" << minDiag3D << " maxDiag=" << maxDiag3D << ")";

    std::vector<CandidateInfo> candidates;

    for (auto& r : raw) {
        if (r.diag > maxDiag3D) {
            log << "\n  skip face[" << r.index << "] oversized diag=" << r.diag;
            continue;
        }

        BRepAdaptor_Surface adapt(r.face, true);
        double uMin = adapt.FirstUParameter(), uMax = adapt.LastUParameter();
        double vMin = adapt.FirstVParameter(), vMax = adapt.LastVParameter();
        double uRange = uMax - uMin;
        double vRange = vMax - vMin;

        gp_Vec sumN(0, 0, 0);
        int nDefined = 0;

        for (int ui = 0; ui < normalSampleGrid; ++ui) {
            for (int vi = 0; vi < normalSampleGrid; ++vi) {
                double u = uMin + uRange * ui / std::max(1, normalSampleGrid - 1);
                double v = vMin + vRange * vi / std::max(1, normalSampleGrid - 1);
                BRepLProp_SLProps props(adapt, u, v, 1, 1e-6);
                if (props.IsNormalDefined()) {
                    sumN += props.Normal();
                    ++nDefined;
                }
            }
        }
        if (nDefined == 0)
            continue;
        sumN.Normalize();

        double curvSum = 0.0;
        int curvDef = 0;
        for (int ui = 0; ui < curvatureSampleGrid; ++ui) {
            for (int vi = 0; vi < curvatureSampleGrid; ++vi) {
                double u = uMin + uRange * ui / std::max(1, curvatureSampleGrid - 1);
                double v = vMin + vRange * vi / std::max(1, curvatureSampleGrid - 1);
                BRepLProp_SLProps props(adapt, u, v, 2, 1e-4);
                if (props.IsCurvatureDefined()) {
                    curvSum += props.MeanCurvature();
                    ++curvDef;
                }
            }
        }
        double curv = (curvDef > 0) ? curvSum / curvDef : 0.0;

        if (std::abs(curv) < 1e-4) {
            log << "\n  skip face[" << r.index << "] near-flat (curv=" << curv << ")";
            continue;
        }

        CandidateInfo ci;
        ci.originalIndex = r.index;
        ci.diag3D = r.diag;
        ci.paramArea = r.pArea;
        ci.avgNormal = sumN;
        ci.avgMeanCurv = curv;
        candidates.push_back(ci);

        log << "\n    face[" << r.index << "] diag=" << r.diag
            << " curv=" << curv;
    }

    if (candidates.size() < 2) {
        if (candidates.size() == 1) {
            log << "\n  (single blade candidate — needs splitting)";
            result.pressureFaceIndex = candidates[0].originalIndex;
            result.suctionFaceIndex = candidates[0].originalIndex;
            result.success = true;
            result.pressureLabel = "blade";
            result.suctionLabel = "blade";
            result.message = log.str();
            return result;
        }
        log << "\n  need at least 1 candidate";
        result.message = log.str();
        return result;
    }

    std::vector<int> group0, group1;
    group0.push_back(0);

    for (size_t i = 1; i < candidates.size(); ++i) {
        double dot0 = std::abs(candidates[i].avgNormal.Dot(candidates[group0.front()].avgNormal));
        double dot1 = 0.0;
        if (!group1.empty())
            dot1 = std::abs(candidates[i].avgNormal.Dot(candidates[group1.front()].avgNormal));

        if (group1.empty()) {
            if (dot0 < 0.5)
                group1.push_back(static_cast<int>(i));
            else
                group0.push_back(static_cast<int>(i));
        } else {
            if (dot0 >= dot1)
                group0.push_back(static_cast<int>(i));
            else
                group1.push_back(static_cast<int>(i));
        }
    }

    if (group1.empty()) {
        if (candidates.size() >= 2) {
            auto pickLargestTwo = [&]() -> std::pair<int,int> {
                std::vector<int> idx(candidates.size());
                for (size_t i = 0; i < candidates.size(); ++i) idx[i] = static_cast<int>(i);
                std::sort(idx.begin(), idx.end(), [&](int a, int b) {
                    return candidates[a].diag3D > candidates[b].diag3D;
                });
                return {idx[0], idx[1]};
            };
            auto [a, b] = pickLargestTwo();
            log << "\n  (all faces share similar normals — returning 2 largest)";
            result.pressureFaceIndex = candidates[a].originalIndex;
            result.suctionFaceIndex = candidates[b].originalIndex;
            result.success = true;
            result.pressureLabel = "faceA";
            result.suctionLabel = "faceB";
            result.message = log.str();
            return result;
        }
        log << "\n  could not cluster into 2 groups";
        result.message = log.str();
        return result;
    }

    auto computeGroupCurv = [&](const std::vector<int>& grp) -> double {
        double sum = 0.0;
        for (int idx : grp)
            sum += candidates[idx].avgMeanCurv;
        return sum / grp.size();
    };

    double curv0 = computeGroupCurv(group0);
    double curv1 = computeGroupCurv(group1);

    const std::vector<int>* pressGrp = nullptr;
    const std::vector<int>* suctGrp = nullptr;

    if (std::abs(curv0) >= std::abs(curv1)) {
        suctGrp = &group0;
        pressGrp = &group1;
    } else {
        suctGrp = &group1;
        pressGrp = &group0;
    }

    auto pickLargest = [&](const std::vector<int>& grp) -> int {
        int best = grp.front();
        for (int idx : grp)
            if (candidates[idx].diag3D > candidates[best].diag3D)
                best = idx;
        return best;
    };

    int pi = pickLargest(*pressGrp);
    int si = pickLargest(*suctGrp);

    result.pressureFaceIndex = candidates[pi].originalIndex;
    result.suctionFaceIndex = candidates[si].originalIndex;
    result.success = true;
    result.pressureLabel = "pressure";
    result.suctionLabel = "suction";

    log << "\n  Group A (" << group0.size() << " faces): avgCurv=" << curv0;
    log << "\n  Group B (" << group1.size() << " faces): avgCurv=" << curv1;
    log << "\n  Pressure face[" << result.pressureFaceIndex
        << "] diag=" << candidates[pi].diag3D;
    log << "\n  Suction  face[" << result.suctionFaceIndex
        << "] diag=" << candidates[si].diag3D;
    result.message = log.str();

    return result;
}

} // namespace simple
