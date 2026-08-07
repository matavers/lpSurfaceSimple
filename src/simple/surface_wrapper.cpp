#include "simple/surface_wrapper.hpp"

#include <GeomLib.hxx>
#include <gp_Vec.hxx>
#include <gp_Dir.hxx>

namespace simple {

SurfaceWrapper::SurfaceWrapper(const Handle(Geom_BSplineSurface)& surf,
                               double uMinOverride, double uMaxOverride,
                               double vMinOverride, double vMaxOverride,
                               bool wrapU)
    : m_surf(surf), m_uMinOverride(uMinOverride), m_uMaxOverride(uMaxOverride),
      m_vMinOverride(vMinOverride), m_vMaxOverride(vMaxOverride), m_wrapU(wrapU) {}

Vec3 SurfaceWrapper::evaluate(double u, double v) const {
    gp_Pnt p;
    m_surf->D0(u, v, p);
    return Vec3(p.X(), p.Y(), p.Z());
}

Vec3 SurfaceWrapper::evaluateWrapped(double t, double v) const {
    auto [umin, umax] = paramDomainU();
    double uRange = (m_wrapU && umax < umin) ? (1.0 - umin + umax) : (umax - umin);
    double u = umin + uRange * t;
    if (m_wrapU && u > 1.0) u -= 1.0;
    return evaluate(u, v);
}

Vec3 SurfaceWrapper::normal(double u, double v) const {
    gp_Pnt S;
    gp_Vec Su, Sv;
    m_surf->D1(u, v, S, Su, Sv);
    gp_Vec n = Su.Crossed(Sv);
    double mag = n.Magnitude();
    if (mag < 1e-10) return Vec3(0, 0, 1);
    return Vec3(n.X() / mag, n.Y() / mag, n.Z() / mag);
}

void SurfaceWrapper::generateMesh(int resU, int resV,
                                   Vec3Arr& vertices,
                                   FaceArr& faces,
                                   Vec2Arr& uvs) const
{
    auto [umin, umax] = paramDomainU();
    auto [vmin, vmax] = paramDomainV();
    int nU = resU + 1, nV = resV + 1;

    vertices.resize(nU * nV);
    uvs.resize(nU * nV);

    double uRange = (m_wrapU && umax < umin) ? (1.0 - umin + umax) : (umax - umin);

    for (int i = 0; i < nU; ++i) {
        double t = (double)i / resU;
        double u = umin + uRange * t;
        if (m_wrapU && u > 1.0) u -= 1.0;
        for (int j = 0; j < nV; ++j) {
            double v = vmin + (vmax - vmin) * j / resV;
            int idx = i * nV + j;
            vertices[idx] = evaluate(u, v);
            uvs[idx] = Vec2(t, v);
        }
    }

    faces.clear();
    for (int i = 0; i < nU - 1; ++i) {
        for (int j = 0; j < nV - 1; ++j) {
            int a = i * nV + j;
            faces.push_back({a, a + 1, a + nV});
            faces.push_back({a + 1, a + nV + 1, a + nV});
        }
    }
}

Handle(Geom_BSplineCurve) SurfaceWrapper::extractIsoCurveU(double u) const {
    Handle(Geom_Curve) curve = m_surf->UIso(u);
    return Handle(Geom_BSplineCurve)::DownCast(curve);
}

Handle(Geom_BSplineCurve) SurfaceWrapper::extractIsoCurveV(double v) const {
    Handle(Geom_Curve) curve = m_surf->VIso(v);
    return Handle(Geom_BSplineCurve)::DownCast(curve);
}

std::pair<double, double> SurfaceWrapper::paramDomainU() const {
    if (m_uMinOverride >= 0.0) {
        return {m_uMinOverride, m_uMaxOverride};
    }
    return {m_surf->UKnot(m_surf->FirstUKnotIndex()),
            m_surf->UKnot(m_surf->LastUKnotIndex())};
}

std::pair<double, double> SurfaceWrapper::paramDomainV() const {
    if (m_vMinOverride >= 0.0) {
        return {m_vMinOverride, m_vMaxOverride};
    }
    return {m_surf->VKnot(m_surf->FirstVKnotIndex()),
            m_surf->VKnot(m_surf->LastVKnotIndex())};
}

Vec3Arr sampleCurve(const Handle(Geom_BSplineCurve)& curve, int nSamples) {
    Vec3Arr pts(nSamples);
    double t0 = curve->FirstParameter();
    double t1 = curve->LastParameter();
    for (int i = 0; i < nSamples; ++i) {
        double t = t0 + (t1 - t0) * i / (nSamples - 1);
        gp_Pnt p;
        curve->D0(t, p);
        pts[i] = Vec3(p.X(), p.Y(), p.Z());
    }
    return pts;
}

Vec3Arr sampleCurveRange(const Handle(Geom_BSplineCurve)& curve,
                         double tStart, double tEnd, int nSamples) {
    Vec3Arr pts(nSamples);
    for (int i = 0; i < nSamples; ++i) {
        double t = tStart + (tEnd - tStart) * i / (nSamples - 1);
        gp_Pnt p;
        curve->D0(t, p);
        pts[i] = Vec3(p.X(), p.Y(), p.Z());
    }
    return pts;
}

Vec3Arr sampleCurveRangeWrap(const Handle(Geom_BSplineCurve)& curve,
                              double uStart, double uEnd, int nSamples) {
    if (uStart <= uEnd || nSamples <= 1)
        return sampleCurveRange(curve, uStart, uEnd, nSamples);

    double totalRange = (1.0 - uStart) + uEnd;
    int n1 = std::max(1, (int)std::round(nSamples * (1.0 - uStart) / totalRange));
    int n2 = nSamples - n1;

    Vec3Arr pts1 = sampleCurveRange(curve, uStart, 1.0, n1);
    Vec3Arr pts2 = sampleCurveRange(curve, 0.0, uEnd, n2);

    Vec3Arr result;
    result.reserve(nSamples);
    result.insert(result.end(), pts1.begin(), pts1.end() - (n2 > 0 ? 0 : 0));
    if (!pts1.empty() && !pts2.empty()) {
        result.insert(result.end(), pts2.begin() + 1, pts2.end());
    } else if (pts2.empty()) {
        result = pts1;
    } else {
        result = pts2;
    }
    return result;
}

} // namespace simple
