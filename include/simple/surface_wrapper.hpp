#pragma once

#include "common.hpp"

#include <array>

#include <Geom_BSplineSurface.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Curve.hxx>
#include <gp_Pnt.hxx>

namespace simple {

using FaceArr = std::vector<std::array<int, 3>>;

class SurfaceWrapper {
public:
    SurfaceWrapper() = default;
    SurfaceWrapper(const Handle(Geom_BSplineSurface)& surf,
                   double uMinOverride = -1.0, double uMaxOverride = -1.0,
                   double vMinOverride = -1.0, double vMaxOverride = -1.0,
                   bool wrapU = false);
    ~SurfaceWrapper() = default;

    const Handle(Geom_BSplineSurface)& surface() const { return m_surf; }

    Vec3 evaluate(double u, double v) const;
    Vec3 evaluateWrapped(double t, double v) const;
    Vec3 normal(double u, double v) const;

    void generateMesh(int resU, int resV,
                      Vec3Arr& vertices, FaceArr& faces,
                      Vec2Arr& uvs) const;

    Handle(Geom_BSplineCurve) extractIsoCurveU(double u) const;
    Handle(Geom_BSplineCurve) extractIsoCurveV(double v) const;

    std::pair<double, double> paramDomainU() const;
    std::pair<double, double> paramDomainV() const;
    bool isWrapU() const { return m_wrapU; }

    int numCtrlU() const { return m_surf->NbUPoles(); }
    int numCtrlV() const { return m_surf->NbVPoles(); }
    int degreeU() const { return m_surf->UDegree(); }
    int degreeV() const { return m_surf->VDegree(); }

private:
    Handle(Geom_BSplineSurface) m_surf;
    double m_uMinOverride = -1.0;
    double m_uMaxOverride = -1.0;
    double m_vMinOverride = -1.0;
    double m_vMaxOverride = -1.0;
    bool m_wrapU = false;
};

Vec3Arr sampleCurve(const Handle(Geom_BSplineCurve)& curve, int nSamples);
Vec3Arr sampleCurveRange(const Handle(Geom_BSplineCurve)& curve,
                         double tStart, double tEnd, int nSamples);
Vec3Arr sampleCurveRangeWrap(const Handle(Geom_BSplineCurve)& curve,
                              double uStart, double uEnd, int nSamples);

} // namespace simple
