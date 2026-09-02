#pragma once

#include "common.hpp"
#include "surface_wrapper.hpp"

#include <array>

#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>

namespace simple {

// 单个格子（参数域矩形）的直纹面拟合结果。
// fitDir 为全局一致的直纹面方向（U 或 V），由网格拟合在每次二分后统一试算确定。
struct RuledCellFit {
    ParamDir fitDir;
    Handle(Geom_BSplineCurve) curveC0;
    Handle(Geom_BSplineCurve) curveC1;
    Vec3Arr curveC0Samples;
    Vec3Arr curveC1Samples;
    Vec3Arr ruledMeshVerts;
    FaceArr ruledMeshFaces;
    double maxError;
    double rmsError;
    double twist;        // 母线法向扭转角最大值(度)，0≈可展
};

// 单格直纹面拟合（指定方向），内部含准线优化（最小二乘 + lambda 正则）。
RuledCellFit fitCellRuled(const SurfaceWrapper& surf,
                          double u0, double u1, double v0, double v1,
                          ParamDir dir,
                          int nUSamples, int nVSamples,
                          int nRibs, double lambda);

Vec3Arr generateRuledMesh(const Vec3Arr& c0Samples,
                           const Vec3Arr& c1Samples,
                           int nAlong, int nAcross,
                           FaceArr& faces);

void optimizeDirectrices(
    const SurfaceWrapper& surf,
    double uSeg0, double uSeg1,
    double vSeg0, double vSeg1,
    ParamDir directrixDir,
    Vec3Arr& c0Samples,
    Vec3Arr& c1Samples,
    int nRibs, double lambda);

bool exportOBJ(const std::string& path,
               const Vec3Arr& verts, const FaceArr& faces);

// 导出优化后的直纹面准线（采样点形式），供 CAM 使用。
bool exportDirectrixTXT(const std::string& path,
                        const Vec3Arr& c0Samples, const Vec3Arr& c1Samples);

} // namespace simple
