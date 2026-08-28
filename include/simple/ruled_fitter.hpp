#pragma once

#include "common.hpp"
#include "surface_wrapper.hpp"

#include <array>

#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>

namespace simple {

// 单个格子（参数域矩形）的直纹面拟合结果。
// fitDir 为直纹面方向：V = 准线沿 u、母线沿 v；U = 准线沿 v、母线沿 u。
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
};

// 等距分段直纹面拟合结果（主程序 / 兼容旧接口使用）。
struct RuledSegment {
    int segmentIndex;
    ParamDir directrixDir;
    Handle(Geom_BSplineCurve) curveC0;
    Handle(Geom_BSplineCurve) curveC1;
    Vec3Arr curveC0Samples;
    Vec3Arr curveC1Samples;
    Vec3Arr ruledMeshVerts;
    FaceArr ruledMeshFaces;
    double maxError;
    double rmsError;
};

struct RuledResult {
    std::string name;
    int version;
    ParamDir splitDir;
    std::vector<RuledSegment> segments;
    Vec3Arr origMeshVerts;
    FaceArr origMeshFaces;
};

// 单格直纹面拟合（指定方向），内部含准线优化（最小二乘 + lambda 正则）。
RuledCellFit fitCellRuled(const SurfaceWrapper& surf,
                          double u0, double u1, double v0, double v1,
                          ParamDir dir,
                          int nUSamples, int nVSamples,
                          int nRibs, double lambda);

// 单格直纹面拟合（自动方向）：U、V 各拟合一次，取 maxError 更小者。
RuledCellFit fitCellRuledAuto(const SurfaceWrapper& surf,
                              double u0, double u1, double v0, double v1,
                              int nUSamples, int nVSamples,
                              int nRibs, double lambda);

// 等距分段直纹面拟合：沿 splitDir 均分 numSegments 段，逐段调用 fitCellRuled。
RuledResult fitRuledSegments(const SurfaceWrapper& surf,
                              int numSegments,
                              ParamDir splitDir,
                              const std::vector<ParamDir>& directrixDirs,
                              int nUSamples,
                              int nVSamples,
                              int version,
                              const std::string& name);

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

bool exportErrorsCSV(const std::string& path,
                     const std::vector<RuledResult>& results);

bool exportMetaJSON(const std::string& path,
                    const std::string& file1, const std::string& file2,
                    const std::vector<RuledResult>& results);

} // namespace simple
