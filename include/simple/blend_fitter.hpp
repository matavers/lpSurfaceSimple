#pragma once

#include "common.hpp"
#include "surface_wrapper.hpp"
#include "grid_fitter.hpp"

namespace simple {

// 缩回 + 过渡面（blend）后处理配置
struct BlendConfig {
    double fMax = 0.3;        // 单格最大内缩比例（相对该边跨度）
    int nAlong = 40;          // 沿接缝方向采样数
    int nAcross = 10;         // 跨界方向采样数
    int nCellGrid = 16;       // trimmed 格子误差采样网格（每向）
    int nMeshRes = 20;        // trimmed 格子导出网格分辨率（每向）
    int nSearch = 24;         // δ 一维搜索网格点数
    int maxPasses = 4;        // 坐标下降轮数
    bool useMaxError = true;  // true=用 max 误差，false=用 RMS
};

// 过渡条（一条内部接缝上的 Hermite/Coons 过渡面）
struct BlendStrip {
    int r = 0, c = 0;
    bool isVertical = true;   // true=竖直缝(u=const)，false=水平缝(v=const)
    bool valid = true;        // false=两侧 f 均为 0，退化零宽，跳过
    Vec3Arr meshVerts;
    FaceArr meshFaces;
    double maxError = 0.0;
    double rmsError = 0.0;
};

// 角点过渡 patch（内部网格顶点处，双线性 Coons）
struct BlendCorner {
    int r = 0, c = 0;         // 顶点 (c,r)
    Vec3Arr meshVerts;
    FaceArr meshFaces;
    double maxError = 0.0;
    double rmsError = 0.0;
};

// trimmed 格子（缩回后的直纹面）
struct TrimmedCell {
    int row = 0, col = 0;
    double f = 0.0;
    double u0 = 0, u1 = 0, v0 = 0, v1 = 0;   // trimmed 参数域
    Vec3Arr meshVerts;
    FaceArr meshFaces;
    double maxError = 0.0;
    double rmsError = 0.0;
};

struct BlendResult {
    std::vector<TrimmedCell> cells;
    std::vector<BlendStrip> strips;
    std::vector<BlendCorner> corners;
    double totalMaxError = 0.0;
    double totalRmsError = 0.0;
};

// 按给定每格内缩量 f 构建「trimmed 格子 + 过渡条」复合面
BlendResult buildBlend(const SurfaceWrapper& surf, const GridResult& gr,
                       const std::vector<double>& f, const BlendConfig& cfg);

// 优化统一内缩量 f，最小化复合面拟合误差（一维搜索）
std::vector<double> optimizeTrim(const SurfaceWrapper& surf, const GridResult& gr,
                                 const BlendConfig& cfg);

// 导出 trimmed 格子与过渡条 OBJ
bool exportBlendOBJs(const std::string& outDir, const std::string& prefix,
                     const BlendResult& br);

} // namespace simple
