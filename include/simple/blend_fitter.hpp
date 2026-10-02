#pragma once

#include "common.hpp"
#include "surface_wrapper.hpp"
#include "grid_fitter.hpp"

namespace simple {

// 过渡面（单位分解硬连续）后处理配置。
// 每个格胞对应一个直纹面块 R_i，在其「邻域」上定义 C^2 的基函数 phi_i
// （恒一区域 = 整个格胞，向外以过渡带宽度平滑衰减到 0），再用单位分解
//   S(u,v) = sum_i phi_i(u,v) * R_i(u,v) / sum_i phi_i(u,v)
// 把所有直纹面块混合成全局 C^2 连续的曲面。连续性由构造保证（硬连续）。
//
// 说明：该构造产生的是「全局」曲面 S，并非一块独立的「过渡面」——
// 在未被其他邻域覆盖的内缩核心区内 S 严格等于该格直纹面 R_i，仅在
// 网格线两侧的窄过渡带内才是相邻直纹面的混合。因此可视化不单独画整张
// 混合曲面（与直纹格几乎重合），而是画出内缩核心区边界线来表征过渡带。
struct BlendConfig {
    double bandWidth = 0.15;   // 过渡带宽度：相对最小格胞宽度的比例（外推量，docx 建议 10%~15%）
    int continuity = 2;        // 连续性阶数 C^continuity（0=C⁰线性, 1=C¹, 2=C², 3=C³）
    int nSampleU = 160;        // 拟合误差采样 U 向点数
    int nSampleV = 80;         // 拟合误差采样 V 向点数
};

// 全局混合曲面的拟合误差（相对原曲面）
struct BlendResult {
    double maxError = 0.0;     // 与原曲面最大误差（mm）
    double rmsError = 0.0;     // 与原曲面 RMS 误差（mm）
};

// 在参数 (u,v) 处求全局混合曲面。内缩核心区严格等于直纹面，过渡区内为相邻直纹面的 C^continuity 混合。
Vec3 evalBlend(const GridResult& gr, const BlendConfig& cfg, double u, double v);

// 计算全局混合曲面相对原曲面的拟合误差。
BlendResult computeBlendError(const SurfaceWrapper& surf, const GridResult& gr,
                              const BlendConfig& cfg);

// 导出各格「内缩核心区」边界折线到 VTK。
// 内缩核心区内 S 严格等于该格直纹面；核心边界线与网格线(seam)之间的条带即为过渡带。
bool exportTransitionBandVTK(const std::string& path,
                             const GridResult& gr,
                             const BlendConfig& cfg,
                             int nPerEdge = 40);

// 导出「完整曲面模型」到 JSON，供外部（Python/CasADi）精确重建 C^k 曲面：
// 网格 (uEdges/vEdges/fitDir)、每格直纹准线 B 样条（控制点+节点+次数）、
// 单位分解函数参数（bandWidth/continuity）。保证曲面被完整、可微地导出。
bool exportSurfaceModelJson(const std::string& path,
                            const GridResult& gr,
                            const BlendConfig& cfg);

// 导出全局 C^k 曲面网格 OBJ（供 UI 可视化拟合结果）。
bool exportBlendSurfaceOBJ(const std::string& path,
                           const GridResult& gr,
                           const BlendConfig& cfg,
                           int nU = 160, int nV = 80);

} // namespace simple
