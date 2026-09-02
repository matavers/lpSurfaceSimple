#pragma once

#include "common.hpp"
#include "surface_wrapper.hpp"
#include "ruled_fitter.hpp"
#include "planar_fitter.hpp"

namespace simple {

// 网格中的一个格子（表格单元）。当前参数域为轴对齐矩形，
// 未来倾斜分割时改为 4 个 (u,v) 角点，拟合/采样逻辑不变。
struct GridCell {
    int row = 0;
    int col = 0;
    double u0 = 0, u1 = 0, v0 = 0, v1 = 0;
    ParamDir fitDir = ParamDir::U;   // 直纹面：全局一致方向；平面：无效
    PlaneCellFit plane;
    RuledCellFit ruled;
    double maxError = 0.0;
    double rmsError = 0.0;
    double twist = 0.0;      // 母线法向扭转角最大值(度)
    bool developable = true; // twist ≤ devTol，可侧铣；否则卷曲区，点铣
};

struct GridConfig {
    int nSplitU = 2;          // U 向分割次数（竖直分割）→ 列数 = nSplitU+1
    int nSplitV = 2;          // V 向分割次数（水平分割）→ 行数 = nSplitV+1
    double tolerance = 0.1;   // 容差（目标最大误差，单位 mm）
    int maxDepth = 1000;      // 最大细分步数（安全上限，正常由 tolerance 控制）
    int maxCells = 10000;     // 格子数量安全上限
    int nUSamples = 50;
    int nVSamples = 10;
    int nRibs = 20;           // 直纹面准线优化的 rib 采样数
    double lambda = 1.0;      // 直纹面准线正则强度
    double devTol = 2.0;      // 可展性阈值(度)：仅用于统计 developableCount（twist ≤ devTol），不参与细分
    bool noRefine = false;    // true=固定分片（不自适应细分），只按 nSplitU×nSplitV 均匀网格
};

struct GridResult {
    std::string name;
    int nRows = 0;
    int nCols = 0;
    std::vector<double> uEdges;   // nCols+1
    std::vector<double> vEdges;   // nRows+1
    std::vector<GridCell> cells;  // 行优先，size = nRows*nCols
    ParamDir fitDir = ParamDir::V;   // 全局一致直纹面方向（每次二分后重试算）
    double maxError = 0.0;        // 细分后实际达到的最大误差（mm）
    bool toleranceMet = true;     // 是否满足 tolerance
    double maxTwist = 0.0;        // 细分后最大母线扭转角（度）
    int developableCount = 0;     // 可展格子数（twist ≤ devTol）
};

// 二维网格分区拟合主入口：先均匀 M×N，再按容差整行/整列自适应细分。
GridResult fitGridRuled(const SurfaceWrapper& surf,
                        const GridConfig& cfg,
                        const std::string& name);

GridResult fitGridPlanar(const SurfaceWrapper& surf,
                         const GridConfig& cfg,
                         const std::string& name);

// 导出：每个格子一个 OBJ（平面或直纹面），平面附 _desc.txt、直纹面附 _params.txt
bool exportGridOBJs(const std::string& outDir, const std::string& prefix,
                    const GridResult& gr, bool planar);

// 导出网格分割线（水平 + 竖直）为 VTK PolyData（LINES），供可视化树 Grid 节点显示
bool exportGridLinesVTK(const std::string& path,
                        const SurfaceWrapper& surf,
                        const GridResult& gr,
                        int nSamplesPerLine = 60);

// 生成 meta.json 内容（供 DLL 的 metaJson 字段与 GUI 读取）
std::string buildGridMetaJson(const std::vector<GridResult>& results,
                              const std::string& mode,
                              const std::string& file1,
                              const std::string& file2);

} // namespace simple
