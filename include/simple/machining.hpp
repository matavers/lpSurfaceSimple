#pragma once
#include <string>
#include <vector>
#include "common.hpp"

namespace simple {

struct MachiningConfig {
    double feed = 500.0;          /* mm/min */
    double tool_r = 5.0;          /* 侧铣锥度刀半径 mm */
    double ball_r = 5.0;          /* 球头刀半径 mm */
    double scallop = 0.1;         /* 点铣残留高度 mm */
    double twist_limit = 2.0;     /* 可展判定阈值（法向扭转角 deg） */
    double overhead = 4.0;        /* 侧铣每区域进退刀时间 s */
    double point_overhead = 10.0; /* 点铣整体进退刀时间 s */
};

struct MachiningSummary {
    bool ok = false;
    std::string errorMsg;
    int numPatches = 0;
    int flankRegions = 0;
    double flankCut = 0.0;      /* 侧铣切削 s */
    double flankOverhead = 0.0; /* 侧铣非切削 s */
    double flankTotal = 0.0;    /* 侧铣总 s */
    double pointCut = 0.0;      /* 点铣切削 s */
    double pointTotal = 0.0;    /* 点铣总 s */
    double speedup = 0.0;       /* 提速比 */
    double totalArea = 0.0;     /* 拟合面面积 mm² */
    double originalArea = 0.0;  /* 原曲面面积 mm² */
    double maxFlankErr = 0.0;   /* 严谨侧铣点-轴距离最大残差 mm */
    double toolAxisDiscBefore = 0.0; /* 光顺前相邻刀轴最大夹角 deg */
    double toolAxisDiscAfter = 0.0;  /* 光顺后相邻刀轴最大夹角 deg */
    double elapsedSec = 0.0;    /* 刀轨计算耗时 s */
};

/* 读输入目录下的 *_params.txt + blade*_mesh.obj，生成侧铣/点铣刀轨
   （VTK/CSV/DXF），写 summary.json，返回加工时间对比。 */
MachiningSummary computeToolpath(const std::string& inputDir,
                                 const std::string& outputDir,
                                 const MachiningConfig& cfg);

} // namespace simple
