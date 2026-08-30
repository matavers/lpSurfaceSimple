# ruledSurfaceFitting — 直纹面/平面分片拟合 API 技术文档

## 1. 算法概述

本库 `ruledSurfaceFitting` 实现基于 NURBS 曲面的**直纹面分片拟合**与**平面分片拟合**算法，用于航空发动机叶片曲面的几何简化与可加工性改造。将自由曲面转化为分片直纹面（或平面）集合后，可用侧刃铣削（侧铣）替代点铣。

### 1.1 支持的导入文件类型

| 格式 | 扩展名 | 标准 |
|------|--------|------|
| **STEP** | `.step`, `.stp` | ISO 10303 AP203/AP214 |
| **IGES** | `.igs`, `.iges` | ANSI US PRO/IPO-100 |

### 1.2 接口总览

动态库导出 **4 个接口**（固定分段与自适应细分相互独立），头文件 `ruledSurfaceFitting.h` 暂时只暴露 2 个简化接口（面向客户交付）：

| 函数 | 暴露位置 | 拟合类型 | 分段方式 |
|------|---------|---------|---------|
| `ruled_fitting(RuledFitConfig)` | 仅 DLL（暂不暴露） | 直纹面 | 井字形网格自适应细分（按容差） |
| `plane_fitting(PlaneFitConfig)` | 仅 DLL（暂不暴露） | 平面 | 自适应细分（按容差，每细分面拟合为平面） |
| `ruled_fitting_simple(dir, dir)` | 头文件 | 直纹面 | 固定 3 等分（不分细分） |
| `plane_fitting_simple(dir, dir)` | 头文件 | 平面 | 固定 3 等分（不分细分） |

### 1.3 固定分段流程（简化接口）

```
输入目录遍历 STEP/IGES
  → 单文件自动识别压力面/吸力面（法向聚类 + 弦向截面分割）
  → 每面沿叶片高度方向（U 向）等距 3 等分
  → 每段：直纹面拟合（准线优化） 或 平面拟合（PCA）
  → 导出 .obj 网格 + .txt 参数
```

### 1.4 自适应细分流程（3 参数接口）

```
单文件 CAD 导入
  → 自动识别压力面/吸力面
  → 直纹面：井字形网格（整行/整列二分）按容差细分
  → 平面：沿高度方向自适应二分，每个细分面拟合为平面
  → 导出 .obj 网格 + .txt 参数
```

### 1.5 直纹面算法原理

对每段子曲面参数域 `[uSeg0,uSeg1]×[vSeg0,vSeg1]`（准线沿 u、母线沿 v）：

1. **初始准线**：取 v=vSeg0、v=vSeg1 两条等参线作为上下准线 C₀、C₁
2. **rib 采样**：沿准线方向固定参数，在母线方向采样 nRibs 个中间点构成"肋条"
3. **最小二乘优化**：对每条母线独立求解 2×2 线性方程组，调整 C₀[i]、C₁[i] 使肋条点到直纹面的距离平方和最小，并用 λ 正则化约束准线不偏离初始边界
4. **直纹面**：`S(p,t) = (1−t)·C₀(p) + t·C₁(p)`，p 沿准线、t 沿母线

### 1.6 平面算法原理

对每段子曲面采样点做主成分分析（PCA）：协方差矩阵最小特征值对应法向量 n，过质心 c，平面方程 `n·(x−c)=0`；误差 = 点到平面距离 `|n·(P−c)|`。

---

## 2. 接口函数

### 2.1 `ruled_fitting` — 直纹面自适应（单文件，仅 DLL）

```c
RULED_API RuledFittingResult* ruled_fitting(const RuledFitConfig* config);
```

加载单个 STEP/IGES 文件 → 自动识别压力面/吸力面 → 每面按容差做**井字形网格自适应细分**。仅输出 `.obj` 与 `.txt`。

### 2.2 `plane_fitting` — 平面自适应（单文件，仅 DLL）

```c
RULED_API RuledFittingResult* plane_fitting(const PlaneFitConfig* config);
```

加载单个 STEP/IGES 文件 → 自动识别压力面/吸力面 → 每面按容差**自适应细分**，每个细分面拟合为平面。仅输出 `.obj` 与 `.txt`。

### 2.3 `ruled_fitting_simple` — 直纹面固定 3 等分（批处理）

```c
RULED_API RuledFittingResult* ruled_fitting_simple(const char* inputDir, const char* outputDir);
```

遍历输入目录下的 STEP/IGES 文件，对每个文件自动识别并**沿叶片高度方向固定 3 等分直纹面拟合**（无容差判断、不细分）。输出文件名以输入文件名 + `_pressure`/`_suction` 为前缀。

### 2.4 `plane_fitting_simple` — 平面固定 3 等分（批处理）

```c
RULED_API RuledFittingResult* plane_fitting_simple(const char* inputDir, const char* outputDir);
```

遍历输入目录下的 STEP/IGES 文件，对每个文件自动识别并**沿叶片高度方向固定 3 等分平面拟合**（无容差判断、不细分）。

### 2.5 `free_result`

```c
RULED_API void free_result(RuledFittingResult* result);
```

释放上述函数分配的内存。

---

## 3. 参数结构体

### 3.1 `RuledFitConfig` / `PlaneFitConfig` — 拟合配置（仅 DLL）

| 字段 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `inputPath` | `const char*` | — | 输入叶片 STEP/IGES 文件路径（单文件） |
| `outputDir` | `const char*` | — | 输出目录 |
| `tolerance` | `double` | 0.1（≤0 时） | 容差 (mm)，最大误差阈值 |

其余拟合参数（采样数 50×10、rib 数 20、正则强度 1.0、分割方向 U=叶片高度、准线方向 V=弦向等）均由算法内部取默认值。

### 3.2 `RuledFittingResult` — 返回结果

| 字段 | 类型 | 说明 |
|------|------|------|
| `errorCode` | `RuledErrorCode` | 错误码，0=成功 |
| `errorMsg[256]` | `char[256]` | 错误描述 |
| `numSurfaces` | `int` | 处理的曲面数（压力面/吸力面，最多 2） |
| `surfaces[2]` | `RuledSurfaceResult[2]` | 每面拟合结果 |
| `metaJson[4096]` | `char[4096]` | 结果 JSON 摘要 |

### 3.3 `RuledSurfaceResult` — 单面结果

| 字段 | 类型 | 说明 |
|------|------|------|
| `name[64]` | `char[64]` | 曲面名称（Pressure / Suction） |
| `maxError` | `double` | 该面整体最大误差（mm） |
| `numSegments` | `int` | 实际面片数（固定 3 或细分后的单元数） |
| `segments[64]` | `RuledSegmentResult[64]` | 各面片误差（前 64 个） |

### 3.4 `RuledSegmentResult` — 单面片结果

| 字段 | 类型 | 说明 |
|------|------|------|
| `segmentIndex` | `int` | 面片索引（0-based） |
| `maxError` | `double` | 最大距离误差 |
| `rmsError` | `double` | 均方根误差 |

### 3.5 `RuledErrorCode` — 错误码

| 值 | 含义 |
|----|------|
| `RULED_OK (0)` | 成功 |
| `RULED_ERR_FILE_NOT_FOUND (1)` | 文件未找到 |
| `RULED_ERR_STEP_READ_FAILED (2)` | CAD 文件读取失败 |
| `RULED_ERR_NO_VALID_FACE (3)` | 未找到有效面 |
| `RULED_ERR_EXPORT_FAILED (4)` | 结果导出失败 |
| `RULED_ERR_INVALID_PARAMS (5)` | 参数无效 |

---

## 4. 输出文件

| 拟合类型 | 文件 | 格式 | 说明 |
|---------|------|------|------|
| 直纹面 | `{prefix}_segN.obj` | Wavefront OBJ | 第 N 个面片的直纹面四边形网格 |
| 直纹面 | `{prefix}_segN_params.txt` | TXT | 第 N 个面片的优化准线采样点 |
| 直纹面 | `{prefix}_segN_bspline.txt` | TXT | 第 N 个面片优化准线的 B 样条参数（次数、节点向量、控制点） |
| 平面 | `{prefix}_planeN.obj` | Wavefront OBJ | 第 N 个面片的平面网格 |
| 平面 | `{prefix}_planeN_desc.txt` | TXT | 第 N 个面片的平面质心 + 法向 |

`{prefix}` 为 `pressure` / `suction`（简化接口为 `<文件名>_pressure` / `<文件名>_suction`）。仅导出 `.obj` 与 `.txt` 两类文件，无其他内部产物。

### 4.1 直纹面准线参数 TXT（`*_params.txt`）

```
[C0]
n = 50
297.725952 12.023527 206.297000
...

[C1]
n = 50
...

[mapping]
description = identity, C0 and C1 share the same parameterization (i-th sample pair forms a ruling)
```

- `[C0]` / `[C1]`：上下两条**优化后**准线的采样点坐标（每行 x y z）
- 直纹面由 `S(i,t) = (1−t)·C0[i] + t·C1[i]` 定义，第 i 对采样点构成一条母线

### 4.2 直纹面准线 B 样条 TXT（`*_bspline.txt`）

将优化后的准线采样点按**共同弦长参数化**拟合为 3 次 B 样条曲线后导出，下游（NX/UG/CAM）可直接读取重建，无需预处理：

```
nSamples = 50

[C0]
degree = 3
rational = no
nbPoles = 50
poles (x y z w):
  322.001752 -9.372975 206.297000 1.000000
  ...
nbKnots = 48
knots: 0.00000000 0.02136009 ... 1.00000000
multiplicities: 4 1 1 ... 1 4

[C1]
...
[mapping]
description = identity, C0(u) and C1(u) share the same parameterization (i-th point/control pair forms a ruling)
```

- `degree`：曲线次数（3）
- `rational`：是否有理（当前均为 `no`，权重全 1）
- `poles (x y z w)`：控制点坐标及权重
- `knots`：节点向量（去重后）；`multiplicities`：对应节点重数（首尾重数为 `degree+1`，为 clamped 曲线）
- `[C0]` / `[C1]` 共享同一节点向量，`S(u,v) = (1−v)·C0(u) + v·C1(u)` 为直纹面参数方程

### 4.3 平面描述 TXT（`*_desc.txt`）

```
centroid = 302.034950 8.112750 236.289690
normal = -0.834570 -0.550860 -0.006960
```

- `centroid`：平面质心坐标
- `normal`：平面单位法向量，平面方程 `n·(x−centroid)=0`

---

## 5. 使用示例

### 5.1 C 语言（简化接口）

```c
#include "ruledSurfaceFitting.h"

int main() {
    RuledFittingResult* res = ruled_fitting_simple("./input", "./output");
    if (res->errorCode == RULED_OK) {
        for (int s = 0; s < res->numSurfaces; ++s)
            printf("%s: %d pieces, maxError=%.4f\n",
                   res->surfaces[s].name,
                   res->surfaces[s].numSegments,
                   res->surfaces[s].maxError);
    }
    free_result(res);
    return 0;
}
```

### 5.2 Python ctypes 调用（简化接口）

```python
import ctypes

lib = ctypes.CDLL("./ruledSurfaceFitting.dll")
lib.ruled_fitting_simple.restype = ctypes.c_void_p
lib.ruled_fitting_simple.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
lib.free_result.argtypes = [ctypes.c_void_p]

ptr = lib.ruled_fitting_simple(b"./input", b"./output")
lib.free_result(ptr)
```

---

## 6. 构建与依赖

- **OpenCASCADE 8.0.0**：静态链接（TKernel、TKMath、TKG2d、TKG3d、TKGeomBase、TKBRep、TKGeomAlgo、TKTopAlgo、TKPrim、TKBO、TKShHealing、TKBool、TKHLR、TKHelix、TKFillet、TKOffset、TKFeat、TKMesh、TKXMesh、TKExpress、TKService、TKV3d、TKOpenGl、TKMeshVS、TKCDF、TKLCAF、TKCAF、TKBinL、TKXmlL、TKBin、TKXml、TKStdL、TKStd、TKTObj、TKBinTObj、TKXmlTObj、TKVCAF、TKDE、TKXSBase、TKDESTEP、TKXCAF、TKDEIGES、TKDESTL、TKDEVRML、TKRWMesh、TKDECascade、TKBinXCAF、TKXmlXCAF、TKDEOBJ、TKDEPLY）
- **Eigen 3**（bundled in VTK 9.4）

```bash
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```

产物：
- `build/Release/ruledSurfaceFitting.dll` + `.lib` — 动态库（OCCT 静态内嵌，无 OCCT 运行时 DLL 依赖）
- `build/Release/simple.exe` — 命令行工具

命令行模式：`ruled`（直纹面固定）、`ruled-adaptive`（直纹面井字形自适应）、`planar`（平面固定）、`planar-adaptive`（平面自适应）。

---

## 7. 误差计算方法

对每个面片参数域，沿准线方向取 `nAlong` 点、母线方向取 `nAcross` 点构成采样点集 {P_k}，在每个面片网格上搜索最近点：

- `maxError = max_k dist(P_k, S)` — 最大局部偏差（用于容差判断）
- `rmsError = sqrt(1/N · Σ_k dist(P_k, S)²)` — 整体拟合质量

---

## 8. 变更记录

| 版本 | 日期 | 变更 |
|------|------|------|
| v4.1 | 2026-08-30 | 直纹面新增 `_bspline.txt` 导出：将优化后准线采样点按共同弦长参数化拟合为 3 次 B 样条并导出次数/节点向量/控制点（修复原 `_bspline.txt` 导出的是未优化等参线且未声明导致编译失败的问题）。GUI 同步适配：树中直纹面段下新增 Directrix C0/C1 节点，可视化优化准线曲线，点击节点展示 B 样条参数。 |
| v4.0 | 2026-08-29 | 固定分段与自适应细分拆分为独立接口：新增 `plane_fitting`（平面自适应）、`plane_fitting_simple`（平面固定 3 段）、`ruled_fitting`（直纹面井字形自适应，纯网格）、`ruled_fitting_simple`（直纹面固定 3 段）。3 段分割默认方向改为叶片高度（U 向）。头文件暂只暴露两个简化接口。 |
| v3.0 | 2026-08-28 | 单文件 API：`ruled_fitting`（3 字段 RuledFitConfig）+ `ruled_fitting_simple`。仅导出 obj/txt。txt 改为导出优化后准线采样点（修复不同段 txt 相同问题）。删除旧 6 个接口。 |
| v2.0 | 2026-08-09 | 新增 4 个全自动函数（pressure/suction × ruled/plane）。 |
| v1.0 | 2026-07 | 初始版本：双文件手动模式。 |
