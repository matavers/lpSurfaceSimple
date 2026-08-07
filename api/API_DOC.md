# ruledSurfaceFitting — 直纹面/平面逼近算法 API 技术文档

## 1. 算法概述

本库 `ruledSurfaceFitting` 实现了基于 NURBS 曲面的**分段直纹面拟合**和**分段平面拟合**两种算法，用于航空发动机叶片曲面的几何简化。

### 1.1 支持的导入文件类型

| 格式 | 扩展名 | 标准 |
|------|--------|------|
| **STEP** | `.step`, `.stp` | ISO 10303 AP203/AP214 |
| **IGES** | `.igs`, `.iges` | ANSI US PRO/IPO-100 |

支持包含多个面的文件——可在可视化 UI 中预览各面，选择任意两个面进行处理。

### 1.2 算法流程

```
CAD文件导入(STEP/IGES) → 面裁剪域提取 → 参数域等距分段 →
  [直纹面模式]   等参线准线提取 → rib最小二乘优化 → 直纹面网格 + NURBS参数TXT导出
  [平面模式]     PCA最佳拟合平面  → 平面四边形     + 描述JSON导出
→ 误差CSV报告 → 元数据JSON摘要
```

### 1.3 两种模式

| 模式 | DLL 函数 | 描述 |
|------|---------|------|
| 直纹面逼近 | `ruled_surface_fitting` | 将曲面沿参数域分段，每段用双准线+rib优化的直纹面拟合 |
| 平面逼近 | `plane_surface_fitting` | 将曲面沿参数域分段，每段用PCA最佳拟合平面逼近 |

### 1.4 直纹面算法原理

对于每段子曲面，在参数域裁剪范围 `[uSeg0, uSeg1] × [vSeg0, vSeg1]` 内：

1. **初始准线**：沿分割方向取等参线作为上下两条边界曲线 C₀、C₁
2. **rib采样**：沿准线方向固定参数值，在交叉方向上采样 nRibs 个中间点构成"肋条"
3. **最小二乘优化**：对每条ruling线独立求解 6 维线性方程组，调整采样点 C₀[i] 和 C₁[i] 使肋条点到直纹面的距离平方和最小，同时用 λ 正则化约束准线不偏离原始边界曲线
4. **直纹面**：`S(p, t) = (1-t)·C₀(p) + t·C₁(p)`，p 沿准线方向，t 沿交叉方向

### 1.5 方向组合（直纹面模式）

| 分割方向 | 准线方向 | 适用场景 |
|---------|---------|---------|
| V | V | 叶片高度方向分段，水平直纹（默认） |
| V | U | 高度方向分段，纵向直纹 |
| U | U | 周向分段，纵向直纹 |
| U | V | 周向分段，水平直纹 |

### 1.6 平面逼近算法原理

对每段子曲面，在参数域裁剪范围 `[uSeg0, uSeg1] × [vSeg0, vSeg1]` 内：

1. **均匀采样**：取 N = nUSamples × nVSamples 个曲面点 {P_k}
2. **PCA 协方差分析**：计算质心 c = (1/N) Σ P_k，协方差矩阵 C = (1/N) Σ (P_k−c)(P_k−c)^T
3. **特征分解**：C 为 3×3 实对称矩阵，最小特征值 λ_min 对应的单位特征向量 n 即为最佳拟合平面的法向
4. **平面方程**：`n · (x − centroid) = 0`，即过质心、法向为 n 的平面
5. **误差定义**：对每个采样点 P_k，`dist(P_k) = |n · (P_k − centroid)|`，取最大距离 maxError 和均方根 rmsError
6. **网格生成**：将分段的 4 个角点投影到平面上构成四边形

---

## 2. 接口函数

### 2.1 `ruled_surface_fitting` — 直纹面逼近

```c
RULED_API RuledFittingResult* ruled_surface_fitting(const RuledConfig* config);
```

执行完整直纹面拟合管线：加载 STEP/IGES 文件 → 面裁剪 → V/U 方向等距分段 → 等参线准线提取 → rib 最小二乘优化 → OBJ 网格 + NURBS 参数 TXT 导出。

### 2.2 `plane_surface_fitting` — 平面逼近

```c
RULED_API RuledFittingResult* plane_surface_fitting(const PlanarConfig* config);
```

执行完整平面拟合管线：加载 STEP/IGES 文件 → 面裁剪 → V/U 方向等距分段 → PCA 最佳拟合平面 → OBJ 网格 + 平面描述 JSON 导出。每段拟合平面由其质心 `centroid` 和单位法向量 `normal` 定义，平面方程 `normal · (x − centroid) = 0`。

### 2.3 `free_result`

```c
RULED_API void free_result(RuledFittingResult* result);
```

释放上述两个函数分配的内存。

---

## 3. 参数结构体

### 3.1 `RuledConfig` — 直纹面逼近配置

| 字段 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `stepFile1` | `const char*` | — | 第一个 CAD 文件路径（STEP/IGES） |
| `stepFile2` | `const char*` | — | 第二个 CAD 文件路径 |
| `outputDir` | `const char*` | `"."` | 输出目录 |
| `nUSamples` | `int` | 50 | 沿准线方向采样点数 |
| `nVSamples` | `int` | 10 | 交叉方向采样点数 |
| `nRibs` | `int` | 20 | 每条 ruling 线肋条点数 |
| `lambda` | `double` | 1.0 | 正则化强度 |
| `splitDir[2]` | `RuledDirection[2]` | `{V,V}` | 每面分割方向 |
| `directrixDirs[2][10]` | `int[2][10]` | `{{V,V,V},{V,V,V}}` | 每段准线方向 |
| `numDirectrixDirs[2]` | `int[2]` | `{3,3}` | 每面段数 |
| `faceIdx[2]` | `int[2]` | `{-1,-1}` | 面索引，-1 自动选最大面 |

### 3.2 `PlanarConfig` — 平面逼近配置

| 字段 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `stepFile1` | `const char*` | — | 第一个 CAD 文件路径 |
| `stepFile2` | `const char*` | — | 第二个 CAD 文件路径 |
| `outputDir` | `const char*` | `"."` | 输出目录 |
| `nUSamples` | `int` | 50 | U 方向采样点数 |
| `nVSamples` | `int` | 10 | V 方向采样点数 |
| `splitDir[2]` | `RuledDirection[2]` | `{V,V}` | 每面分割方向 |
| `faceIdx[2]` | `int[2]` | `{-1,-1}` | 面索引，-1 自动选最大面 |

### 3.3 `RuledDirection` — 方向枚举

| 值 | 含义 |
|----|------|
| `RULED_DIR_U (0)` | U 参数方向 |
| `RULED_DIR_V (1)` | V 参数方向 |

### 3.4 `RuledFittingResult` — 返回结果

| 字段 | 类型 | 说明 |
|------|------|------|
| `errorCode` | `RuledErrorCode` | 错误码，0=成功 |
| `errorMsg[256]` | `char[256]` | 错误描述 |
| `numSurfaces` | `int` | 处理的曲面数（固定为2） |
| `surfaces[2]` | `RuledSurfaceResult[2]` | 每面拟合结果 |
| `metaJson[2048]` | `char[2048]` | 结果 JSON 摘要 |

### 3.5 `RuledSurfaceResult` — 单面结果

| 字段 | 类型 | 说明 |
|------|------|------|
| `name[64]` | `char[64]` | 曲面名称 |
| `numSegments` | `int` | 段数 |
| `segments[10]` | `RuledSegmentResult[10]` | 每段结果 |

### 3.6 `RuledSegmentResult` — 单段结果

| 字段 | 类型 | 说明 |
|------|------|------|
| `segmentIndex` | `int` | 段索引（0-based） |
| `maxError` | `double` | 最大距离误差 |
| `rmsError` | `double` | 均方根误差 |

### 3.7 `RuledErrorCode` — 错误码

| 值 | 含义 |
|----|------|
| `RULED_OK (0)` | 成功 |
| `RULED_ERR_FILE_NOT_FOUND (1)` | 文件路径为空 |
| `RULED_ERR_STEP_READ_FAILED (2)` | CAD 文件读取失败 |
| `RULED_ERR_NO_VALID_FACE (3)` | 未找到有效面 |
| `RULED_ERR_EXPORT_FAILED (4)` | 结果导出失败 |
| `RULED_ERR_INVALID_PARAMS (5)` | 参数为空 |

---

## 4. 输出文件总览

| 文件 | 格式 | 模式 | 说明 |
|------|------|------|------|
| `blade1_mesh.obj` | Wavefront OBJ | 通用 | 面 1 原始裁剪三角网格 |
| `blade2_mesh.obj` | Wavefront OBJ | 通用 | 面 2 原始裁剪三角网格 |
| `blade1_segN.obj` | Wavefront OBJ | 直纹面 | 第 N 段直纹面四边形网格 |
| `blade1_segN_params.txt` | TXT | 直纹面 | 第 N 段 NURBS 参数方程（度数+控制点+节点向量+权重+映射） |
| `blade1_planeN.obj` | Wavefront OBJ | 平面 | 第 N 段拟合平面四边形 |
| `blade1_planeN_desc.txt` | TXT | 平面 | 第 N 段平面质心 + 法向量 |
| `errors.csv` | CSV | 通用 | 所有段误差指标表 |
| `meta.json` | JSON | 通用 | 元数据与误差摘要 |

### 4.1 准线参数方程 TXT（`*_params.txt`）— 仅直纹面模式

```
nSamples = 50

[C0]
degree = 3
rational = no
nbPoles = 25
poles (x y z w):
  322.19 -8.96 206.30 1.0
  ...
knots: 0.0 0.09375 0.125 ... 1.0
multiplicities: 4 1 1 ... 4

[C1]
...

[mapping]
description = identity, C0(u) and C1(u) share same U-parameterization from surface iso-curves
```

- `degree`：NURBS 曲线次数
- `poles`：控制点 (x, y, z) 及权重 w，每行一个
- `knots`：节点向量
- `multiplicities`：各节点的重复度
- `mapping`：C₀→C₁ 的参数对应关系为恒等映射

### 4.2 平面描述 TXT（`*_desc.txt`）— 仅平面模式

```
centroid = cx cy cz
normal = nx ny nz
corner0 = x0 y0 z0
corner1 = x1 y1 z1
corner2 = x2 y2 z2
corner3 = x3 y3 z3
```

拟合平面方程为 `normal · (x − centroid) = 0`，由4个角点围成边界四边形。

---

## 5. 使用示例

### 5.1 C 语言 — 直纹面逼近

```c
#include "ruledSurfaceFitting.h"

int main() {
    RuledConfig cfg = {0};
    cfg.stepFile1 = "Blade-raw1.STEP";
    cfg.stepFile2 = "Blade-raw2.STEP";
    cfg.outputDir = "./output";
    cfg.nUSamples = 50; cfg.nVSamples = 10;
    cfg.nRibs = 20; cfg.lambda = 1.0;
    cfg.splitDir[0] = RULED_DIR_V;
    cfg.splitDir[1] = RULED_DIR_V;
    cfg.numDirectrixDirs[0] = 3;
    cfg.directrixDirs[0][0] = RULED_DIR_V;
    cfg.directrixDirs[0][1] = RULED_DIR_V;
    cfg.directrixDirs[0][2] = RULED_DIR_V;
    cfg.directrixDirs[1][0] = RULED_DIR_V;
    cfg.directrixDirs[1][1] = RULED_DIR_V;
    cfg.directrixDirs[1][2] = RULED_DIR_V;

    RuledFittingResult* res = ruled_surface_fitting(&cfg);
    if (res->errorCode == RULED_OK) {
        printf("Seg 1 maxErr: %.4f\n",
               res->surfaces[0].segments[1].maxError);
    }
    free_result(res);
    return 0;
}
```

### 5.2 C 语言 — 平面逼近

```c
PlanarConfig cfg = {0};
cfg.stepFile1 = "Blade-raw1.STEP";
cfg.stepFile2 = "Blade-raw2.STEP";
cfg.outputDir = "./output";
cfg.nUSamples = 50; cfg.nVSamples = 10;
cfg.splitDir[0] = RULED_DIR_V;
cfg.splitDir[1] = RULED_DIR_V;

RuledFittingResult* res = plane_surface_fitting(&cfg);
// ... same result handling ...
free_result(res);
```

### 5.3 Python ctypes 调用

```python
import ctypes

lib = ctypes.CDLL("./ruledSurfaceFitting.dll")
lib.ruled_surface_fitting.restype = ctypes.c_void_p
lib.free_result.argtypes = [ctypes.c_void_p]

class RuledConfig(ctypes.Structure):
    _fields_ = [
        ("stepFile1", ctypes.c_char_p),
        ("stepFile2", ctypes.c_char_p),
        ("outputDir", ctypes.c_char_p),
        ("nUSamples", ctypes.c_int),
        ("nVSamples", ctypes.c_int),
        ("nRibs", ctypes.c_int),
        ("lambda", ctypes.c_double),
        ("splitDir", ctypes.c_int * 2),
        ("directrixDirs", (ctypes.c_int * 10) * 2),
        ("numDirectrixDirs", ctypes.c_int * 2),
        ("faceIdx", ctypes.c_int * 2),
    ]

cfg = RuledConfig()
cfg.stepFile1 = b"Blade.igs"
cfg.stepFile2 = b"Blade.igs"
cfg.outputDir = b"./output"
cfg.nUSamples = 50

ptr = lib.ruled_surface_fitting(ctypes.byref(cfg))
# ... read result ...
lib.free_result(ptr)
```

---

## 6. 构建与依赖

### 6.1 依赖库

- **OpenCASCADE 8.0.0** (TKernel, TKMath, TKGeomBase, TKGeomAlgo, TKG2d, TKG3d, TKMesh, TKBRep, TKTopAlgo, TKBO, TKService, TKXSBase, TKDESTEP, TKDEIGES, TKShHealing)
- **Eigen 3** (bundled in VTK 9.4)

### 6.2 构建

```bash
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```

产物：
- `build/Release/ruledSurfaceFitting.dll` + `.lib` — 动态库
- `build/Release/simple.exe` — 命令行工具

### 6.3 部署

使用 DLL 时需确保以下 OCCT 运行时 DLL 在 PATH 中：
`TKernel.dll`, `TKMath.dll`, `TKGeomBase.dll`, `TKGeomAlgo.dll`, `TKG2d.dll`, `TKG3d.dll`, `TKMesh.dll`, `TKBRep.dll`, `TKTopAlgo.dll`, `TKBO.dll`, `TKService.dll`, `TKXSBase.dll`, `TKDESTEP.dll`, `TKDEIGES.dll`, `TKShHealing.dll`, `freetype.dll`, `FreeImage.dll`, `tbb12.dll`, `tbbmalloc.dll`, `jemalloc.dll`

---

## 7. 测试

### 7.1 Python 测试（推荐）

```bash
D:\anaconda\envs\simple\python.exe tests/test_api.py
D:\anaconda\envs\simple\python.exe tests/test_api.py --gui
```

### 7.2 C 测试

```bash
cl /EHsc /I"api" tests\test_api.c /link /LIBPATH:"build\Release" ruledSurfaceFitting.lib
test_api.exe Blade-raw1.STEP Blade-raw2.STEP
```

---

## 8. 误差计算方法

### 8.1 直纹面逼近误差

对每段子曲面参数域 `[uSeg0, uSeg1] × [vSeg0, vSeg1]`：

1. **采样网格**：沿准线方向取 `nAlong` 个点，交叉方向取 `nAcross` 个点，共 N 个采样点 {P_k}
2. **点到直纹面距离**：在直纹面网格上暴力搜索最近点
3. **误差指标**：
   - `maxError = max_k dist(P_k, S)` — 最大局部偏差
   - `rmsError = sqrt(1/N · Σ_k dist(P_k, S)²)` — 整体拟合质量

### 8.2 平面逼近误差

1. **PCA 平面拟合**：协方差矩阵最小特征值对应法向 n
2. **点到平面距离**：`dist(P_k) = |n·(P_k − centroid)|`
3. **误差指标**：同上

---

## 9. 导出文件格式详述

### 9.1 原始网格 OBJ（`*_mesh.obj`）

Wavefront OBJ 格式，包含面裁剪域内的三角网格。

### 9.2 直纹面网格 OBJ（`*_segN.obj`）

四边形网格，由两条准线的采样点经线性插值生成。

### 9.3 准线参数方程 TXT（`*_segN_params.txt`）— 仅直纹面模式

INI格式纯文本，包含完整 NURBS 数学表示：
- `nSamples`：采样点数
- `[C0]` / `[C1]`：上下准线的 degree、rational、nbPoles、poles(x y z w)、knots、multiplicities
- `[mapping]`：参数对应关系（恒等：C₀ 和 C₁ 共享同一 U 参数化）

### 9.4 平面网格 OBJ（`*_planeN.obj`）— 仅平面模式

4 个顶点的四边形（2 个三角形），将分段参数域 4 角点投影到 PCA 最佳拟合平面所得。

### 9.5 平面描述 TXT（`*_planeN_desc.txt`）— 仅平面模式

```
centroid = 317.360873 -12.783161 236.289687
normal = -0.733008 -0.665962 -0.138541
corner0 = 320.774019 -10.244814 206.029191
corner1 = 314.577559 -15.952813 266.252473
corner2 = 310.909618 -11.908336 266.217559
corner3 = 317.347649 -6.551311 206.403261
```
平面方程：`normal · (x − centroid) = 0`，normal 为单位法向量。4个corner为分段参数域四角点到平面上的投影坐标，围成平面边界四边形。

### 9.6 误差表 CSV（`errors.csv`）

```csv
surface,version,segment,mode,maxError,rmsError
Blade-1,0,0,ruled,0.06369,0.13380
```

### 9.7 元数据 JSON（`meta.json`）

```json
{"mode":"ruled","files":["Blade-raw1.STEP","Blade-raw2.STEP"],"surfaces":[...]}
```

> **注意**：可视化 UI 仅加载 `.obj` 文件用于 3D 渲染。`.txt`、`.json` 和 `.csv` 文件不会被 UI 读取，仅供外部数据分析使用。
