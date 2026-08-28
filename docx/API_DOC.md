# ruledSurfaceFitting — 直纹面拟合 API 技术文档

## 1. 算法概述

本库 `ruledSurfaceFitting` 实现基于 NURBS 曲面的**直纹面分片拟合**算法，用于航空发动机叶片曲面的几何简化与可加工性改造。将自由曲面转化为分片直纹面集合后，可用侧刃铣削（侧铣）替代点铣。

### 1.1 支持的导入文件类型

| 格式 | 扩展名 | 标准 |
|------|--------|------|
| **STEP** | `.step`, `.stp` | ISO 10303 AP203/AP214 |
| **IGES** | `.igs`, `.iges` | ANSI US PRO/IPO-100 |

### 1.2 算法流程（主接口）

```
单文件 CAD 导入(STEP/IGES)
  → 法向聚类自动识别压力面/吸力面
  → 弦向截面分割，提取叶盆/叶背参数区间
  → 对每面：等距 3 等分直纹面拟合
       → 若最大误差 < tolerance：直接输出 3 段
       → 若最大误差 ≥ tolerance：井字形网格自适应细分（整行/整列二分，直至满足容差）
  → 每段/每格导出 .obj 网格 + .txt 准线参数
```

### 1.3 直纹面算法原理

对每段子曲面参数域 `[uSeg0,uSeg1]×[vSeg0,vSeg1]`（准线沿 u、母线沿 v）：

1. **初始准线**：取 v=vSeg0、v=vSeg1 两条等参线作为上下准线 C₀、C₁
2. **rib 采样**：沿准线方向固定参数，在母线方向采样 nRibs 个中间点构成"肋条"
3. **最小二乘优化**：对每条母线独立求解 2×2 线性方程组，调整 C₀[i]、C₁[i] 使肋条点到直纹面的距离平方和最小，并用 λ 正则化约束准线不偏离初始边界
4. **直纹面**：`S(p,t) = (1−t)·C₀(p) + t·C₁(p)`，p 沿准线、t 沿母线

---

## 2. 接口函数

### 2.1 `ruled_fitting` — 主接口（单文件自动处理）

```c
RULED_API RuledFittingResult* ruled_fitting(const RuledFitConfig* config);
```

加载单个 STEP/IGES 文件 → 自动识别压力面/吸力面 → 每面先 3 等分拟合；若最大误差小于容差直接输出，否则按容差做井字形网格自适应细分。仅输出 `.obj` 与 `.txt` 文件。

### 2.2 `ruled_fitting_simple` — 简化接口（批处理，固定 3 等分）

```c
RULED_API RuledFittingResult* ruled_fitting_simple(const char* inputDir, const char* outputDir);
```

遍历输入目录下的 STEP/IGES 文件，对每个文件自动识别并**固定 3 等分拟合**（无容差判断、不细分），输出到输出目录。输出文件名以输入文件名 + `_pressure`/`_suction` 为前缀。

### 2.3 `free_result`

```c
RULED_API void free_result(RuledFittingResult* result);
```

释放上述函数分配的内存。

---

## 3. 参数结构体

### 3.1 `RuledFitConfig` — 拟合配置

| 字段 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `inputPath` | `const char*` | — | 输入叶片 STEP/IGES 文件路径（单文件） |
| `outputDir` | `const char*` | — | 输出目录 |
| `tolerance` | `double` | 0.1（≤0 时） | 容差 (mm)，最大误差阈值 |

其余拟合参数（采样数、rib 数、正则强度、分割方向、准线方向等）均由算法内部取默认值。

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
| `numSegments` | `int` | 实际面片数（3 或网格单元数） |
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

| 文件 | 格式 | 说明 |
|------|------|------|
| `{prefix}_segN.obj` | Wavefront OBJ | 第 N 个面片的直纹面四边形网格 |
| `{prefix}_segN_params.txt` | TXT | 第 N 个面片的优化准线参数（采样点） |

`{prefix}` 为 `pressure` / `suction`（简化接口为 `<文件名>_pressure` / `<文件名>_suction`）。仅导出 `.obj` 与 `.txt` 两类文件，无其他内部产物。

### 4.1 准线参数 TXT（`*_params.txt`）

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

---

## 5. 使用示例

### 5.1 C 语言

```c
#include "ruledSurfaceFitting.h"

int main() {
    RuledFitConfig cfg = {0};
    cfg.inputPath = "Blade.step";
    cfg.outputDir = "./output";
    cfg.tolerance = 0.5;

    RuledFittingResult* res = ruled_fitting(&cfg);
    if (res->errorCode == RULED_OK) {
        for (int s = 0; s < res->numSurfaces; ++s) {
            printf("%s: %d pieces, maxError=%.4f\n",
                   res->surfaces[s].name,
                   res->surfaces[s].numSegments,
                   res->surfaces[s].maxError);
        }
    }
    free_result(res);
    return 0;
}
```

### 5.2 Python ctypes 调用

```python
import ctypes

lib = ctypes.CDLL("./ruledSurfaceFitting.dll")
lib.ruled_fitting.restype = ctypes.c_void_p
lib.free_result.argtypes = [ctypes.c_void_p]

class RuledFitConfig(ctypes.Structure):
    _fields_ = [
        ("inputPath", ctypes.c_char_p),
        ("outputDir", ctypes.c_char_p),
        ("tolerance", ctypes.c_double),
    ]

cfg = RuledFitConfig()
cfg.inputPath = b"Blade.step"
cfg.outputDir = b"./output"
cfg.tolerance = 0.5

ptr = lib.ruled_fitting(ctypes.byref(cfg))
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

---

## 7. 误差计算方法

对每个面片参数域，沿准线方向取 `nAlong` 点、母线方向取 `nAcross` 点构成采样点集 {P_k}，在每个面片网格上搜索最近点：

- `maxError = max_k dist(P_k, S)` — 最大局部偏差（用于容差判断）
- `rmsError = sqrt(1/N · Σ_k dist(P_k, S)²)` — 整体拟合质量

---

## 8. 变更记录

| 版本 | 日期 | 变更 |
|------|------|------|
| v3.0 | 2026-08-28 | 单文件 API：`ruled_fitting`（3 字段 RuledFitConfig）+ `ruled_fitting_simple`。两段式流程（3 等分 → 容差判断 → 井字形细分）。仅导出 obj/txt。txt 改为导出优化后准线采样点（修复不同段 txt 相同问题）。删除旧 6 个接口。 |
| v2.0 | 2026-08-09 | 新增 4 个全自动函数（pressure/suction × ruled/plane）。 |
| v1.0 | 2026-07 | 初始版本：双文件手动模式。 |
