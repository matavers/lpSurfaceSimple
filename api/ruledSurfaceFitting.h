#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32) || defined(_WIN64)
#ifdef RULED_DLL_EXPORTS
#define RULED_API __declspec(dllexport)
#else
#define RULED_API __declspec(dllimport)
#endif
#else
#define RULED_API __attribute__((visibility("default")))
#endif

/* ─── Error codes ──────────────────────────────────────────── */
typedef enum {
    RULED_OK = 0,
    RULED_ERR_FILE_NOT_FOUND = 1,
    RULED_ERR_STEP_READ_FAILED = 2,
    RULED_ERR_NO_VALID_FACE = 3,
    RULED_ERR_EXPORT_FAILED = 4,
    RULED_ERR_INVALID_PARAMS = 5
} RuledErrorCode;

/* ─── Config ───────────────────────────────────────────────── */
typedef struct {
    const char* inputPath;   /* 输入叶片 STEP/IGES 文件路径（单文件） */
    const char* outputDir;   /* 输出目录 */
    double      tolerance;   /* 容差 (mm)，<=0 时取默认 0.1 */
} RuledFitConfig;

/* ─── Per-segment result ───────────────────────────────────── */
typedef struct {
    int segmentIndex;
    double maxError;
    double rmsError;
} RuledSegmentResult;

/* ─── Per-surface result ───────────────────────────────────── */
typedef struct {
    char name[64];
    double maxError;                /* 该面整体最大误差（用于判断是否满足容差） */
    RuledSegmentResult segments[64]; /* 各面片误差（前 64 个） */
    int numSegments;                /* 实际面片数（可能超过 64，详见 metaJson） */
} RuledSurfaceResult;

/* ─── Overall result ───────────────────────────────────────── */
typedef struct {
    RuledErrorCode errorCode;
    char errorMsg[256];
    int numSurfaces;
    RuledSurfaceResult surfaces[2];
    char metaJson[4096];
} RuledFittingResult;

/* ─── API functions ────────────────────────────────────────── */

/* 主接口：单文件自动识别压力面/吸力面并拟合。
   流程：先 3 等分拟合；若最大误差 < tolerance 直接输出；
   否则按 tolerance 做井字形网格自适应细分。仅输出 .obj 与 .txt。 */
RULED_API RuledFittingResult* ruled_fitting(const RuledFitConfig* config);

/* 简化接口：批量处理输入目录下的 STEP/IGES 文件，固定 3 等分拟合（无容差判断）。 */
RULED_API RuledFittingResult* ruled_fitting_simple(const char* inputDir, const char* outputDir);

RULED_API void free_result(RuledFittingResult* result);

#ifdef __cplusplus
}
#endif
