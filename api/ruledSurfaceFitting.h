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

/* ─── Parameter / fit direction (per cell, auto-selected) ──── */
typedef enum {
    RULED_DIR_U = 0,
    RULED_DIR_V = 1
} RuledDirection;

/* ─── Config (ruled) ───────────────────────────────────────── */
typedef struct {
    const char* stepFile1;
    const char* stepFile2;
    const char* outputDir;
    int nUSamples;          /* default 50 */
    int nVSamples;          /* default 10 */
    int nRibs;              /* default 20 */
    double lambda;          /* default 1.0 */
    int nSplitU;            /* U-direction (vertical) splits -> columns   default 2 */
    int nSplitV;            /* V-direction (horizontal) splits -> rows    default 2 */
    double tolerance;       /* default 0.1 */
    int maxDepth;           /* default 20 */
    int faceIdx[2];         /* -1 = auto-pick largest */
} RuledConfig;

/* ─── Config (planar) ──────────────────────────────────────── */
typedef struct {
    const char* stepFile1;
    const char* stepFile2;
    const char* outputDir;
    int nUSamples;          /* default 50 */
    int nVSamples;          /* default 10 */
    int nSplitU;            /* default 2 */
    int nSplitV;            /* default 2 */
    double tolerance;       /* default 0.1 */
    int maxDepth;           /* default 20 */
    int faceIdx[2];         /* -1 = auto-pick largest */
} PlanarConfig;

/* ─── Per-cell result ──────────────────────────────────────── */
typedef struct {
    int index;
    int row;
    int col;
    int fitDir;             /* RuledDirection: optimizer-selected (U or V); planar = U */
    double maxError;
    double rmsError;
    double maxOvercut;      /* 最大过切量 mm（signed<0 的幅度，拟合面切入材料侧） */
    double maxUndercut;     /* 最大欠切量 mm（signed>0 的幅度，拟合面留料） */
    double meanSigned;      /* 平均符号误差 mm（>0 欠切，<0 过切） */
} RuledCellResult;

/* ─── Per-surface result ───────────────────────────────────── */
typedef struct {
    char name[64];
    RuledCellResult cells[512];
    int numCells;
    double maxError;        /* 该面最大误差 mm */
    double rmsError;        /* 该面整体 RMS 误差 mm */
} RuledSurfaceResult;

/* ─── Overall result ───────────────────────────────────────── */
typedef struct {
    RuledErrorCode errorCode;
    char errorMsg[256];
    int numSurfaces;
    RuledSurfaceResult surfaces[2];
    char metaJson[4096];
    double fitSec;          /* 分区拟合耗时 s */
} RuledFittingResult;

/* ─── Toolpath (machining) result ─────────────────────────── */
typedef struct {
    RuledErrorCode errorCode;
    char errorMsg[256];
    int numPatches;
    int flankRegions;
    double flankCut;      /* s */
    double flankOverhead; /* s */
    double flankTotal;    /* s */
    double pointCut;      /* s */
    double pointTotal;    /* s */
    double speedup;       /* pointTotal / flankTotal */
    double totalArea;     /* mm^2 */
    double originalArea;  /* mm^2 */
    double elapsedSec;    /* 刀轨计算耗时 s */
} ToolpathResult;

/* ─── API functions ────────────────────────────────────────── */

/* Two-file input, manual face selection */
RULED_API RuledFittingResult* ruled_surface_fitting(const RuledConfig* config);
RULED_API RuledFittingResult* plane_surface_fitting(const PlanarConfig* config);

/* Single-file auto blade processing (auto-identify + split + grid fit) */
RULED_API RuledFittingResult* pressure_ruled_fitting(const RuledConfig* config);
RULED_API RuledFittingResult* pressure_plane_fitting(const PlanarConfig* config);
RULED_API RuledFittingResult* suction_ruled_fitting(const RuledConfig* config);
RULED_API RuledFittingResult* suction_plane_fitting(const PlanarConfig* config);

/* Toolpath (machining): read exported *_params.txt + meta.json + blade*_mesh.obj
   from inputPath, compute flank/point toolpaths + machining time comparison,
   write toolpath_*.vtk/csv/dxf + summary.json to outputPath. */
RULED_API ToolpathResult* compute_toolpath(const char* inputPath,
                                           const char* outputPath);

RULED_API void free_result(RuledFittingResult* result);
RULED_API void free_toolpath_result(ToolpathResult* result);

#ifdef __cplusplus
}
#endif
