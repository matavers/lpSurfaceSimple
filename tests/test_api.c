/**
 * test_api.c — standalone C test for simple.dll
 *
 * Build (MSVC):
 *   cl /EHsc /I"..\api" test_api.c /link /LIBPATH:"..\build\Release" simple.lib
 *
 * Run:
 *   test_api.exe Blade-raw1.STEP Blade-raw2.STEP
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ruledSurfaceFitting.h"

int main(int argc, char* argv[]) {
    SimpleConfig cfg;
    memset(&cfg, 0, sizeof(cfg));

    cfg.stepFile1 = (argc > 1) ? argv[1] : "Blade-raw1.STEP";
    cfg.stepFile2 = (argc > 2) ? argv[2] : "Blade-raw2.STEP";
    cfg.outputDir  = "./output_test_c";
    cfg.nUSamples  = 50;
    cfg.nVSamples  = 10;
    cfg.nRibs      = 20;

    /* Use setter-like approach to work around 'lambda' being a reserved
       word in some compilers when compiling as C++. In pure C it's fine. */
    *(double*)((char*)&cfg + offsetof(SimpleConfig, lambda)) = 1.0;

    /* Surface 1: V split, V directrix (3 segments) */
    cfg.splitDir[0] = SIMPLE_DIR_V;
    cfg.numDirectrixDirs[0] = 3;
    cfg.directrixDirs[0][0] = SIMPLE_DIR_V;
    cfg.directrixDirs[0][1] = SIMPLE_DIR_V;
    cfg.directrixDirs[0][2] = SIMPLE_DIR_V;

    /* Surface 2: V split, V directrix (3 segments) */
    cfg.splitDir[1] = SIMPLE_DIR_V;
    cfg.numDirectrixDirs[1] = 3;
    cfg.directrixDirs[1][0] = SIMPLE_DIR_V;
    cfg.directrixDirs[1][1] = SIMPLE_DIR_V;
    cfg.directrixDirs[1][2] = SIMPLE_DIR_V;

    printf("=== Simple API C Test ===\n");
    printf("  File1: %s\n  File2: %s\n  Output: %s\n",
           cfg.stepFile1, cfg.stepFile2, cfg.outputDir);

    SimpleResult* res = simple_run_ruled_fitting(&cfg);
    if (!res) {
        printf("DLL returned NULL\n");
        return 1;
    }

    if (res->errorCode != SIMPLE_OK) {
        printf("[ERROR %d] %s\n", res->errorCode, res->errorMsg);
        simple_free_result(res);
        return 1;
    }

    printf("[OK] %d surfaces processed\n\n", res->numSurfaces);
    for (int si = 0; si < res->numSurfaces; ++si) {
        SimpleSurfaceResult* srf = &res->surfaces[si];
        printf("  %s  (split=%s, %d segments):\n",
               srf->name,
               srf->splitDir == SIMPLE_DIR_U ? "U" : "V",
               srf->numSegments);
        for (int j = 0; j < srf->numSegments; ++j) {
            SimpleSegmentResult* seg = &srf->segments[j];
            printf("    Seg %d  (dirx=%s)  maxErr=%.5f  rmsErr=%.5f\n",
                   seg->segmentIndex,
                   seg->directrixDir == SIMPLE_DIR_U ? "U" : "V",
                   seg->maxError, seg->rmsError);
        }
    }

    printf("\n  Meta JSON:\n%s\n", res->metaJson);

    simple_free_result(res);

    printf("\nDone. Output files in: %s\n", cfg.outputDir);
    return 0;
}
