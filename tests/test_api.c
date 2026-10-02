/**
 * test_api.c — standalone C test for ruledSurfaceFitting.dll
 *
 * Build (MSVC):
 *   cl /EHsc /I"..\api" test_api.c /link /LIBPATH:"..\build\Release" ruledSurfaceFitting.lib
 *
 * Run:
 *   test_api.exe Blade-raw1.STEP Blade-raw2.STEP
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ruledSurfaceFitting.h"

int main(int argc, char* argv[]) {
    RuledConfig cfg;
    memset(&cfg, 0, sizeof(cfg));

    cfg.stepFile1 = (argc > 1) ? argv[1] : "Blade-raw1.STEP";
    cfg.stepFile2 = (argc > 2) ? argv[2] : "Blade-raw2.STEP";
    cfg.outputDir  = "./output_test_c";
    cfg.nUSamples  = 50;
    cfg.nVSamples  = 10;
    cfg.nRibs      = 20;
    cfg.lambda     = 1.0;
    cfg.nSplitU    = 2;
    cfg.nSplitV    = 2;
    cfg.tolerance  = 0.1;
    cfg.maxDepth   = 20;
    cfg.faceIdx[0] = -1;
    cfg.faceIdx[1] = -1;

    printf("=== Simple API C Test ===\n");
    printf("  File1: %s\n  File2: %s\n  Output: %s\n",
           cfg.stepFile1, cfg.stepFile2, cfg.outputDir);

    RuledFittingResult* res = ruled_surface_fitting(&cfg);
    if (!res) {
        printf("DLL returned NULL\n");
        return 1;
    }

    if (res->errorCode != RULED_OK) {
        printf("[ERROR %d] %s\n", res->errorCode, res->errorMsg);
        free_result(res);
        return 1;
    }

    printf("[OK] %d surfaces processed\n\n", res->numSurfaces);
    for (int si = 0; si < res->numSurfaces; ++si) {
        RuledSurfaceResult* srf = &res->surfaces[si];
        printf("  %s  (%d cells):\n", srf->name, srf->numCells);
        for (int j = 0; j < srf->numCells; ++j) {
            RuledCellResult* c = &srf->cells[j];
            printf("    Cell %d (r%d,c%d) [dir=%s]  maxErr=%.5f  rmsErr=%.5f\n",
                   c->index, c->row, c->col,
                   c->fitDir == RULED_DIR_U ? "U" : "V",
                   c->maxError, c->rmsError);
        }
    }

    printf("\n  Meta JSON:\n%s\n", res->metaJson);

    free_result(res);

    printf("\nDone. Output files in: %s\n", cfg.outputDir);
    return 0;
}
