/**
 * test_api.c — standalone C test for ruledSurfaceFitting.dll (API v3)
 *
 * Build (MSVC):
 *   cl /EHsc /I"..\api" test_api.c /link /LIBPATH:"..\build\Release" ruledSurfaceFitting.lib
 *
 * Run:
 *   test_api.exe Blade.igs [output_dir] [tolerance]
 */
#include <stdio.h>
#include <string.h>
#include "ruledSurfaceFitting.h"

int main(int argc, char* argv[]) {
    RuledFitConfig cfg;
    memset(&cfg, 0, sizeof(cfg));

    cfg.inputPath = (argc > 1) ? argv[1] : "Blade.igs";
    cfg.outputDir = (argc > 2) ? argv[2] : "./output_test_c";
    cfg.tolerance = (argc > 3) ? atof(argv[3]) : 0.5;

    printf("=== ruledSurfaceFitting API v3 C Test ===\n");
    printf("  Input: %s\n  Output: %s\n  Tolerance: %.3f\n",
           cfg.inputPath, cfg.outputDir, cfg.tolerance);

    RuledFittingResult* res = ruled_fitting(&cfg);
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
        printf("  %s  (maxError=%.5f, %d segments):\n",
               srf->name, srf->maxError, srf->numSegments);
        int n = srf->numSegments > 64 ? 64 : srf->numSegments;
        for (int j = 0; j < n; ++j) {
            RuledSegmentResult* seg = &srf->segments[j];
            printf("    Seg %d  maxErr=%.5f  rmsErr=%.5f\n",
                   seg->segmentIndex, seg->maxError, seg->rmsError);
        }
    }

    printf("\n  Meta JSON:\n%s\n", res->metaJson);

    free_result(res);

    printf("\nDone. Output files in: %s\n", cfg.outputDir);
    return 0;
}
