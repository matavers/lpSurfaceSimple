/**
 * test_api.c — standalone C test for ruledSurfaceFitting.dll (API v3, 简化接口)
 *
 * Build (MSVC):
 *   cl /EHsc /I"..\api" test_api.c /link /LIBPATH:"..\build\Release" ruledSurfaceFitting.lib
 *
 * Run:
 *   test_api.exe <input_dir> <output_dir> <ruled|plane>
 */
#include <stdio.h>
#include <string.h>
#include "ruledSurfaceFitting.h"

int main(int argc, char* argv[]) {
    const char* inputDir = (argc > 1) ? argv[1] : "./input";
    const char* outputDir = (argc > 2) ? argv[2] : "./output_test_c";
    const char* mode = (argc > 3) ? argv[3] : "ruled";

    printf("=== ruledSurfaceFitting API v3 (simplified) C Test ===\n");
    printf("  Input: %s\n  Output: %s\n  Mode: %s\n", inputDir, outputDir, mode);

    RuledFittingResult* res;
    if (strcmp(mode, "plane") == 0) {
        res = plane_fitting_simple(inputDir, outputDir);
    } else {
        res = ruled_fitting_simple(inputDir, outputDir);
    }
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
    }

    printf("\n  Meta JSON:\n%s\n", res->metaJson);

    free_result(res);
    printf("\nDone. Output files in: %s\n", outputDir);
    return 0;
}
