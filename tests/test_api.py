"""
test_api.py — use ctypes to call ruledSurfaceFitting.dll (API v3, 简化接口).

Usage:
    D:\anaconda\envs\simple\python.exe tests/test_api.py [input_dir] [--outdir DIR] [--mode ruled|plane]
"""
import sys, os, json, ctypes
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent
DLL_PATH = PROJECT_DIR / "build" / "Release" / "ruledSurfaceFitting.dll"
DEFAULT_INPUT_DIR = PROJECT_DIR / "input"
OUT_DIR = PROJECT_DIR / "output"


# ── ctypes structures (must match ruledSurfaceFitting.h layout) ──

class RuledSegmentResult(ctypes.Structure):
    _fields_ = [
        ("segmentIndex", ctypes.c_int),
        ("maxError",      ctypes.c_double),
        ("rmsError",      ctypes.c_double),
    ]


class RuledSurfaceResult(ctypes.Structure):
    _fields_ = [
        ("name",        ctypes.c_char * 64),
        ("maxError",    ctypes.c_double),
        ("segments",    RuledSegmentResult * 64),
        ("numSegments", ctypes.c_int),
    ]


class RuledFittingResult(ctypes.Structure):
    _fields_ = [
        ("errorCode",   ctypes.c_int),
        ("errorMsg",    ctypes.c_char * 256),
        ("numSurfaces", ctypes.c_int),
        ("surfaces",    RuledSurfaceResult * 2),
        ("metaJson",    ctypes.c_char * 4096),
    ]


def _load_lib():
    if not DLL_PATH.exists():
        raise FileNotFoundError(
            f"DLL not found: {DLL_PATH}. Build first: cmake --build build --config Release")
    os.add_dll_directory(str(DLL_PATH.parent))
    lib = ctypes.CDLL(str(DLL_PATH))
    lib.free_result.argtypes = [ctypes.POINTER(RuledFittingResult)]
    return lib


def run_simple(fn_name, input_dir, output_dir):
    lib = _load_lib()
    fn = getattr(lib, fn_name)
    fn.restype = ctypes.POINTER(RuledFittingResult)
    fn.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
    res_ptr = fn(str(input_dir).encode(), str(output_dir).encode())
    if not res_ptr:
        return None
    result = res_ptr.contents
    lib.free_result(res_ptr)
    return result


def print_summary(result):
    if result.errorCode != 0:
        print(f"[ERROR] code={result.errorCode}  {result.errorMsg.decode()}")
        return
    print(f"[OK]  {result.numSurfaces} surfaces processed\n")
    for si in range(result.numSurfaces):
        srf = result.surfaces[si]
        print(f"  {srf.name.decode().strip()}  "
              f"(maxError={srf.maxError:.5f}, {srf.numSegments} segments)")
    try:
        meta = json.loads(result.metaJson.decode())
        print(f"\n  meta: {json.dumps(meta, indent=2)}")
    except Exception:
        pass


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="ruledSurfaceFitting API v3 (simplified) test")
    ap.add_argument("input", nargs="?", default=str(DEFAULT_INPUT_DIR),
                    help="Input directory of STEP/IGES files")
    ap.add_argument("--outdir", default=str(OUT_DIR), help="Output directory")
    ap.add_argument("--mode", default="ruled", choices=["ruled", "plane"],
                    help="ruled=直纹面 fitting, plane=平面 fitting")
    args = ap.parse_args()

    input_dir = Path(args.input)
    if not input_dir.is_dir():
        print(f"[ERROR] Input directory not found: {input_dir}")
        sys.exit(1)

    import shutil
    out = Path(args.outdir)
    if out.exists():
        shutil.rmtree(str(out))

    fn_name = "ruled_fitting_simple" if args.mode == "ruled" else "plane_fitting_simple"
    print("=== ruledSurfaceFitting API v3 (simplified) Test ===")
    print(f"  DLL:     {DLL_PATH}")
    print(f"  Input:   {input_dir}")
    print(f"  Output:  {args.outdir}")
    print(f"  Mode:    {args.mode} ({fn_name})")
    print()

    result = run_simple(fn_name, input_dir, args.outdir)
    if result is None:
        print("DLL call returned NULL")
        sys.exit(1)
    print_summary(result)
