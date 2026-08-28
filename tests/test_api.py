"""
test_api.py — use ctypes to call ruledSurfaceFitting.dll (API v3).

Usage:
    D:\anaconda\envs\simple\python.exe tests/test_api.py [input.step] [--outdir DIR] [--tolerance T] [--gui]
"""
import sys, os, json, ctypes
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent
DLL_PATH = PROJECT_DIR / "build" / "Release" / "ruledSurfaceFitting.dll"
DEFAULT_STEP = PROJECT_DIR / "Blade.igs"
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


class RuledFitConfig(ctypes.Structure):
    _fields_ = [
        ("inputPath", ctypes.c_char_p),
        ("outputDir", ctypes.c_char_p),
        ("tolerance", ctypes.c_double),
    ]


def make_config(input_path, output_dir, tolerance):
    cfg = RuledFitConfig()
    cfg.inputPath = str(input_path).encode()
    cfg.outputDir = str(output_dir).encode()
    cfg.tolerance = tolerance
    return cfg


def run(config):
    """Call the DLL and return a copied RuledFittingResult (pointer freed)."""
    if not DLL_PATH.exists():
        raise FileNotFoundError(
            f"DLL not found: {DLL_PATH}. Build first: cmake --build build --config Release")

    os.add_dll_directory(str(DLL_PATH.parent))

    lib = ctypes.CDLL(str(DLL_PATH))
    lib.ruled_fitting.restype = ctypes.POINTER(RuledFittingResult)
    lib.ruled_fitting.argtypes = [ctypes.POINTER(RuledFitConfig)]
    lib.free_result.argtypes = [ctypes.POINTER(RuledFittingResult)]

    res_ptr = lib.ruled_fitting(ctypes.byref(config))
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
              f"(maxError={srf.maxError:.5f}, {srf.numSegments} segments):")
        for j in range(min(srf.numSegments, 64)):
            seg = srf.segments[j]
            print(f"    Seg {seg.segmentIndex}  "
                  f"maxErr={seg.maxError:.5f}  rmsErr={seg.rmsError:.5f}")

    try:
        meta = json.loads(result.metaJson.decode())
        print(f"\n  meta: {json.dumps(meta, indent=2)}")
    except Exception:
        pass


def launch_gui():
    gui_script = PROJECT_DIR / "python" / "gui_main.py"
    if not gui_script.exists():
        print(f"GUI script not found: {gui_script}")
        return
    import subprocess
    py = r"D:\anaconda\envs\simple\python.exe"
    subprocess.Popen([py, str(gui_script)], cwd=str(gui_script.parent))


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="ruledSurfaceFitting API v3 test")
    ap.add_argument("input", nargs="?", default=str(DEFAULT_STEP),
                    help="Input blade STEP/IGES file (default: Blade.igs)")
    ap.add_argument("--outdir", default=str(OUT_DIR), help="Output directory")
    ap.add_argument("--tolerance", type=float, default=0.5, help="Tolerance (mm)")
    ap.add_argument("--gui", action="store_true", help="Launch GUI after run")
    args = ap.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"[ERROR] Input file not found: {input_path}")
        sys.exit(1)

    cfg = make_config(input_path, args.outdir, args.tolerance)

    print("=== ruledSurfaceFitting API v3 Test ===")
    print(f"  DLL:     {DLL_PATH}")
    print(f"  Input:   {input_path}")
    print(f"  Output:  {args.outdir}")
    print(f"  Tolerance: {args.tolerance}")
    print()

    import shutil
    out = Path(args.outdir)
    if out.exists():
        shutil.rmtree(str(out))

    result = run(cfg)
    if result is None:
        print("DLL call returned NULL")
        sys.exit(1)

    print_summary(result)

    if args.gui:
        print("\nLaunching GUI...")
        launch_gui()
