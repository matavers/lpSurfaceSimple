"""
test_api.py — use ctypes to call simple.dll, then visualise with the GUI.

Usage:
    D:\anaconda\envs\simple\python.exe tests/test_api.py [--gui]
"""
import sys, os, json, ctypes
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent
DLL_PATH = PROJECT_DIR / "build" / "Release" / "ruledSurfaceFitting.dll"
STEP1 = PROJECT_DIR / "Blade-raw1.STEP"
STEP2 = PROJECT_DIR / "Blade-raw2.STEP"
OUT_DIR = PROJECT_DIR / "output"


# ── ctypes structures (must match simple_api.h layout) ──────────

class SimpleSegmentResult(ctypes.Structure):
    _fields_ = [
        ("segmentIndex",  ctypes.c_int),
        ("directrixDir",  ctypes.c_int),
        ("maxError",      ctypes.c_double),
        ("rmsError",      ctypes.c_double),
    ]

class SimpleSurfaceResult(ctypes.Structure):
    _fields_ = [
        ("name",           ctypes.c_char * 64),
        ("version",        ctypes.c_int),
        ("splitDir",       ctypes.c_int),
        ("numSegments",    ctypes.c_int),
        ("nCurveSamples",  ctypes.c_int),
        ("nAcrossSamples", ctypes.c_int),
        ("segments",       SimpleSegmentResult * 10),
    ]

class SimpleConfig(ctypes.Structure):
    _fields_ = [
        ("stepFile1",         ctypes.c_char_p),
        ("stepFile2",         ctypes.c_char_p),
        ("outputDir",         ctypes.c_char_p),
        ("nUSamples",         ctypes.c_int),
        ("nVSamples",         ctypes.c_int),
        ("nRibs",             ctypes.c_int),
        ("lambda",            ctypes.c_double),
        ("splitDir",          ctypes.c_int * 2),
        ("directrixDirs",     (ctypes.c_int * 10) * 2),
        ("numDirectrixDirs",  ctypes.c_int * 2),
    ]

class SimpleResult(ctypes.Structure):
    _fields_ = [
        ("errorCode",   ctypes.c_int),
        ("errorMsg",    ctypes.c_char * 256),
        ("numSurfaces", ctypes.c_int),
        ("surfaces",    SimpleSurfaceResult * 2),
        ("metaJson",    ctypes.c_char * 2048),
    ]


# ── helper to build config ──────────────────────────────────────

def make_config(surf1_split='V', surf1_dirx='V,V,V',
                surf2_split='V', surf2_dirx='V,V,V',
                nUSamples=50, nVSamples=10, nRibs=20, lambda_=1.0):
    cfg = SimpleConfig()
    cfg.stepFile1 = str(STEP1).encode()
    cfg.stepFile2 = str(STEP2).encode()
    cfg.outputDir = str(OUT_DIR).encode()
    cfg.nUSamples = nUSamples
    cfg.nVSamples = nVSamples
    cfg.nRibs = nRibs
    setattr(cfg, 'lambda', lambda_)

    def parse_dir(d):
        return {'U': 0, 'V': 1}[d.upper()]

    cfg.splitDir[0] = parse_dir(surf1_split)
    cfg.splitDir[1] = parse_dir(surf2_split)

    dirs1 = [parse_dir(d.strip()) for d in surf1_dirx.split(',')]
    dirs2 = [parse_dir(d.strip()) for d in surf2_dirx.split(',')]
    cfg.numDirectrixDirs[0] = len(dirs1)
    cfg.numDirectrixDirs[1] = len(dirs2)
    for i, d in enumerate(dirs1): cfg.directrixDirs[0][i] = d
    for i, d in enumerate(dirs2): cfg.directrixDirs[1][i] = d
    return cfg


def run(config):
    """Call the DLL and return SimpleResult."""
    if not DLL_PATH.exists():
        raise FileNotFoundError(f"DLL not found: {DLL_PATH}. Build first: cmake --build build --config Release")

    os.add_dll_directory(str(DLL_PATH.parent))

    lib = ctypes.CDLL(str(DLL_PATH))
    lib.simple_run_ruled_fitting.restype = ctypes.POINTER(SimpleResult)
    lib.simple_run_ruled_fitting.argtypes = [ctypes.POINTER(SimpleConfig)]
    lib.simple_free_result.argtypes = [ctypes.POINTER(SimpleResult)]

    res_ptr = lib.simple_run_ruled_fitting(ctypes.byref(config))
    if not res_ptr:
        return None

    result = res_ptr.contents
    return result


def print_summary(result):
    """Pretty-print the result."""
    if result.errorCode != 0:
        print(f"[ERROR] code={result.errorCode}  {result.errorMsg.decode()}")
        return

    print(f"[OK]  {result.numSurfaces} surfaces processed\n")
    for si in range(result.numSurfaces):
        srf = result.surfaces[si]
        sd = 'U' if srf.splitDir == 0 else 'V'
        print(f"  {srf.name.decode().strip()}  (split={sd}, {srf.numSegments} segments):")
        for j in range(srf.numSegments):
            seg = srf.segments[j]
            dd = 'U' if seg.directrixDir == 0 else 'V'
            print(f"    Seg {seg.segmentIndex}  (dirx={dd})  "
                  f"maxErr={seg.maxError:.5f}  rmsErr={seg.rmsError:.5f}")

    try:
        meta = json.loads(result.metaJson.decode())
        print(f"\n  meta: {json.dumps(meta, indent=2)}")
    except:
        pass


def launch_gui():
    """Launch the PyQt5+pyvistaqt GUI to visualise results."""
    gui_script = PROJECT_DIR / "python" / "gui_main.py"
    if not gui_script.exists():
        print(f"GUI script not found: {gui_script}")
        return
    import subprocess
    py = r"D:\anaconda\envs\simple\python.exe"
    subprocess.Popen([py, str(gui_script)], cwd=str(gui_script.parent))


# ── main ────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Simple API DLL test")
    ap.add_argument("--gui", action="store_true", help="Launch GUI after run")
    ap.add_argument("--split1", default="V", help="Split dir surface 1 (default: V)")
    ap.add_argument("--split2", default="V", help="Split dir surface 2 (default: V)")
    ap.add_argument("--dirx1", default="V,V,V", help="Directrix dirs surface 1 (default: V,V,V)")
    ap.add_argument("--dirx2", default="V,V,V", help="Directrix dirs surface 2 (default: V,V,V)")
    ap.add_argument("--nu", type=int, default=50, help="U-samples")
    ap.add_argument("--nv", type=int, default=10, help="V-samples")
    ap.add_argument("--ribs", type=int, default=20, help="Rib points")
    ap.add_argument("--reg", type=float, default=1.0, help="Regularization strength")
    ap.add_argument("--mode", type=str, default="ruled", help="Algorithm mode: ruled|planar")
    args = ap.parse_args()

    cfg = make_config(
        surf1_split=args.split1, surf1_dirx=args.dirx1,
        surf2_split=args.split2, surf2_dirx=args.dirx2,
        nUSamples=args.nu, nVSamples=args.nv,
        nRibs=args.ribs, lambda_=args.reg)

    print("=== Simple API DLL Test ===")
    print(f"  DLL:  {DLL_PATH}")
    print(f"  File1: {STEP1}")
    print(f"  File2: {STEP2}")
    print(f"  Output: {OUT_DIR}")
    print(f"  Split: {args.split1}/{args.split2}  Dirx: {args.dirx1} / {args.dirx2}")
    print()

    import shutil
    if OUT_DIR.exists():
        shutil.rmtree(str(OUT_DIR))

    result = run(cfg)
    if result is None:
        print("DLL call returned NULL")
        sys.exit(1)

    print_summary(result)

    if args.gui:
        print("\nLaunching GUI...")
        launch_gui()
