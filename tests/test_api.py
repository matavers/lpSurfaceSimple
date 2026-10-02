"""
test_api.py — use ctypes to call ruledSurfaceFitting.dll, then visualise with the GUI.

Usage:
    D:\anaconda\envs\simple\python.exe tests/test_api.py [--gui] [--mode ruled|planar]
"""
import sys, os, json, ctypes
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent
DLL_PATH = PROJECT_DIR / "build" / "Release" / "ruledSurfaceFitting.dll"
STEP1 = PROJECT_DIR / "Blade-raw1.STEP"
STEP2 = PROJECT_DIR / "Blade-raw2.STEP"
OUT_DIR = PROJECT_DIR / "output"


# ── ctypes structures (must match ruledSurfaceFitting.h layout) ──

class RuledCellResult(ctypes.Structure):
    _fields_ = [
        ("index",     ctypes.c_int),
        ("row",       ctypes.c_int),
        ("col",       ctypes.c_int),
        ("fitDir",    ctypes.c_int),
        ("maxError",  ctypes.c_double),
        ("rmsError",  ctypes.c_double),
    ]


class RuledSurfaceResult(ctypes.Structure):
    _fields_ = [
        ("name",     ctypes.c_char * 64),
        ("cells",    RuledCellResult * 512),
        ("numCells", ctypes.c_int),
    ]


class RuledConfig(ctypes.Structure):
    _fields_ = [
        ("stepFile1", ctypes.c_char_p),
        ("stepFile2", ctypes.c_char_p),
        ("outputDir", ctypes.c_char_p),
        ("nUSamples", ctypes.c_int),
        ("nVSamples", ctypes.c_int),
        ("nRibs",     ctypes.c_int),
        ("lambda",    ctypes.c_double),
        ("nSplitU",   ctypes.c_int),
        ("nSplitV",   ctypes.c_int),
        ("tolerance", ctypes.c_double),
        ("maxDepth",  ctypes.c_int),
        ("faceIdx",   ctypes.c_int * 2),
    ]


class PlanarConfig(ctypes.Structure):
    _fields_ = [
        ("stepFile1", ctypes.c_char_p),
        ("stepFile2", ctypes.c_char_p),
        ("outputDir", ctypes.c_char_p),
        ("nUSamples", ctypes.c_int),
        ("nVSamples", ctypes.c_int),
        ("nSplitU",   ctypes.c_int),
        ("nSplitV",   ctypes.c_int),
        ("tolerance", ctypes.c_double),
        ("maxDepth",  ctypes.c_int),
        ("faceIdx",   ctypes.c_int * 2),
    ]


class RuledFittingResult(ctypes.Structure):
    _fields_ = [
        ("errorCode",   ctypes.c_int),
        ("errorMsg",    ctypes.c_char * 256),
        ("numSurfaces", ctypes.c_int),
        ("surfaces",    RuledSurfaceResult * 2),
        ("metaJson",    ctypes.c_char * 4096),
    ]


def make_ruled_config(nU=50, nV=10, nRibs=20, reg=1.0,
                      nSplitU=2, nSplitV=2, tolerance=0.1, maxDepth=20):
    cfg = RuledConfig()
    cfg.stepFile1 = str(STEP1).encode()
    cfg.stepFile2 = str(STEP2).encode()
    cfg.outputDir = str(OUT_DIR).encode()
    cfg.nUSamples = nU
    cfg.nVSamples = nV
    cfg.nRibs = nRibs
    setattr(cfg, 'lambda', reg)
    cfg.nSplitU = nSplitU
    cfg.nSplitV = nSplitV
    cfg.tolerance = tolerance
    cfg.maxDepth = maxDepth
    cfg.faceIdx[0] = -1
    cfg.faceIdx[1] = -1
    return cfg


def make_planar_config(nU=50, nV=10, nSplitU=2, nSplitV=2,
                       tolerance=0.3, maxDepth=20):
    cfg = PlanarConfig()
    cfg.stepFile1 = str(STEP1).encode()
    cfg.stepFile2 = str(STEP2).encode()
    cfg.outputDir = str(OUT_DIR).encode()
    cfg.nUSamples = nU
    cfg.nVSamples = nV
    cfg.nSplitU = nSplitU
    cfg.nSplitV = nSplitV
    cfg.tolerance = tolerance
    cfg.maxDepth = maxDepth
    cfg.faceIdx[0] = -1
    cfg.faceIdx[1] = -1
    return cfg


def run(config, func_name):
    if not DLL_PATH.exists():
        raise FileNotFoundError(f"DLL not found: {DLL_PATH}. Build first: cmake --build build --config Release")

    os.add_dll_directory(str(DLL_PATH.parent))
    lib = ctypes.CDLL(str(DLL_PATH))

    fn = getattr(lib, func_name)
    fn.restype = ctypes.POINTER(RuledFittingResult)
    fn.argtypes = [ctypes.POINTER(type(config))]
    lib.free_result.argtypes = [ctypes.POINTER(RuledFittingResult)]

    res_ptr = fn(ctypes.byref(config))
    if not res_ptr:
        return None
    return res_ptr.contents


def print_summary(result):
    if result.errorCode != 0:
        print(f"[ERROR] code={result.errorCode}  {result.errorMsg.decode()}")
        return

    print(f"[OK]  {result.numSurfaces} surfaces processed\n")
    for si in range(result.numSurfaces):
        srf = result.surfaces[si]
        print(f"  {srf.name.decode().strip()}  ({srf.numCells} cells):")
        for j in range(srf.numCells):
            c = srf.cells[j]
            d = 'U' if c.fitDir == 0 else 'V'
            print(f"    Cell {c.index} (r{c.row},c{c.col}) [dir={d}] "
                  f"maxErr={c.maxError:.5f}  rmsErr={c.rmsError:.5f}")

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
    ap = argparse.ArgumentParser(description="Simple API DLL test")
    ap.add_argument("--gui", action="store_true", help="Launch GUI after run")
    ap.add_argument("--mode", type=str, default="ruled", help="ruled|planar")
    ap.add_argument("--nu", type=int, default=50)
    ap.add_argument("--nv", type=int, default=10)
    ap.add_argument("--ribs", type=int, default=20)
    ap.add_argument("--reg", type=float, default=1.0)
    ap.add_argument("--nsplit-u", type=int, default=2)
    ap.add_argument("--nsplit-v", type=int, default=2)
    ap.add_argument("--tolerance", type=float, default=0.1)
    ap.add_argument("--max-depth", type=int, default=20)
    args = ap.parse_args()

    import shutil
    if OUT_DIR.exists():
        shutil.rmtree(str(OUT_DIR))

    print("=== Simple API DLL Test ===")
    print(f"  DLL:    {DLL_PATH}")
    print(f"  Mode:   {args.mode}")
    print(f"  Output: {OUT_DIR}")

    if args.mode == "planar":
        cfg = make_planar_config(nU=args.nu, nV=args.nv,
                                 nSplitU=args.nsplit_u, nSplitV=args.nsplit_v,
                                 tolerance=args.tolerance, maxDepth=args.max_depth)
        func = "plane_surface_fitting"
    else:
        cfg = make_ruled_config(nU=args.nu, nV=args.nv, nRibs=args.ribs, reg=args.reg,
                                nSplitU=args.nsplit_u, nSplitV=args.nsplit_v,
                                tolerance=args.tolerance, maxDepth=args.max_depth)
        func = "ruled_surface_fitting"

    result = run(cfg, func)
    if result is None:
        print("DLL call returned NULL")
        sys.exit(1)

    print_summary(result)

    if args.gui:
        print("\nLaunching GUI...")
        launch_gui()
