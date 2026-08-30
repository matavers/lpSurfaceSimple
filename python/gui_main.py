"""
Simple - Ruled Surface Fitting GUI
Simplified from distillation gui_main.py

Tree structure:
  Surface 1
    ├── Original Mesh
    └── Version 0
          ├── Ruled Seg 0
          ├── Ruled Seg 1
          └── Ruled Seg 2
  Surface 2
    ├── Original Mesh
    └── Version 0
          ├── Ruled Seg 0
          ├── Ruled Seg 1
          └── Ruled Seg 2
"""
import sys, os, json, re, subprocess, tempfile
from pathlib import Path
from datetime import datetime

import numpy as np

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QComboBox, QCheckBox, QFileDialog, QTextEdit,
    QTreeWidget, QTreeWidgetItem, QSplitter, QFormLayout,
    QSpinBox, QDoubleSpinBox, QMessageBox, QStackedWidget, QDialog,
)
from PyQt5.QtCore import Qt, QTimer, pyqtSignal, QThread
from PyQt5.QtGui import QFont, QTextCursor, QPixmap

try:
    from pyvistaqt import QtInteractor
    HAS_PYVISTA = True
except ImportError:
    HAS_PYVISTA = False

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
BUILD_EXE = PROJECT_DIR / "build" / "Release" / "simple.exe"
CONFIG_PATH = SCRIPT_DIR / ".simple_config.json"

sys.path.insert(0, str(PROJECT_DIR))
try:
    from research import compute_machining as _cm
    from research import sweep as _sweep
    HAS_MACHINING = True
except Exception:
    HAS_MACHINING = False

BLADE_COLORS = [
    [0.30, 0.60, 0.95],
    [0.95, 0.50, 0.30],
]

SEG_COLORS = [
    [0.20, 0.85, 0.40],
    [0.85, 0.70, 0.15],
    [0.85, 0.30, 0.50],
]

PREVIEW_COLORS = [
    [0.88, 0.72, 0.53],
    [0.62, 0.78, 0.88],
    [0.82, 0.62, 0.72],
    [0.68, 0.83, 0.63],
    [0.78, 0.68, 0.83],
    [0.88, 0.83, 0.58],
]

HIGHLIGHT_COLOR = [1.0, 0.84, 0.0]

CONSOLE_CSS = """
QTextEdit {
    background: #1e1e1e; color: #d4d4d4;
    font-family: Consolas; font-size: 13px; padding: 4px;
}
"""


class ProcRunner(QThread):
    output_signal = pyqtSignal(str)
    finished_signal = pyqtSignal(int)

    def __init__(self, cmd, cwd=None):
        super().__init__()
        self.cmd = cmd
        self.cwd = cwd

    def run(self):
        try:
            self._proc = subprocess.Popen(
                self.cmd, shell=True, cwd=self.cwd,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1
            )
            for line in iter(self._proc.stdout.readline, ''):
                self.output_signal.emit(line.rstrip())
            self._proc.wait()
            self.finished_signal.emit(self._proc.returncode)
        except Exception as e:
            self.output_signal.emit(f"[Error] {e}")
            self.finished_signal.emit(1)


class SweepWorker(QThread):
    done = pyqtSignal(object)
    failed = pyqtSignal(str)

    def __init__(self, blade, dll, tolerances, params):
        super().__init__()
        self._blade = blade
        self._dll = dll
        self._tolerances = tolerances
        self._params = params
        self._stop = False

    def run(self):
        try:
            rows = []
            for tol in self._tolerances:
                if self._stop:
                    break
                out_dir = os.path.join(tempfile.gettempdir(), f"sweep_gui_{tol:.3f}")
                if os.path.isdir(out_dir):
                    shutil.rmtree(out_dir)
                code = _sweep.run_fitting(self._dll, self._blade, out_dir, tol)
                if code != 0:
                    continue
                patches = _cm.compute(out_dir, self._params)
                if not patches:
                    continue
                s = _cm.summarize(patches, self._params)
                rows.append({
                    "tolerance": tol,
                    "num_patches": s["num_patches"],
                    "developable": s["developable"],
                    "non_developable": s["non_developable"],
                    "max_twist": max((p.twist for p in patches), default=0.0),
                    "flank_total": s["flank"]["total"],
                    "point_total": s["point"]["total"],
                    "speedup": s["point"]["total"] / s["flank"]["total"] if s["flank"]["total"] > 0 else 0.0,
                })
            self.done.emit(rows)
        except Exception as e:
            self.failed.emit(str(e))

    def request_stop(self):
        self._stop = True


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Simple - Ruled Surface Fitting")
        self.resize(1400, 900)

        self._file1 = str(PROJECT_DIR / "Blade-raw1.STEP")
        self._file2 = str(PROJECT_DIR / "Blade-raw2.STEP")
        self._out_dir = str(PROJECT_DIR / "output")
        self._nSamplesU = 50
        self._nSamplesV = 10
        self._numSegments = 3
        self._tolerance = 0.1
        self._surfaceSegCounts = [3, 3]
        self._currentVersion = 0
        self._mode = "ruled"
        self._single_file_mode = False
        self._uRange1 = None
        self._uRange2 = None
        self._mesh_shared = False
        self._preview_colors = {}
        self._preview_face_ids = set()

        self._proc = None
        self._loaded_files = set()

        self._setup_menu()
        self._setup_ui()
        self._load_config()

    def _setup_menu(self):
        m = self.menuBar().addMenu("File")
        m.addAction("Import Existing Output...", self._on_import_output)
        m.addAction("Clear Results", self._on_clear_results)
        m.addSeparator()
        m.addAction("Exit", self.close)

        wb = self.menuBar().addMenu("工作台")
        self._act_wb_fit = wb.addAction("拟合", lambda: self._switch_workbench(0))
        self._act_wb_fit.setCheckable(True)
        self._act_wb_mach = wb.addAction("加工仿真", lambda: self._switch_workbench(1))
        self._act_wb_mach.setCheckable(True)
        self._act_wb_fit.setChecked(True)

    def _switch_workbench(self, idx):
        self._workbench = idx
        self._act_wb_fit.setChecked(idx == 0)
        self._act_wb_mach.setChecked(idx == 1)
        self._wb_stack.setCurrentIndex(idx)

    def _on_clear_results(self):
        self._tree.clear()
        self._loaded_files.clear()
        self._clear_3d()
        self._console.clear()
        self._uRange1 = None
        self._uRange2 = None
        self._mesh_shared = False
        self._log("Results cleared.")

    def _on_import_output(self):
        d = QFileDialog.getExistingDirectory(
            self, "Select Output Directory", self._out_dir)
        if not d:
            return
        self._out_dir = d
        self._txt_outdir.setText(d)
        self._tree.clear()
        self._loaded_files.clear()
        self._clear_3d()
        self._console.clear()

        self._log(f"Importing data from: {d}")
        meta = os.path.join(d, "meta.json")
        if os.path.exists(meta):
            self._on_file_written(os.path.normpath(meta))
        else:
            self._log("[Warning] No meta.json found in selected directory.")

    def _setup_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        ml = QVBoxLayout(central)
        ml.setContentsMargins(4, 4, 4, 4)

        splitter = QSplitter(Qt.Horizontal)

        if HAS_PYVISTA:
            self._plotter = QtInteractor(self, shape=(1, 1))
            self._plotter.set_background('lightblue')
            splitter.addWidget(self._plotter)
        else:
            splitter.addWidget(QLabel("pyvistaqt not installed. Run: pip install pyvistaqt"))

        right = QWidget()
        rl = QVBoxLayout(right)
        rl.setContentsMargins(4, 0, 4, 0)

        self._wb_stack = QStackedWidget()

        fit_page = QWidget()
        fl = QVBoxLayout(fit_page)
        fl.setContentsMargins(0, 0, 0, 0)
        self._setup_params(fl)
        self._setup_tree(fl)
        self._wb_stack.addWidget(fit_page)

        mach_page = QWidget()
        mach_layout = QVBoxLayout(mach_page)
        mach_layout.setContentsMargins(0, 0, 0, 0)
        self._setup_machining(mach_layout)
        self._wb_stack.addWidget(mach_page)

        rl.addWidget(self._wb_stack)
        splitter.addWidget(right)
        splitter.setSizes([900, 500])
        ml.addWidget(splitter, 1)

        self._console = QTextEdit()
        self._console.setReadOnly(True)
        self._console.setMaximumHeight(220)
        self._console.setStyleSheet(CONSOLE_CSS)
        self._console.setFont(QFont("Consolas", 12))
        ml.addWidget(self._console)

    def _setup_params(self, layout):
        form = QFormLayout()
        form.setSpacing(6)

        self._chk_single = QCheckBox("Single File Mode (multi-face import)")
        self._chk_single.toggled.connect(self._on_toggle_single)
        form.addRow(self._chk_single)

        row1 = QHBoxLayout()
        self._txt_file1 = QLineEdit(self._file1)
        row1.addWidget(self._txt_file1)
        btn1 = QPushButton("...")
        btn1.setMaximumWidth(30)
        btn1.clicked.connect(lambda: self._on_browse_file(1))
        row1.addWidget(btn1)
        form.addRow("Blade 1 STEP:", row1)

        row2 = QHBoxLayout()
        self._txt_file2 = QLineEdit(self._file2)
        row2.addWidget(self._txt_file2)
        btn2 = QPushButton("...")
        btn2.setMaximumWidth(30)
        btn2.clicked.connect(lambda: self._on_browse_file(2))
        row2.addWidget(btn2)
        form.addRow("Blade 2 STEP:", row2)

        row_pv = QHBoxLayout()
        self._btn_preview = QPushButton("Load Preview")
        self._btn_preview.clicked.connect(self._on_load_preview)
        row_pv.addWidget(self._btn_preview)
        self._btn_auto_id = QPushButton("Auto-Identify")
        self._btn_auto_id.clicked.connect(self._on_auto_identify)
        row_pv.addWidget(self._btn_auto_id)
        self._btn_split = QPushButton("Split Blade")
        self._btn_split.clicked.connect(self._on_split_blade)
        row_pv.addWidget(self._btn_split)
        row_pv.addStretch()
        form.addRow(row_pv)

        self._cmb_face1 = QComboBox()
        self._cmb_face1.addItem("(select face)")
        self._cmb_face1.currentIndexChanged.connect(self._on_face_selection_changed)
        form.addRow("Face for Surface 1:", self._cmb_face1)

        self._cmb_face2 = QComboBox()
        self._cmb_face2.addItem("(select face)")
        self._cmb_face2.currentIndexChanged.connect(self._on_face_selection_changed)
        form.addRow("Face for Surface 2:", self._cmb_face2)

        self._hide_single_controls()

        from PyQt5.QtWidgets import QListWidget, QWidget
        self._split_panel = QWidget()
        sp_layout = QVBoxLayout(self._split_panel)
        sp_layout.setContentsMargins(0, 0, 0, 0)
        split_label = QLabel("Split Results")
        sp_layout.addWidget(split_label)
        self._split_list = QListWidget()
        self._split_list.setMaximumHeight(120)
        self._split_list.itemDoubleClicked.connect(self._on_split_item_double_clicked)
        sp_layout.addWidget(self._split_list)
        id_bar = QHBoxLayout()
        self._btn_identify = QPushButton("Identify Pressure/Suction")
        self._btn_identify.clicked.connect(self._on_identify)
        self._btn_identify.setEnabled(False)
        id_bar.addWidget(self._btn_identify)
        id_bar.addStretch()
        sp_layout.addLayout(id_bar)
        form.addRow(self._split_panel)
        self._split_panel.setVisible(False)

        row3 = QHBoxLayout()
        self._txt_outdir = QLineEdit(self._out_dir)
        row3.addWidget(self._txt_outdir)
        btn3 = QPushButton("...")
        btn3.setMaximumWidth(30)
        btn3.clicked.connect(self._on_browse_outdir)
        row3.addWidget(btn3)
        form.addRow("Output Dir:", row3)

        self._spn_u = QSpinBox()
        self._spn_u.setRange(10, 200)
        self._spn_u.setValue(self._nSamplesU)
        form.addRow("U Samples:", self._spn_u)

        self._spn_v = QSpinBox()
        self._spn_v.setRange(5, 50)
        self._spn_v.setValue(self._nSamplesV)
        form.addRow("V Samples:", self._spn_v)

        self._spn_seg = QSpinBox()
        self._spn_seg.setRange(1, 10)
        self._spn_seg.setValue(self._numSegments)
        form.addRow("Segments:", self._spn_seg)

        self._spn_tol = QDoubleSpinBox()
        self._spn_tol.setRange(0.001, 10.0)
        self._spn_tol.setDecimals(3)
        self._spn_tol.setSingleStep(0.01)
        self._spn_tol.setValue(self._tolerance)
        form.addRow("Tolerance (mm):", self._spn_tol)

        self._cmb_split1 = QComboBox()
        self._cmb_split1.addItems(["U", "V"])
        form.addRow("Split Dir (Surface 1):", self._cmb_split1)

        self._cmb_split2 = QComboBox()
        self._cmb_split2.addItems(["U", "V"])
        form.addRow("Split Dir (Surface 2):", self._cmb_split2)

        self._cmb_dirx1 = []
        form.addRow(QLabel("Directrix Dir (Surface 1):"))
        for seg in range(3):
            cmb = QComboBox()
            cmb.addItems(["V", "U"])
            row = QHBoxLayout()
            row.addWidget(QLabel(f"  Seg {seg+1}:"))
            row.addWidget(cmb)
            row.addStretch()
            form.addRow(row)
            self._cmb_dirx1.append(cmb)

        self._cmb_dirx2 = []
        form.addRow(QLabel("Directrix Dir (Surface 2):"))
        for seg in range(3):
            cmb = QComboBox()
            cmb.addItems(["V", "U"])
            row = QHBoxLayout()
            row.addWidget(QLabel(f"  Seg {seg+1}:"))
            row.addWidget(cmb)
            row.addStretch()
            form.addRow(row)
            self._cmb_dirx2.append(cmb)

        bar = QHBoxLayout()
        self._btn_run = QPushButton("Run")
        self._btn_run.setStyleSheet(
            "QPushButton{background:#4CAF50;color:white;font-size:18px;padding:8px}")
        self._btn_run.clicked.connect(self._on_run)
        bar.addWidget(self._btn_run, 2)

        self._cmb_mode = QComboBox()
        self._cmb_mode.addItems(["Ruled (Fixed)", "Ruled (Adaptive)", "Planar (Fixed)", "Planar (Adaptive)"])
        self._cmb_mode.currentIndexChanged.connect(self._on_mode_changed)
        bar.addWidget(self._cmb_mode, 1)

        self._btn_stop = QPushButton("Stop")
        self._btn_stop.setEnabled(False)
        self._btn_stop.setStyleSheet(
            "QPushButton{background:#888;color:#ddd;font-size:18px;padding:8px}")
        self._btn_stop.clicked.connect(self._stop)
        bar.addWidget(self._btn_stop)
        form.addRow(bar)

        layout.addLayout(form)

    def _setup_tree(self, layout):
        layout.addWidget(QLabel("Visualization"))
        self._tree = QTreeWidget()
        self._tree.setHeaderHidden(True)
        self._tree.itemChanged.connect(self._on_check)
        layout.addWidget(self._tree)

    # ── 加工仿真工作台 ──────────────────────────────────────────
    def _setup_machining(self, layout):
        if not HAS_MACHINING:
            layout.addWidget(QLabel("未找到 research/compute_machining.py"))
            return

        mcfg = _cm.load_config()
        form = QFormLayout()
        form.setSpacing(6)

        self._spn_m_feed = QDoubleSpinBox()
        self._spn_m_feed.setRange(1.0, 100000.0)
        self._spn_m_feed.setValue(mcfg.get("feed", 500.0))
        form.addRow("进给率 (mm/min):", self._spn_m_feed)

        self._spn_m_ballr = QDoubleSpinBox()
        self._spn_m_ballr.setRange(0.1, 100.0)
        self._spn_m_ballr.setValue(mcfg.get("ball_r", 5.0))
        form.addRow("球头刀半径 (mm):", self._spn_m_ballr)

        self._spn_m_scallop = QDoubleSpinBox()
        self._spn_m_scallop.setRange(0.001, 10.0)
        self._spn_m_scallop.setDecimals(3)
        self._spn_m_scallop.setValue(mcfg.get("scallop", 0.01))
        form.addRow("残留高度 (mm):", self._spn_m_scallop)

        self._spn_m_twist = QDoubleSpinBox()
        self._spn_m_twist.setRange(0.01, 90.0)
        self._spn_m_twist.setDecimals(2)
        self._spn_m_twist.setValue(mcfg.get("twist_limit", 1.0))
        form.addRow("可展阈值 (度):", self._spn_m_twist)

        self._spn_m_overhead = QDoubleSpinBox()
        self._spn_m_overhead.setRange(0.0, 100.0)
        self._spn_m_overhead.setValue(mcfg.get("overhead", 2.0))
        form.addRow("进退刀开销 (s/块):", self._spn_m_overhead)

        self._spn_m_poverhead = QDoubleSpinBox()
        self._spn_m_poverhead.setRange(0.0, 100.0)
        self._spn_m_poverhead.setValue(mcfg.get("point_overhead", 10.0))
        form.addRow("点铣进退刀 (s):", self._spn_m_poverhead)

        layout.addLayout(form)

        row1 = QHBoxLayout()
        self._btn_tool = QPushButton("刀具计算")
        self._btn_tool.clicked.connect(self._on_compute_tool)
        row1.addWidget(self._btn_tool)
        self._btn_tool_stop = QPushButton("急停")
        self._btn_tool_stop.setEnabled(False)
        self._btn_tool_stop.clicked.connect(self._on_stop_tool)
        row1.addWidget(self._btn_tool_stop)
        layout.addLayout(row1)

        row2 = QHBoxLayout()
        self._btn_clear_tool = QPushButton("清除刀轨结果")
        self._btn_clear_tool.clicked.connect(self._on_clear_toolpath)
        row2.addWidget(self._btn_clear_tool)
        self._btn_switch_viz = QPushButton("切换可视化内容")
        self._btn_switch_viz.clicked.connect(self._on_switch_viz)
        row2.addWidget(self._btn_switch_viz)
        layout.addLayout(row2)

        layout.addWidget(QLabel("仿真结果"))
        self._txt_mach = QTextEdit()
        self._txt_mach.setReadOnly(True)
        self._txt_mach.setMaximumHeight(240)
        layout.addWidget(self._txt_mach)

        self._viz_mode = 0          # 0=刀轨, 1=扫描图表
        self._tool_actors = []
        self._mach_patches = None
        self._mach_summary = None
        self._sweep_worker = None

    def _mach_args(self):
        import argparse
        return argparse.Namespace(
            feed=self._spn_m_feed.value(), ball_r=self._spn_m_ballr.value(),
            scallop=self._spn_m_scallop.value(), twist_limit=self._spn_m_twist.value(),
            overhead=self._spn_m_overhead.value(),
            point_overhead=self._spn_m_poverhead.value())

    def _on_compute_tool(self):
        if not HAS_MACHINING or not HAS_PYVISTA:
            QMessageBox.warning(self, "Error", "加工仿真模块或 pyvista 不可用。")
            return
        out_dir = self._out_dir
        if not os.path.isdir(out_dir) or not list(Path(out_dir).glob("*_params.txt")):
            QMessageBox.warning(self, "Error", "请先在拟合工作台运行拟合，生成 *_params.txt。")
            return
        args = self._mach_args()
        patches = _cm.compute(out_dir, args)
        if not patches:
            QMessageBox.warning(self, "Error", "未能解析准线参数文件。")
            return
        self._mach_patches = patches
        self._mach_summary = _cm.summarize(patches, args)
        self._show_mach_report(patches, self._mach_summary)
        self._viz_mode = 0
        self._render_toolpath()

    def _show_mach_report(self, patches, summary):
        lines = []
        f, p = summary['flank'], summary['point']
        lines.append(f"面片数: {summary['num_patches']}（可展 {summary['developable']}，不可展 {summary['non_developable']}）")
        lines.append(f"侧铣: 切削 {f['cut']:.1f}s + 非切削 {f['overhead']:.1f}s = {f['total']:.1f}s")
        lines.append(f"点铣: 切削 {p['cut']:.1f}s + 非切削 {p['overhead']:.1f}s = {p['total']:.1f}s")
        if f['total'] > 0:
            lines.append(f"侧铣相对点铣提速: {p['total'] / f['total']:.1f}x")
        lines.append("")
        for pp in patches:
            lines.append(f"{pp.name}: 扭转 {pp.twist:.2f}° "
                         f"{'可展' if pp.developable else '不可展'} "
                         f"侧铣 {pp.flank_time:.1f}s / 点铣 {pp.point_time:.1f}s")
        self._txt_mach.setPlainText("\n".join(lines))
        self._log(f"[Machining] {summary['num_patches']} patches, "
                  f"flank {f['total']:.1f}s vs point {p['total']:.1f}s")

    def _render_toolpath(self):
        self._on_clear_toolpath()
        if not self._mach_patches or not HAS_PYVISTA:
            return
        import numpy as np
        args = self._mach_args()
        stepover = 2.0 * math.sqrt(max(0.0, 2.0 * args.ball_r * args.scallop - args.scallop ** 2))
        for pp in self._mach_patches:
            for line in _cm.flank_toolpath_lines(pp.C0, pp.C1):
                pts = np.array(line)
                if len(pts) >= 2:
                    a = self._plotter.add_lines(pts, color='#d62728', width=3,
                                                name=f"flank_{pp.name}_{len(self._tool_actors)}")
                    self._tool_actors.append(a)
            for line in _cm.point_toolpath_lines(pp.C0, pp.C1, stepover):
                pts = np.array(line)
                if len(pts) >= 2:
                    a = self._plotter.add_lines(pts, color='#1f77b4', width=1,
                                                name=f"point_{pp.name}_{len(self._tool_actors)}")
                    self._tool_actors.append(a)
        self._plotter.render()

    def _on_clear_toolpath(self):
        if not HAS_PYVISTA:
            return
        for a in self._tool_actors:
            try:
                self._plotter.remove_actor(a)
            except Exception:
                pass
        self._tool_actors = []
        self._plotter.render()

    def _on_stop_tool(self):
        if self._sweep_worker and self._sweep_worker.isRunning():
            self._sweep_worker.request_stop()
            self._log("[Machining] 请求停止扫描...")

    def _on_switch_viz(self):
        if not HAS_MACHINING:
            return
        if self._viz_mode == 0:
            self._viz_mode = 1
            self._run_sweep()
        else:
            self._viz_mode = 0
            self._render_toolpath()

    def _run_sweep(self):
        blade = self._file1
        if not os.path.isfile(blade):
            QMessageBox.warning(self, "Error", "请先在拟合工作台选择叶片文件。")
            self._viz_mode = 0
            return
        tolerances = _cm.load_config().get("tolerances", [0.05, 0.1, 0.2, 0.5, 1.0])
        dll = str(PROJECT_DIR / "build" / "Release" / "ruledSurfaceFitting.dll")
        self._btn_tool_stop.setEnabled(True)
        self._log("[Machining] 开始容差扫描...")
        self._sweep_worker = SweepWorker(blade, dll, tolerances, self._mach_args())
        self._sweep_worker.done.connect(self._on_sweep_done)
        self._sweep_worker.failed.connect(lambda e: self._log(f"[Machining] 扫描失败: {e}"))
        self._sweep_worker.finished.connect(lambda: self._btn_tool_stop.setEnabled(False))
        self._sweep_worker.start()

    def _on_sweep_done(self, rows):
        if not rows:
            self._log("[Machining] 扫描无结果")
            return
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        tols = [r['tolerance'] for r in rows]
        flank = [r['flank_total'] for r in rows]
        point = [r['point_total'] for r in rows]
        n_patches = [r['num_patches'] for r in rows]

        fig, ax1 = plt.subplots(figsize=(6, 4), dpi=110)
        ax1.plot(tols, flank, 'o-', color='#d62728', label='侧铣')
        ax1.plot(tols, point, 's-', color='#1f77b4', label='点铣')
        ax1.set_xlabel('容差 (mm)')
        ax1.set_ylabel('加工时间 (s)')
        ax1.legend(loc='upper left')
        ax2 = ax1.twinx()
        ax2.plot(tols, n_patches, 'x--', color='#2ca02c', label='面片数')
        ax2.set_ylabel('面片数')
        ax2.legend(loc='upper right')
        fig.tight_layout()

        png = os.path.join(tempfile.gettempdir(), 'mach_sweep.png')
        fig.savefig(png)
        plt.close(fig)

        dlg = QDialog(self)
        dlg.setWindowTitle("参数扫描结果")
        dlg.resize(680, 480)
        vl = QVBoxLayout(dlg)
        lbl = QLabel()
        lbl.setPixmap(QPixmap(png))
        lbl.setScaledContents(True)
        vl.addWidget(lbl)
        dlg.exec_()

        self._log("[Machining] 扫描完成，图表已显示。")

    def _on_browse_file(self, idx):
        path, _ = QFileDialog.getOpenFileName(
            self, f"Select Blade {idx} STEP File",
            str(PROJECT_DIR), "CAD Files (*.step *.stp *.igs *.iges *.STEP *.STP *.IGS *.IGES);;STEP (*.step *.stp);;IGES (*.igs *.iges);;All (*)")
        if path:
            if idx == 1:
                self._file1 = path
                self._txt_file1.setText(path)
            else:
                self._file2 = path
                self._txt_file2.setText(path)

    def _on_browse_outdir(self):
        d = QFileDialog.getExistingDirectory(
            self, "Select Output Directory", self._out_dir)
        if d:
            self._out_dir = d
            self._txt_outdir.setText(d)

    def _get_dirx_str(self, surface_idx):
        cmbs = self._cmb_dirx1 if surface_idx == 1 else self._cmb_dirx2
        return ",".join(c.currentText().lower() for c in cmbs)

    def _log(self, text):
        ts = datetime.now().strftime("[%H:%M:%S] ")
        self._console.append(ts + text)
        c = self._console.textCursor()
        c.movePosition(QTextCursor.End)
        self._console.setTextCursor(c)

    def _load_config(self):
        try:
            with open(CONFIG_PATH) as f:
                cfg = json.load(f)
            self._file1 = cfg.get("file1", self._file1)
            self._file2 = cfg.get("file2", self._file2)
            self._out_dir = cfg.get("out_dir", self._out_dir)
            self._nSamplesU = cfg.get("nUSamples", 50)
            self._nSamplesV = cfg.get("nVSamples", 10)
            self._numSegments = cfg.get("numSegments", 3)
            self._txt_file1.setText(self._file1)
            self._txt_file2.setText(self._file2)
            self._txt_outdir.setText(self._out_dir)
            self._spn_u.setValue(self._nSamplesU)
            self._spn_v.setValue(self._nSamplesV)
            self._spn_seg.setValue(self._numSegments)
        except Exception:
            pass

    def _save_config(self):
        cfg = {
            "file1": self._txt_file1.text(),
            "file2": self._txt_file2.text(),
            "out_dir": self._txt_outdir.text(),
            "nUSamples": self._spn_u.value(),
            "nVSamples": self._spn_v.value(),
            "numSegments": self._spn_seg.value(),
        }
        try:
            with open(CONFIG_PATH, 'w') as f:
                json.dump(cfg, f, indent=2)
        except Exception:
            pass

    def _on_mode_changed(self, idx):
        if idx == 0:
            self._mode = "ruled"
        elif idx == 1:
            self._mode = "ruled-adaptive"
        elif idx == 2:
            self._mode = "planar"
        else:
            self._mode = "planar-adaptive"

    def _hide_single_controls(self):
        for w in [self._btn_preview, self._btn_auto_id, self._btn_split,
                  self._cmb_face1, self._cmb_face2]:
            w.setVisible(False)
        if hasattr(self, '_split_panel'):
            self._split_panel.setVisible(False)

    def _lock_face_controls(self):
        self._cmb_face1.setEnabled(False)
        self._cmb_face2.setEnabled(False)
        self._btn_split.setEnabled(False)

    def _unlock_face_controls(self):
        self._cmb_face1.setEnabled(True)
        self._cmb_face2.setEnabled(True)
        self._btn_split.setEnabled(True)

    def _on_toggle_single(self, checked):
        self._single_file_mode = checked
        if checked:
            self._txt_file2.setVisible(False)
            self._btn_preview.setVisible(True)
            self._btn_auto_id.setVisible(True)
            self._btn_split.setVisible(True)
            self._cmb_face1.setVisible(True)
            self._cmb_face2.setVisible(True)
            if hasattr(self, '_split_panel'):
                self._split_panel.setVisible(True)
            self._cmb_face1.clear(); self._cmb_face1.addItem("(click Load Preview first)")
            self._cmb_face2.clear(); self._cmb_face2.addItem("(click Load Preview first)")
        else:
            self._txt_file2.setVisible(True)
            self._hide_single_controls()

    def _on_load_preview(self):
        exe = str(BUILD_EXE)
        if not os.path.exists(exe):
            self._log(f"[Error] simple.exe not found: {exe}")
            return
        filepath = self._txt_file1.text()
        if not os.path.exists(filepath):
            self._log(f"[Error] File not found: {filepath}")
            return

        import subprocess, json
        self._clear_3d()
        self._log(f"Loading preview: {filepath}")
        try:
            result = subprocess.run(
                [exe, filepath, '--mode', 'info'],
                cwd=str(PROJECT_DIR), capture_output=True, text=True,
                encoding='utf-8', errors='replace', timeout=30)
            raw = result.stdout.strip()
            idx = raw.find('{')
            if idx >= 0:
                raw = raw[idx:]
            info = json.loads(raw)
            self._log(f"  Found {info['numFaces']} faces")

            self._cmb_face1.clear()
            self._cmb_face2.clear()
            if not HAS_PYVISTA:
                return

            self._labels_pts = {}
            labels_pts = []
            labels_text = []

            for face in info['faces']:
                idx = face['index']
                u0, u1 = face.get('uMin', 0), face.get('uMax', 0)
                v0, v1 = face.get('vMin', 0), face.get('vMax', 0)
                diag = face.get('diag', 0)
                label = f"F{idx}  diag={diag:.1f}mm  UV=[{u0:.2f},{u1:.2f}]x[{v0:.2f},{v1:.2f}]"
                self._cmb_face1.addItem(label, idx)
                self._cmb_face2.addItem(label, idx)

                tmp_path = os.path.join(tempfile.gettempdir(), f"simple_preview_f{idx}.obj")
                try:
                    subprocess.run(
                        [exe, filepath, '--mode', 'face-obj',
                         '--face-idx', str(idx), '--face-out', tmp_path],
                        cwd=str(PROJECT_DIR), capture_output=True, text=True,
                        encoding='utf-8', errors='replace', timeout=20)
                except subprocess.TimeoutExpired:
                    self._log(f"    Face {idx}: mesh generation timed out, skipping")
                    continue
                except Exception as e:
                    self._log(f"    Face {idx}: mesh generation failed ({e})")
                    continue

                if os.path.exists(tmp_path) and os.path.getsize(tmp_path) > 100:
                    import pyvista as pv
                    m = pv.read(tmp_path)
                    ci = idx % len(PREVIEW_COLORS)
                    color = PREVIEW_COLORS[ci]
                    self._preview_colors[idx] = color
                    self._preview_face_ids.add(idx)
                    self._plotter.add_mesh(
                        m, name=f"preview_face_{idx}",
                        color=color, opacity=0.92,
                        smooth_shading=True,
                        show_edges=False,
                        pbr=True, metallic=0.05, roughness=0.4)
                    ctr = m.center
                    self._labels_pts[idx] = [ctr[0], ctr[1], ctr[2]]
                    labels_pts.append([ctr[0], ctr[1], ctr[2]])
                    labels_text.append(f"F{idx}")
                    self._log(f"    Face {idx}: diag={diag:.1f}mm ({m.n_points} pts)")
                else:
                    self._log(f"    Face {idx}: empty mesh, skipping")

            if labels_pts:
                self._plotter.add_point_labels(
                    np.array(labels_pts), labels_text,
                    font_size=14, bold=True,
                    text_color='black',
                    shape='rounded_rect', shape_color='white',
                    shape_opacity=0.85, margin=4,
                    point_size=1, always_visible=True,
                    name='face_labels')
                self._log(f"  Added {len(labels_pts)} face labels")

            self._plotter.view_isometric()
            self._plotter.render()
            self._setup_mouse_picking()
            self._log("  Preview loaded. Click faces or use dropdowns to select.")

        except Exception as e:
            self._log(f"[Error] Preview failed: {e}")

    def _on_face_selection_changed(self):
        if not HAS_PYVISTA or not self._preview_colors:
            return
        sel = set()
        d1 = self._cmb_face1.currentData()
        d2 = self._cmb_face2.currentData()
        if d1 is not None:
            sel.add(d1 if not isinstance(d1, dict) else d1.get('face_idx', -1))
        if d2 is not None:
            sel.add(d2 if not isinstance(d2, dict) else d2.get('face_idx', -1))
        try:
            for fid in self._preview_face_ids:
                name = f"preview_face_{fid}"
                actor = None
                try:
                    actor = self._plotter.actors.get(name)
                except Exception:
                    for a in self._plotter.renderer._actors:
                        if hasattr(a, '_name') and a._name == name:
                            actor = a
                            break
                if actor is None:
                    continue
                if fid in sel:
                    actor.GetProperty().SetColor(HIGHLIGHT_COLOR)
                    actor.GetProperty().SetOpacity(1.0)
                else:
                    orig = self._preview_colors.get(fid, [0.7, 0.7, 0.7])
                    actor.GetProperty().SetColor(orig)
                    actor.GetProperty().SetOpacity(0.92)
            self._plotter.render()
        except Exception:
            pass

    def _setup_mouse_picking(self):
        if not HAS_PYVISTA or not self._preview_face_ids:
            return
        try:
            import vtk
        except ImportError:
            return
        self._pick_slot = 0
        self._picker = vtk.vtkPropPicker()

        def on_click(obj, evt):
            if not self._single_file_mode:
                return
            x, y = obj.GetEventPosition()
            self._picker.Pick(x, y, 0, self._plotter.renderer)
            actor = self._picker.GetActor()
            if actor:
                name = getattr(actor, '_name', '')
                if name.startswith('preview_face_'):
                    fid = int(name.rsplit('_', 1)[-1])
                    self._on_face_picked(fid)
        self._pick_observer = self._plotter.iren.interactor.AddObserver(
            "LeftButtonPressEvent", on_click)

    def _on_face_picked(self, fid):
        self._last_picked_fid = fid
        slot = self._pick_slot
        self._pick_slot = 1 - self._pick_slot
        cmb = self._cmb_face1 if slot == 0 else self._cmb_face2
        for i in range(cmb.count()):
            if cmb.itemData(i) == fid:
                cmb.setCurrentIndex(i)
                return

    def _on_split_blade(self):
        exe = str(BUILD_EXE)
        if not os.path.exists(exe):
            self._log(f"[Error] simple.exe not found: {exe}")
            return
        filepath = self._txt_file1.text()
        if not os.path.exists(filepath):
            self._log(f"[Error] File not found: {filepath}")
            return

        faceIdx = self._last_picked_fid
        if faceIdx is None:
            self._log("Click a face in 3D view first, then Split Blade")
            return

        self._log(f"Splitting Face[{faceIdx}]...")
        try:
            sr = subprocess.run(
                [exe, filepath, '--mode', 'split-blade', '--face-idx', str(faceIdx)],
                cwd=str(PROJECT_DIR), capture_output=True, text=True,
                encoding='utf-8', errors='replace', timeout=30)
            sraw = sr.stdout.strip()
            s2 = sraw.find('{')
            if s2 >= 0: sraw = sraw[s2:]
            sinfo = json.loads(sraw)
            if sinfo.get('success'):
                regions = sinfo['regions']
                splitDir = sinfo.get('dir', 'U')
                self._log(f"  Got {len(regions)} regions from Face[{faceIdx}] dir={splitDir}")
                self._hide_actor(f"preview_face_{faceIdx}")

                self._split_list.clear()
                self._split_items = []
                self._split_face_idx = faceIdx
                bp_set = set()
                for reg in regions:
                    us = reg['uStart']; ue = reg['uEnd']
                    if us > 0.001: bp_set.add(us)
                    if ue < 0.999: bp_set.add(ue)
                bp_list = sorted(bp_set)

                for k, reg in enumerate(regions):
                    us = reg['uStart']
                    ue = reg['uEnd']
                    label = reg.get('label', '')
                    name = f"split_{faceIdx}_{k}"
                    sub_path = os.path.join(tempfile.gettempdir(), f"simple_split_{faceIdx}_{k}.obj")
                    if splitDir == 'V':
                        args = ['--u-range1', '0.25,0.75', '--v-range1', f'{us},{ue}']
                    else:
                        args = ['--u-range1', f'{us},{ue}', '--v-range1', '0.25,0.75']
                    subprocess.run(
                        [exe, filepath, '--mode', 'face-obj-range',
                         '--face-idx', str(faceIdx),
                         *args,
                         '--face-out', sub_path],
                        cwd=str(PROJECT_DIR), capture_output=True, text=True,
                        encoding='utf-8', errors='replace', timeout=15)
                    if os.path.exists(sub_path) and os.path.getsize(sub_path) > 100:
                        import pyvista as pv
                        m = pv.read(sub_path)
                        ci = k % len(PREVIEW_COLORS)
                        color = PREVIEW_COLORS[ci]
                        self._plotter.add_mesh(
                            m, name=name, color=color, opacity=0.92,
                            smooth_shading=True, show_edges=False,
                            pbr=True, metallic=0.05, roughness=0.4)
                        self._split_items.append({
                            "face_idx": faceIdx, "sub_idx": k,
                            "u_start": us, "u_end": ue,
                            "label": label, "name": name,
                            "dir": splitDir
                        })
                        item_text = f"Face{faceIdx}-split-{k}  {label}"
                        self._split_list.addItem(item_text)
                        self._log(f"    [{k}] U=[{us:.3f},{ue:.3f}] {label}")

                self._plotter.render()
                self._btn_identify.setEnabled(len(self._split_items) >= 2)
                self._refocus_camera()
                self._log(f"  Split complete. {len(self._split_items)} regions in list.")

                guide_path = os.path.join(tempfile.gettempdir(),
                    f"simple_guide_{faceIdx}.obj")
                subprocess.run(
                    [exe, filepath, '--mode', 'iso-curve',
                     '--face-idx', str(faceIdx),
                     '--v-range1', '0.5,0.5',
                     '--face-out', guide_path],
                    cwd=str(PROJECT_DIR), capture_output=True, text=True,
                    encoding='utf-8', errors='replace', timeout=10)
                if os.path.exists(guide_path) and os.path.getsize(guide_path) > 100:
                    import pyvista as pv
                    gl = pv.read(guide_path)
                    self._plotter.add_mesh(
                        gl, name=f"guide_{faceIdx}",
                        color='blue', line_width=5, opacity=0.8)
                    self._log("  Blue line = span direction (root-to-tip guide)")
                self._plotter.render()
            else:
                self._log(f"  Split failed: {sinfo.get('message','')}")
        except Exception as e:
            self._log(f"[Error] Split failed: {e}")

    def _on_split_item_double_clicked(self, item):
        idx = self._split_list.row(item)
        if 0 <= idx < len(self._split_items):
            name = self._split_items[idx].get('name', '')
            try:
                for fid in self._preview_face_ids:
                    pname = f"preview_face_{fid}"
                    try:
                        a = self._plotter.actors.get(pname)
                    except Exception:
                        a = None
                    if a:
                        a.SetVisibility(False)
                for a in self._plotter.renderer._actors:
                    if hasattr(a, '_name') and a._name == name:
                        a.GetProperty().SetColor(HIGHLIGHT_COLOR)
                        a.GetProperty().SetOpacity(1.0)
                    elif hasattr(a, '_name') and a._name.startswith('split_'):
                        orig = self._preview_colors.get(fid, [0.7, 0.7, 0.7])
                        a.GetProperty().SetColor(orig)
                        a.GetProperty().SetOpacity(0.92)
                self._plotter.render()
            except Exception:
                pass

    def _on_identify(self):
        if not self._split_items:
            return
        self._log("Identifying pressure/suction from split regions...")
        curv_data = [(i, si['label'], si.get('u_start', 0), si.get('u_end', 1))
                     for i, si in enumerate(self._split_items)]
        max_c = float('-inf')
        best_i = -1
        for c in curv_data:
            if c[2] + c[3] > max_c:
                max_c = c[2] + c[3]
        for i, label, us, ue in curv_data:
            if label == 'edge':
                self._hide_actor(self._split_items[i]['name'])
            elif label == 'suction':
                self._split_items[i]['label'] = 'suction'
            elif label == 'pressure':
                self._split_items[i]['label'] = 'pressure'
        self._cmb_face1.clear()
        self._cmb_face2.clear()
        kept = []
        for si in self._split_items:
            if si['label'] in ('suction', 'pressure'):
                kept.append(si)
        if len(kept) >= 2:
            self._cmb_face1.addItem(f"Face{kept[0]['face_idx']}-split-{kept[0]['sub_idx']} {kept[0]['label']}",
                                    kept[0])
            self._cmb_face2.addItem(f"Face{kept[1]['face_idx']}-split-{kept[1]['sub_idx']} {kept[1]['label']}",
                                    kept[1])
            self._log(f"  Surface 1: {self._cmb_face1.currentText()}")
            self._log(f"  Surface 2: {self._cmb_face2.currentText()}")
        elif len(kept) == 1:
            self._cmb_face1.addItem(f"Face{kept[0]['face_idx']}-split-{kept[0]['sub_idx']} {kept[0]['label']}",
                                    kept[0])
            self._log(f"  Only one side found: {self._cmb_face1.currentText()}")
        self._refocus_camera()
        self._log(f"  Ready. Click Run to process.")

    def _hide_actor(self, name):
        try:
            a = self._plotter.actors.get(name)
        except Exception:
            a = None
        if a:
            a.SetVisibility(False)
            return
        for a in self._plotter.renderer._actors:
            if hasattr(a, '_name') and a._name == name:
                a.SetVisibility(False)
                return

    def _clear_split_results(self):
        for si in self._split_items:
            self._hide_actor(si['name'])
        self._split_items = []
        self._split_list.clear()
        self._btn_identify.setEnabled(False)

    def _update_labels(self):
        if not HAS_PYVISTA or not hasattr(self, '_labels_pts'):
            return
        try:
            self._plotter.remove_actor('face_labels')
        except Exception:
            pass
        pts = []
        txt = []
        for fid in self._preview_face_ids:
            try:
                a = self._plotter.actors.get(f"preview_face_{fid}")
            except Exception:
                a = None
            if a is None:
                continue
            if not a.GetVisibility():
                continue
            if fid in self._labels_pts:
                pts.append(self._labels_pts[fid])
                txt.append(f"F{fid}")
        if pts:
            self._plotter.add_point_labels(
                np.array(pts), txt,
                font_size=14, bold=True,
                text_color='black',
                shape='rounded_rect', shape_color='white',
                shape_opacity=0.85, margin=4,
                point_size=1, always_visible=True,
                name='face_labels')

    def _refocus_camera(self):
        if not HAS_PYVISTA:
            return
        import numpy as np
        all_pts = []
        try:
            for name, a in self._plotter.actors.items():
                if not a.GetVisibility():
                    continue
                if hasattr(a, 'GetMapper') and a.GetMapper():
                    inp = a.GetMapper().GetInput()
                    if inp and hasattr(inp, 'GetPoints') and inp.GetPoints():
                        pts = inp.GetPoints()
                        for k in range(pts.GetNumberOfPoints()):
                            all_pts.append(pts.GetPoint(k))
        except Exception:
            pass
        if all_pts:
            arr = np.array(all_pts)
            ctr = arr.mean(axis=0)
            self._plotter.camera.SetFocalPoint(ctr[0], ctr[1], ctr[2])
            self._plotter.render()

    def _on_auto_identify(self):
        exe = str(BUILD_EXE)
        if not os.path.exists(exe):
            self._log(f"[Error] simple.exe not found: {exe}")
            return
        filepath = self._txt_file1.text()
        if not os.path.exists(filepath):
            self._log(f"[Error] File not found: {filepath}")
            return

        self._log(f"Auto-identifying blade faces: {filepath}")
        try:
            result = subprocess.run(
                [exe, filepath, '--mode', 'auto-identify'],
                cwd=str(PROJECT_DIR), capture_output=True, text=True,
                encoding='utf-8', errors='replace', timeout=30)
            raw = result.stdout.strip()
            idx = raw.find('{')
            if idx >= 0:
                raw = raw[idx:]
            info = json.loads(raw)
            if info.get('success'):
                pi = info['pressureIndex']
                si = info['suctionIndex']
                pd = info.get('pressureDiag', 0)
                sd = info.get('suctionDiag', 0)

                self._uRange1 = None
                self._uRange2 = None

                self._cmb_face1.clear()
                self._cmb_face2.clear()
                self._cmb_face1.addItem(f"Face {pi} diag={pd:.1f}mm (pressure)", pi)
                self._cmb_face2.addItem(f"Face {si} diag={sd:.1f}mm (suction)", si)
                self._cmb_face1.setCurrentIndex(0)
                self._cmb_face2.setCurrentIndex(0)

                self._log(f"  Identified {info.get('message','')}")
                self._log(f"  Pressure face[{pi}] diag={pd:.1f}mm")
                self._log(f"  Suction  face[{si}] diag={sd:.1f}mm")
                self._log(f"  Select faces above, optionally Split Blade, then Run.")

                blade_set = {pi, si}
                for fid in self._preview_face_ids:
                    if fid not in blade_set:
                        self._hide_actor(f"preview_face_{fid}")
                self._update_labels()
                self._refocus_camera()
                self._last_picked_fid = pi if pd > sd else si
                self._log(f"  Non-blade faces hidden. Click a blade face, then Split.")
            else:
                self._log(f"  Identification failed: {info.get('message','')}")
                QMessageBox.warning(self, "Auto-Identify", info.get('message', 'Failed'))
        except Exception as e:
            self._log(f"[Error] Auto-identify failed: {e}")

    def _on_run(self):
        if self._proc and self._proc.isRunning():
            QMessageBox.information(self, "Running", "Already running.")
            return

        self._file1 = self._txt_file1.text()
        if self._single_file_mode:
            self._file2 = ""  # will be auto-resolved from sidecar
        else:
            self._file2 = self._txt_file2.text()
        self._out_dir = self._txt_outdir.text()
        self._nSamplesU = self._spn_u.value()
        self._nSamplesV = self._spn_v.value()
        self._numSegments = self._spn_seg.value()
        self._tolerance = self._spn_tol.value()

        if not os.path.exists(self._file1):
            QMessageBox.warning(self, "Error", "File not found.")
            return
        if self._file2 and not os.path.exists(self._file2):
            QMessageBox.warning(self, "Error", "File 2 not found.")
            return

        self._save_config()

        import shutil
        if os.path.isdir(self._out_dir):
            try:
                shutil.rmtree(self._out_dir)
            except Exception:
                pass
        os.makedirs(self._out_dir, exist_ok=True)

        self._tree.clear()
        self._loaded_files.clear()
        self._surfaceSegCounts = [3, 3]
        self._clear_3d()
        self._console.clear()

        self._log(f"File 1: {self._file1}")
        self._log(f"File 2: {self._file2}")
        self._log(f"Output: {self._out_dir}")

        exe = str(BUILD_EXE)
        if not os.path.exists(exe):
            self._log(f"[Error] simple.exe not found: {exe}")
            QMessageBox.warning(self, "Error", f"simple.exe not found at:\n{exe}")
            return

        cmd = (
            f'"{exe}" '
            f'"{self._file1}" "{self._file2}" '
            f'--mode {self._mode} '
            f'--outdir "{self._out_dir}" '
            f'--nusamples {self._nSamplesU} '
            f'--nvsamples {self._nSamplesV} '
            f'--numsegments {self._numSegments} '
            f'--tolerance {self._tolerance} '
            f'--split-dir1 {self._cmb_split1.currentText().lower()} '
            f'--split-dir2 {self._cmb_split2.currentText().lower()} '
            f'--dirx-dir1 {self._get_dirx_str(1)} '
            f'--dirx-dir2 {self._get_dirx_str(2)}'
        )
        if self._single_file_mode:
            f1 = self._cmb_face1.currentData()
            f2 = self._cmb_face2.currentData()
            if f1 is not None:
                if isinstance(f1, dict):
                    cmd += f' --face-idx1 {f1["face_idx"]}'
                    rkey1 = '--v-range1' if f1.get('dir') == 'V' else '--u-range1'
                    cmd += f' {rkey1} {f1["u_start"]:.6f},{f1["u_end"]:.6f}'
                else:
                    cmd += f' --face-idx1 {f1}'
                    if self._uRange1:
                        cmd += f' --u-range1 {self._uRange1[0]:.6f},{self._uRange1[1]:.6f}'
            if f2 is not None:
                if isinstance(f2, dict):
                    cmd += f' --face-idx2 {f2["face_idx"]}'
                    rkey2 = '--v-range2' if f2.get('dir') == 'V' else '--u-range2'
                    cmd += f' {rkey2} {f2["u_start"]:.6f},{f2["u_end"]:.6f}'
                else:
                    cmd += f' --face-idx2 {f2}'
                    if self._uRange2:
                        cmd += f' --u-range2 {self._uRange2[0]:.6f},{self._uRange2[1]:.6f}'
        self._log(f"Running: {cmd}")

        self._btn_run.setEnabled(False)
        self._btn_stop.setEnabled(True)
        self._btn_stop.setStyleSheet(
            "QPushButton{background:#c33;color:white;font-size:18px;padding:8px}")
        self._lock_face_controls()

        self._proc = ProcRunner(cmd, str(PROJECT_DIR))
        self._proc.output_signal.connect(self._on_output)
        self._proc.finished_signal.connect(self._on_done)
        self._proc.start()

    def _on_output(self, line):
        self._log(line)
        try:
            if 'wrote ' in line.lower():
                m = re.search(r'wrote\s+(\S+)', line)
                if m:
                    self._on_file_written(m.group(1))
        except Exception as e:
            self._log(f"[GUI Error] _on_output: {e}")

    def _on_file_written(self, path):
        try:
            path = os.path.normpath(path)
            if path in self._loaded_files:
                return
            self._loaded_files.add(path)
            fn = os.path.basename(path)
            self._log(f"[GUI] Detected file: {fn}")

            if fn == 'meta.json':
                self._load_meta(path)
                self._build_tree()
                self._apply_visibility()
        except Exception as e:
            self._log(f"[GUI Error] _on_file_written: {e}")

    def _build_tree(self):
        self._tree.clear()
        self._log(f"[Tree] surfaceSegCounts={self._surfaceSegCounts}")

        seg_label = "Plane Seg" if self._mode.startswith("planar") else "Ruled Seg"
        seg_prefix = "plane" if self._mode.startswith("planar") else "seg"

        for bi in range(2):
            blade_name = f"Surface {bi + 1}"
            bnode = QTreeWidgetItem([blade_name])
            bnode.setFlags(bnode.flags() | Qt.ItemIsUserCheckable)
            bnode.setCheckState(0, Qt.Checked)
            bnode.setExpanded(True)
            self._tree.addTopLevelItem(bnode)

            mesh_item = QTreeWidgetItem(["Original Mesh"])
            mesh_item.setFlags(mesh_item.flags() | Qt.ItemIsUserCheckable)
            mesh_item.setCheckState(0, Qt.Unchecked)
            mesh_item.setData(1, Qt.UserRole, f"blade{bi + 1}")
            mesh_item.setData(2, Qt.UserRole, "mesh")
            bnode.addChild(mesh_item)

            ver_node = QTreeWidgetItem([f"Version {self._currentVersion}"])
            ver_node.setFlags(ver_node.flags() | Qt.ItemIsUserCheckable)
            ver_node.setCheckState(0, Qt.Checked)
            ver_node.setExpanded(True)
            bnode.addChild(ver_node)

            for seg in range(self._surfaceSegCounts[bi]):
                seg_item = QTreeWidgetItem([f"{seg_label} {seg + 1}"])
                seg_item.setFlags(seg_item.flags() | Qt.ItemIsUserCheckable)
                seg_item.setCheckState(0, Qt.Checked)
                seg_item.setData(1, Qt.UserRole, f"blade{bi + 1}_{seg_prefix}{seg}")
                seg_item.setData(2, Qt.UserRole, "ruled")
                ver_node.addChild(seg_item)

    def _load_meta(self, meta_path):
        try:
            with open(meta_path) as f:
                raw = f.read()
            self._log(f"[Meta] raw_len={len(raw)}")
            meta = json.load(raw)
            self._log(f"[Meta] surfaces={len(meta.get('surfaces', []))}")
            if "mode" in meta:
                self._mode = meta["mode"]
                if self._mode == "ruled":
                    self._cmb_mode.setCurrentIndex(0)
                elif self._mode == "ruled-adaptive":
                    self._cmb_mode.setCurrentIndex(1)
                elif self._mode == "planar":
                    self._cmb_mode.setCurrentIndex(2)
                else:
                    self._cmb_mode.setCurrentIndex(3)
                self._log(f"  Mode: {self._mode}")
            for i, s in enumerate(meta.get("surfaces", [])):
                ns = len(s.get('segments', []))
                if i < 2:
                    self._surfaceSegCounts[i] = ns
                if ns > self._numSegments:
                    self._numSegments = ns
                self._log(f"  {s['name']}: {ns} segments")
        except Exception as e:
            self._log(f"[Meta Error] {e}")

        self._load_all_objs()

    def _load_all_objs(self):
        out_dir = self._out_dir
        self._mesh_shared = False
        if HAS_PYVISTA:
            self._plotter.disable_render = True

        ruled_files = []
        mesh_files = []
        for fn in sorted(os.listdir(out_dir)):
            if not fn.endswith('.obj'):
                continue
            if fn.endswith('_mesh.obj'):
                mesh_files.append(fn)
            else:
                ruled_files.append(fn)

        # 根据实际输出文件数量更新面片数（不依赖 meta.json 的 segments 字段，避免可视化面片数与实际不符）
        cnt = [0, 0]
        for fn in ruled_files:
            if fn.startswith('blade1_'):
                cnt[0] += 1
            elif fn.startswith('blade2_'):
                cnt[1] += 1
        if cnt[0] > 0:
            self._surfaceSegCounts[0] = cnt[0]
        if cnt[1] > 0:
            self._surfaceSegCounts[1] = cnt[1]

        for fn in ruled_files:
            path = os.path.normpath(os.path.join(out_dir, fn))
            name = fn.replace('.obj', '')
            self._load_obj(path, name, "ruled")

        self._log(f"[Load] ruled_files={len(ruled_files)} mesh_files={len(mesh_files)}")

        mesh_paths = [os.path.normpath(os.path.join(out_dir, fn)) for fn in mesh_files]
        if len(mesh_paths) == 2:
            try:
                s1 = os.path.getsize(mesh_paths[0])
                s2 = os.path.getsize(mesh_paths[1])
                if s1 == s2 and s1 > 0:
                    self._mesh_shared = True
            except Exception:
                pass

        for fn in mesh_files:
            path = os.path.normpath(os.path.join(out_dir, fn))
            name = fn.replace('.obj', '')
            if self._mesh_shared and 'blade2' in name:
                continue
            self._load_obj(path, name, "mesh")

        if HAS_PYVISTA:
            self._plotter.disable_render = False
            self._apply_visibility()
            self._plotter.render()

    def _load_obj(self, path, name, tag):
        if not HAS_PYVISTA:
            return
        try:
            import pyvista as pv
            m = pv.read(path)

            if tag == "mesh":
                if "blade1" in name:
                    color = BLADE_COLORS[0]
                else:
                    color = BLADE_COLORS[1]
                self._plotter.add_mesh(
                    m, name=name, color=color, opacity=0.85,
                    show_edges=True, edge_color='darkgray')
                self._log(f"  Loaded mesh: {name}")
            elif tag == "ruled":
                sidx = int(re.search(r'(?:seg|plane)(\d+)', name).group(1))
                total = self._surfaceSegCounts[0] if "blade1" in name else self._surfaceSegCounts[1]
                import colorsys
                hue = (sidx / max(1, total)) % 1.0
                r, g, b = colorsys.hsv_to_rgb(hue, 0.65, 0.95)
                color = [r, g, b]
                opacity = 0.80 if "plane" in name else 0.35
                self._plotter.add_mesh(
                    m, name=name, color=color, opacity=opacity,
                    show_edges=True, edge_color='dimgray')
                self._log(f"  Loaded: {name}")

        except Exception as e:
            self._log(f"  Load error {name}: {e}")

    def _on_check(self, item, col):
        self._apply_visibility()

    def _apply_visibility(self):
        if not HAS_PYVISTA:
            return

        visible_set = set()

        def walk(node, inherited_vis):
            vis = inherited_vis and node.checkState(0) == Qt.Checked
            data_name = node.data(1, Qt.UserRole)
            tag = node.data(2, Qt.UserRole)
            if vis and data_name:
                if tag == "mesh":
                    blade = int(data_name.replace('blade', ''))
                    name = f"blade{blade}_mesh"
                    visible_set.add(name)
                    if self._mesh_shared and blade == 2:
                        visible_set.add("blade1_mesh")
                elif tag == "ruled":
                    visible_set.add(data_name)
            for i in range(node.childCount()):
                walk(node.child(i), vis)

        for i in range(self._tree.topLevelItemCount()):
            top = self._tree.topLevelItem(i)
            walk(top, top.checkState(0) == Qt.Checked)

        if self._mesh_shared:
            for tag in list(visible_set):
                if tag == "blade2_mesh":
                    visible_set.discard("blade2_mesh")
                    visible_set.add("blade1_mesh")

        try:
            actor_count = 0
            for a in self._plotter.renderer._actors:
                if hasattr(a, '_name'):
                    actor_count += 1
                    a.SetVisibility(a._name in visible_set)
            self._log(f"[Vis] visible_set={len(visible_set)} actors={actor_count}")
            self._plotter.render()
        except Exception:
            pass

    def _clear_3d(self):
        self._preview_colors = {}
        self._preview_face_ids = set()
        self._last_picked_fid = None
        self._split_items = []
        self._split_list.clear()
        self._btn_identify.setEnabled(False)
        if HAS_PYVISTA:
            self._plotter.clear()
            self._plotter.set_background('lightblue')

    def _on_done(self, code):
        self._btn_run.setEnabled(True)
        self._btn_stop.setEnabled(False)
        self._btn_stop.setStyleSheet(
            "QPushButton{background:#888;color:#ddd;font-size:18px;padding:8px}")
        self._unlock_face_controls()
        if code == 0:
            self._log("[Done] Algorithm finished successfully.")
        else:
            self._log(f"[Error] Algorithm failed with code {code}.")

    def _stop(self):
        if self._proc:
            try:
                pid = self._proc._proc.pid
                subprocess.run(
                    f'taskkill /F /T /PID {pid}',
                    shell=True, capture_output=True)
            except Exception:
                pass
            self._proc.wait(1000)
            self._proc = None

        self._btn_run.setEnabled(True)
        self._btn_stop.setEnabled(False)
        self._btn_stop.setStyleSheet(
            "QPushButton{background:#888;color:#ddd;font-size:18px;padding:8px}")
        self._unlock_face_controls()
        self._log("[Stopped]")

    def closeEvent(self, event):
        self._stop()
        event.accept()


def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    w = MainWindow()
    w.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
