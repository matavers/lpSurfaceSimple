"""Generate API_DOC.docx matching API_DOC.md — ruledSurfaceFitting full reference."""
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
import os

DOC_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_PATH = os.path.join(DOC_DIR, "API_DOC.docx")

TITLE_SZ = Pt(22)
H1_SZ = Pt(16)
H2_SZ = Pt(14)
BODY_SZ = Pt(11)
CODE_SZ = Pt(9)


def run_font(run, name="宋体", sz=BODY_SZ, bold=False):
    run.font.name = name; run.font.size = sz; run.bold = bold
    r = run._element; rPr = r.find(qn('w:rPr'))
    if rPr is None: rPr = r.makeelement(qn('w:rPr'), {}); r.insert(0, rPr)
    rF = rPr.find(qn('w:rFonts'))
    if rF is None: rF = rPr.makeelement(qn('w:rFonts'), {}); rPr.append(rF)
    for a in ['w:eastAsia','w:ascii','w:hAnsi']: rF.set(qn(a), name)


def para(doc, text, sz=BODY_SZ, bold=False, align=None):
    p = doc.add_paragraph()
    if align is not None: p.alignment = align
    run = p.add_run(text); run_font(run, sz=sz, bold=bold)
    return p


def heading(doc, text, sz):
    p = doc.add_paragraph()
    run = p.add_run(text); run_font(run, sz=sz, bold=True)
    pf = p.paragraph_format; pf.space_before = Pt(12); pf.space_after = Pt(4)
    return p


def code(doc, text):
    for ln in text.strip().split('\n'):
        p = doc.add_paragraph(); pf = p.paragraph_format
        pf.left_indent = Cm(1); pf.space_before = Pt(0); pf.space_after = Pt(0)
        r = p.add_run(ln); run_font(r, name="Consolas", sz=CODE_SZ)
    doc.add_paragraph()


def tbl(doc, headers, rows):
    t = doc.add_table(rows=1+len(rows), cols=len(headers)); t.style = 'Light Grid Accent 1'
    for i, h in enumerate(headers):
        c = t.rows[0].cells[i]; c.text = h
        for p in c.paragraphs:
            for r in p.runs: run_font(r, sz=Pt(10), bold=True)
    for ri, row in enumerate(rows):
        for ci, v in enumerate(row):
            c = t.rows[ri+1].cells[ci]; c.text = str(v)
            for p in c.paragraphs:
                for r in p.runs: run_font(r, sz=Pt(10))
    doc.add_paragraph()


def build():
    d = Document()
    s = d.sections[0]
    s.page_width  = Cm(21); s.page_height = Cm(29.7)
    s.left_margin = Cm(2.5); s.right_margin = Cm(2.5)
    s.top_margin  = Cm(2.5); s.bottom_margin  = Cm(2.5)

    # ── Title ──
    ti = d.add_paragraph(); ti.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = ti.add_run("ruledSurfaceFitting — 直纹面/平面逼近算法 API 技术文档")
    run_font(r, sz=TITLE_SZ, bold=True)
    d.add_paragraph()

    # ═══ §1 ═══
    heading(d, "1. 算法概述", H1_SZ)
    para(d, "本库 ruledSurfaceFitting 实现了基于 NURBS 曲面的分段直纹面拟合和分段平面拟合两种算法，用于航空发动机叶片曲面的几何简化。")

    heading(d, "1.1 支持的导入文件类型", H2_SZ)
    tbl(d, ["格式", "扩展名", "标准"],
        [["STEP", ".step, .stp", "ISO 10303 AP203/AP214"],
         ["IGES", ".igs, .iges", "ANSI US PRO/IPO-100"]])
    para(d, "支持包含多个面的文件——可在可视化 UI 中预览各面，选择任意两个面进行处理。")

    heading(d, "1.2 算法流程", H2_SZ)
    code(d, "CAD文件导入(STEP/IGES) → 面裁剪域提取 → 参数域等距分段 →\n"
             "  [直纹面模式]   等参线准线提取 → rib最小二乘优化 → 直纹面网格 + 曲线JSON导出\n"
             "  [平面模式]     PCA最佳拟合平面  → 平面四边形     + 描述JSON导出\n"
             "→ 误差CSV报告 → 元数据JSON摘要")

    heading(d, "1.3 三种调用模式", H2_SZ)
    tbl(d, ["模式", "DLL 函数", "描述"],
        [["全自动（新增）", "pressure_ruled_fitting / pressure_plane_fitting / suction_ruled_fitting / suction_plane_fitting", "单文件输入 → 自动识别叶片面 → 弦向分割叶盆/叶背 → 拟合"],
         ["直纹面逼近", "ruled_surface_fitting", "双文件输入，手动选面，分段直纹面拟合"],
         ["平面逼近",   "plane_surface_fitting", "双文件输入，手动选面，分段 PCA 平面拟合"]])

    heading(d, "1.4 直纹面算法原理", H2_SZ)
    para(d, "对于每段子曲面，在参数域裁剪范围 [uSeg0, uSeg1] × [vSeg0, vSeg1] 内：")
    para(d, "(1) 初始准线：沿分割方向取等参线作为上下两条边界曲线 C₀、C₁。")
    para(d, "(2) rib采样：沿准线方向固定参数值，在交叉方向上采样 nRibs 个中间点构成'肋条'。")
    para(d, "(3) 最小二乘优化：对每条ruling线独立求解6维线性方程组，调整采样点 C₀[i] 和 C₁[i] 使肋条点到直纹面的距离平方和最小，同时用 λ 正则化约束准线不偏离原始边界曲线。")
    para(d, "(4) 直纹面：S(p, t) = (1−t)·C₀(p) + t·C₁(p)，p 沿准线方向，t 沿交叉方向。")

    heading(d, "1.5 方向组合（直纹面模式）", H2_SZ)
    tbl(d, ["分割方向", "准线方向", "适用场景"],
        [["V","V","叶片高度方向分段，水平直纹（默认）"],
         ["V","U","高度方向分段，纵向直纹"],
         ["U","U","周向分段，纵向直纹"],
         ["U","V","周向分段，水平直纹"]])

    heading(d, "1.6 平面逼近算法原理", H2_SZ)
    para(d, "对每段子曲面，在参数域裁剪范围 [uSeg0,uSeg1] × [vSeg0,vSeg1] 内：")
    para(d, "(1) 均匀采样。取 N = nUSamples × nVSamples 个曲面点 {P_k}。")
    para(d, "(2) PCA 协方差分析。计算质心 c = (1/N)ΣP_k，协方差矩阵 C = (1/N)Σ(P_k−c)(P_k−c)^T。")
    para(d, "(3) 特征分解。C 为 3×3 实对称矩阵，最小特征值 λ_min 对应的单位特征向量 n 即为最佳拟合平面的法向。")
    para(d, "(4) 平面方程。n·(x−centroid)=0，即过质心、法向为 n 的平面。")
    para(d, "(5) 误差定义。对每个采样点 P_k，dist(P_k)=|n·(P_k−centroid)|，取 maxError 和 rmsError。")
    para(d, "(6) 网格生成。将分段的 4 个角点投影到平面上构成四边形。")

    # ═══ §2 ═══
    heading(d, "2. 接口函数", H1_SZ)

    heading(d, "2.1 pressure_ruled_fitting — 叶盆直纹面逼近（全自动）", H2_SZ)
    code(d, "RULED_API RuledFittingResult* pressure_ruled_fitting(\n"
             "    const RuledConfig* config);")
    para(d, "全自动单文件处理：加载STEP/IGES → 自动识别叶片面 → 弦向截面法分割叶盆/叶背 → 叶盆V区直纹面拟合 → OBJ网格+NURBS参数TXT导出。调用方仅需提供stepFile1（模型文件路径）和outputDir，无需手动指定 faceIdx 或 V-range。")

    heading(d, "2.2 pressure_plane_fitting — 叶盆平面逼近（全自动）", H2_SZ)
    code(d, "RULED_API RuledFittingResult* pressure_plane_fitting(\n"
             "    const PlanarConfig* config);")
    para(d, "同 pressure_ruled_fitting，使用平面逼近模式。内部自动完成识别→分割→叶盆平面PCA拟合。")

    heading(d, "2.3 suction_ruled_fitting — 叶背直纹面逼近（全自动）", H2_SZ)
    code(d, "RULED_API RuledFittingResult* suction_ruled_fitting(\n"
             "    const RuledConfig* config);")
    para(d, "全自动叶背侧直纹面拟合。")

    heading(d, "2.4 suction_plane_fitting — 叶背平面逼近（全自动）", H2_SZ)
    code(d, "RULED_API RuledFittingResult* suction_plane_fitting(\n"
             "    const PlanarConfig* config);")
    para(d, "全自动叶背侧平面拟合。")

    heading(d, "2.5 ruled_surface_fitting — 直纹面逼近（手动模式）", H2_SZ)
    code(d, "RULED_API RuledFittingResult* ruled_surface_fitting(\n"
             "    const RuledConfig* config);")
    para(d, "执行完整直纹面拟合管线：加载两个STEP/IGES文件 → 面裁剪 → V/U方向等距分段 → 等参线准线提取 → rib最小二乘优化 → OBJ网格+NURBS参数TXT导出。需要调用方指定stepFile1、stepFile2和可选的faceIdx。")

    heading(d, "2.6 plane_surface_fitting — 平面逼近（手动模式）", H2_SZ)
    code(d, "RULED_API RuledFittingResult* plane_surface_fitting(\n"
             "    const PlanarConfig* config);")
    para(d, "执行完整平面拟合管线：加载两个STEP/IGES文件 → 面裁剪 → V/U方向等距分段 → PCA最佳拟合平面 → OBJ网格+平面描述JSON导出。")

    heading(d, "2.7 free_result", H2_SZ)
    code(d, "RULED_API void free_result(RuledFittingResult* result);")
    para(d, "释放上述函数分配的内存。")

    heading(d, "2B. 自动处理流水线（新增）", H2_SZ)
    para(d, "4个新函数内部自动执行以下流水线：加载文件 → 法向聚类识别叶片面+曲率过滤平面端盖 → 弦向等参线(UIso)曲率峰检测+边界曲率阈值跨越 → V=[v1,v2]叶盆/叶背区间提取 → 直纹面/平面分段拟合 → 导出。输出以 pressure_* / suction_* 前缀命名。调用方仅需提供 stepFile1 + outputDir。无需指定 faceIdx 或 UV-range。")

    # ═══ §3 ═══
    heading(d, "3. 参数结构体", H1_SZ)

    heading(d, "3.1 RuledConfig — 直纹面逼近配置", H2_SZ)
    tbl(d, ["字段","类型","默认值","说明"],
        [["stepFile1", "const char*", "—", "第一个CAD文件路径（STEP/IGES）"],
         ["stepFile2", "const char*", "—", "第二个CAD文件路径"],
         ["outputDir", "const char*", '"."', "输出目录"],
         ["nUSamples", "int", "50", "沿准线方向采样点数"],
         ["nVSamples", "int", "10", "交叉方向采样点数"],
         ["nRibs", "int", "20", "每条ruling线肋条点数"],
         ["lambda", "double", "1.0", "正则化强度"],
         ["splitDir[2]", "RuledDirection[2]", "{V,V}", "每面分割方向"],
         ["directrixDirs[2][10]", "int[2][10]", "{{V,V,V},{V,V,V}}", "每段准线方向"],
         ["numDirectrixDirs[2]", "int[2]", "{3,3}", "每面段数"],
         ["faceIdx[2]", "int[2]", "{-1,-1}", "面索引，-1自动选最大面"]])

    heading(d, "3.2 PlanarConfig — 平面逼近配置", H2_SZ)
    tbl(d, ["字段","类型","默认值","说明"],
        [["stepFile1", "const char*", "—", "第一个CAD文件路径"],
         ["stepFile2", "const char*", "—", "第二个CAD文件路径"],
         ["outputDir", "const char*", '"."', "输出目录"],
         ["nUSamples", "int", "50", "U方向采样点数"],
         ["nVSamples", "int", "10", "V方向采样点数"],
         ["splitDir[2]", "RuledDirection[2]", "{V,V}", "每面分割方向"],
         ["faceIdx[2]", "int[2]", "{-1,-1}", "面索引，-1自动选最大面"]])

    heading(d, "3.3 RuledDirection — 方向枚举", H2_SZ)
    tbl(d, ["值","含义"],
        [["RULED_DIR_U (0)","U参数方向"],["RULED_DIR_V (1)","V参数方向"]])

    heading(d, "3.4 RuledFittingResult — 返回结果", H2_SZ)
    tbl(d, ["字段","类型","说明"],
        [["errorCode", "RuledErrorCode", "错误码，0=成功"],
         ["errorMsg[256]", "char[256]", "错误描述"],
         ["numSurfaces", "int", "处理的曲面数（固定为2）"],
         ["surfaces[2]", "RuledSurfaceResult[2]", "每面拟合结果"],
         ["metaJson[2048]", "char[2048]", "结果JSON摘要"]])

    heading(d, "3.5 RuledSurfaceResult — 单面结果", H2_SZ)
    tbl(d, ["字段","类型","说明"],
        [["name[64]", "char[64]", "曲面名称"],
         ["numSegments", "int", "段数"],
         ["segments[10]", "RuledSegmentResult[10]", "每段结果"]])

    heading(d, "3.6 RuledSegmentResult — 单段结果", H2_SZ)
    tbl(d, ["字段","类型","说明"],
        [["segmentIndex", "int", "段索引（0-based）"],
         ["maxError", "double", "最大距离误差"],
         ["rmsError", "double", "均方根误差"]])

    heading(d, "3.7 RuledErrorCode — 错误码", H2_SZ)
    tbl(d, ["值","含义"],
        [["RULED_OK (0)", "成功"],
         ["RULED_ERR_FILE_NOT_FOUND (1)", "文件路径为空"],
         ["RULED_ERR_STEP_READ_FAILED (2)", "CAD文件读取失败"],
         ["RULED_ERR_NO_VALID_FACE (3)", "未找到有效面"],
         ["RULED_ERR_EXPORT_FAILED (4)", "结果导出失败"],
         ["RULED_ERR_INVALID_PARAMS (5)", "参数为空"]])

    # ═══ §4 ═══
    heading(d, "4. 输出文件总览", H1_SZ)
    tbl(d, ["文件","格式","模式","说明"],
        [["blade1_mesh.obj","Wavefront OBJ","通用","面1原始裁剪三角网格"],
         ["blade2_mesh.obj","Wavefront OBJ","通用","面2原始裁剪三角网格"],
         ["blade1_segN.obj","Wavefront OBJ","直纹面","第N段直纹面四边形网格"],
         ["blade1_segN_params.txt","TXT","直纹面","第N段NURBS参数方程（度数+控制点+节点向量+权重+映射）"],
         ["blade1_planeN.obj","Wavefront OBJ","平面","第N段拟合平面四边形"],
         ["blade1_planeN_desc.txt","TXT","平面","第N段平面质心+法向量"],
         ["errors.csv","CSV","通用","所有段误差指标表"],
         ["meta.json","JSON","通用","元数据与误差摘要"]])

    # ═══ §5 ═══
    heading(d, "5. 使用示例", H1_SZ)

    heading(d, "5.1 C语言 — 直纹面逼近", H2_SZ)
    code(d, '#include "ruledSurfaceFitting.h"\n\n'
             'int main() {\n'
             '    RuledConfig cfg = {0};\n'
             '    cfg.stepFile1 = "Blade-raw1.STEP";\n'
             '    cfg.stepFile2 = "Blade-raw2.STEP";\n'
             '    cfg.outputDir = "./output";\n'
             '    cfg.nUSamples = 50; cfg.nVSamples = 10;\n'
             '    cfg.nRibs = 20; cfg.lambda = 1.0;\n'
             '    cfg.splitDir[0] = RULED_DIR_V;\n'
             '    cfg.splitDir[1] = RULED_DIR_V;\n'
             '    cfg.numDirectrixDirs[0] = 3;\n'
             '    cfg.directrixDirs[0][0] = RULED_DIR_V;\n'
             '    // ... 设置其他 directrixDirs\n'
             '    RuledFittingResult* res = ruled_surface_fitting(&cfg);\n'
             '    if (res->errorCode == RULED_OK) {\n'
             '        printf("Seg1 maxErr: %.4f\\n",\n'
             '               res->surfaces[0].segments[1].maxError);\n'
             '    }\n'
             '    free_result(res);\n'
             '    return 0;\n'
             '}')

    heading(d, "5.2 C语言 — 平面逼近", H2_SZ)
    code(d, 'PlanarConfig cfg = {0};\n'
             'cfg.stepFile1 = "Blade.igs";\n'
             'cfg.stepFile2 = "Blade.igs";\n'
             'cfg.outputDir = "./output";\n'
             'cfg.nUSamples = 50; cfg.nVSamples = 10;\n'
             'cfg.splitDir[0] = RULED_DIR_V;\n'
             'cfg.splitDir[1] = RULED_DIR_V;\n'
             'cfg.faceIdx[0] = 0; cfg.faceIdx[1] = 1;\n'
             '\n'
             'RuledFittingResult* res = plane_surface_fitting(&cfg);\n'
             '// ... 同上处理 res->surfaces ...\n'
             'free_result(res);')

    heading(d, "5.3 Python ctypes 调用", H2_SZ)
    code(d, 'import ctypes\n\n'
             'lib = ctypes.CDLL("./ruledSurfaceFitting.dll")\n'
             'lib.ruled_surface_fitting.restype = ctypes.c_void_p\n'
             'lib.free_result.argtypes = [ctypes.c_void_p]\n\n'
             'class RuledConfig(ctypes.Structure):\n'
             '    _fields_ = [\n'
             '        ("stepFile1", ctypes.c_char_p),\n'
             '        ("stepFile2", ctypes.c_char_p),\n'
             '        ("outputDir", ctypes.c_char_p),\n'
             '        ("nUSamples", ctypes.c_int),\n'
             '        ("nVSamples", ctypes.c_int),\n'
             '        ("nRibs", ctypes.c_int),\n'
             '        ("lambda", ctypes.c_double),\n'
             '        ("splitDir", ctypes.c_int * 2),\n'
             '        ("directrixDirs", (ctypes.c_int * 10) * 2),\n'
             '        ("numDirectrixDirs", ctypes.c_int * 2),\n'
             '        ("faceIdx", ctypes.c_int * 2),\n'
             '    ]\n\n'
             'cfg = RuledConfig()\n'
             'cfg.stepFile1 = b"Blade.igs"\n'
             'cfg.outputDir = b"./output"\n'
             'ptr = lib.ruled_surface_fitting(ctypes.byref(cfg))\n'
             'lib.free_result(ptr)')

    # ═══ §6 ═══
    heading(d, "6. 构建与依赖", H1_SZ)

    heading(d, "6.1 依赖库", H2_SZ)
    para(d, "OpenCASCADE 8.0.0（TKernel、TKMath、TKGeomBase、TKGeomAlgo、TKG2d、TKG3d、TKMesh、TKBRep、TKTopAlgo、TKBO、TKService、TKXSBase、TKDESTEP、TKDEIGES、TKShHealing）")
    para(d, "Eigen 3（随 VTK 9.4 内置）")

    heading(d, "6.2 构建", H2_SZ)
    code(d, 'cmake -S . -B build -G "Visual Studio 17 2022" -A x64\n'
             'cmake --build build --config Release')
    para(d, "产物：")
    para(d, "  build/Release/ruledSurfaceFitting.dll + .lib — 动态库")
    para(d, "  build/Release/simple.exe — 命令行工具")

    heading(d, "6.3 部署", H2_SZ)
    para(d, "使用DLL时需确保以下OCCT运行时DLL在PATH中：")
    code(d, "TKernel.dll  TKMath.dll    TKGeomBase.dll  TKGeomAlgo.dll\n"
             "TKG2d.dll    TKG3d.dll     TKMesh.dll      TKBRep.dll\n"
             "TKTopAlgo.dll TKBO.dll     TKService.dll   TKXSBase.dll\n"
             "TKDESTEP.dll TKDEIGES.dll  TKShHealing.dll\n"
             "freetype.dll FreeImage.dll tbb12.dll tbbmalloc.dll jemalloc.dll")

    # ═══ §7 ═══
    heading(d, "7. 测试", H1_SZ)

    heading(d, "7.1 Python测试（推荐）", H2_SZ)
    code(d, 'D:\\anaconda\\envs\\simple\\python.exe tests\\test_api.py\n'
             'D:\\anaconda\\envs\\simple\\python.exe tests\\test_api.py --gui')

    heading(d, "7.2 C测试", H2_SZ)
    code(d, 'cl /EHsc /I"api" tests\\test_api.c /link /LIBPATH:"build\\Release" ruledSurfaceFitting.lib\n'
             'test_api.exe Blade-raw1.STEP Blade-raw2.STEP')

    # ═══ §8 ═══
    heading(d, "8. 误差计算方法", H1_SZ)

    heading(d, "8.1 直纹面逼近误差", H2_SZ)
    para(d, "对每段子曲面参数域 [uSeg0, uSeg1] × [vSeg0, vSeg1]：")
    para(d, "(1) 采样网格。沿准线方向取 nAlong 个点，交叉方向取 nAcross 个点，共 N = nAlong × nAcross 个采样点 {P_k}。")
    para(d, "(2) 点到直纹面距离。在直纹面网格上暴力搜索最近点。")
    para(d, "(3) 误差指标。maxError = max_k dist(P_k, S)，rmsError = sqrt(1/N·Σ_k dist(P_k, S)²)。其中 maxError 反映最大局部偏差，rmsError 反映整体拟合质量。")

    heading(d, "8.2 平面逼近误差", H2_SZ)
    para(d, "(1) PCA平面拟合。计算点集质心 c 和 3×3 协方差矩阵 C = (1/N)Σ(P_k−c)(P_k−c)^T。对 C 进行特征分解，最小特征值对应的特征向量即为最佳拟合平面的单位法向 n。平面方程为 n·(x−c)=0。")
    para(d, "(2) 点到平面距离。dist(P_k) = |n·(P_k − c)|。")
    para(d, "(3) 误差指标。与直纹面相同，采用 maxError 和 rmsError。")

    heading(d, "8.3 误差报告", H2_SZ)
    para(d, "运行后在 outputDir/errors.csv 中输出每段的 maxError 和 rmsError。同时在 RuledFittingResult.metaJson 中以JSON格式返回误差摘要。两个误差指标均为欧氏距离，单位与输入模型的几何单位一致（通常为毫米）。")

    # ═══ §9 ═══
    heading(d, "9. 导出文件格式详述", H1_SZ)

    heading(d, "9.1 原始网格 OBJ", H2_SZ)
    para(d, "文件名 blade1_mesh.obj / blade2_mesh.obj，Wavefront OBJ 格式，包含面裁剪域内的三角网格。")
    code(d, "v 322.17 -9.02 206.30\nv 322.17 -9.42 207.52\n...\nf 1 2 3\nf 2 4 3\n...")

    heading(d, "9.2 直纹面网格 OBJ", H2_SZ)
    para(d, "文件名 blade1_segN.obj（N=0,1,2），四边形网格，由两条准线的采样点经线性插值生成。")

    heading(d, "9.3 准线参数方程 TXT（仅直纹面模式）", H2_SZ)
    para(d, "文件名 blade1_segN_params.txt。INI格式纯文本，包含上下准线的完整 NURBS 数学表示：")
    code(d, 'nSamples = 50\n\n'
             '[C0]\ndegree = 3\nrational = no\nnbPoles = 25\npoles (x y z w):\n'
             '  322.19 -8.96 206.30 1.0\n  ...\n'
             'knots: 0.0 0.09375 0.125 ... 1.0\nmultiplicities: 4 1 1 ... 4\n\n'
             '[C1]\n...\n\n[mapping]\n'
             'description = identity, C0(u) and C1(u) share same U-parameterization')
    para(d, "degree为NURBS次数；poles为控制点(x,y,z)及权重w；knots为节点向量；multiplicities为节点重复度。mapping描述C₀到C₁的参数对应关系——恒等映射（两者共享同一U参数化）。")

    heading(d, "9.4 平面网格 OBJ（仅平面模式）", H2_SZ)
    para(d, "文件名 blade1_planeN.obj。包含4个顶点的四边形（2个三角形），是将分段参数域4角点投影到PCA最佳拟合平面所得。")

    heading(d, "9.5 平面描述 TXT（仅平面模式）", H2_SZ)
    para(d, "文件名 blade1_planeN_desc.txt。包含拟合平面的数学描述及边界四边形角点：")
    code(d, 'centroid = 317.360873 -12.783161 236.289687\n'
             'normal = -0.733008 -0.665962 -0.138541\n'
             'corner0 = 320.774019 -10.244814 206.029191\n'
             'corner1 = 314.577559 -15.952813 266.252473\n'
             'corner2 = 310.909618 -11.908336 266.217559\n'
             'corner3 = 317.347649 -6.551311 206.403261')
    para(d, "centroid为拟合平面质心坐标，normal为单位法向量。4个corner为分段参数域四角点到平面上的投影坐标，围成平面边界四边形。平面方程：normal·(x−centroid)=0。")

    heading(d, "9.6 误差表 CSV", H2_SZ)
    para(d, "文件名 errors.csv，每行记录一段的误差指标：")
    code(d, "surface,version,segment,mode,maxError,rmsError\n"
             "Blade-1,0,0,ruled,0.06369,0.13380\n"
             "Blade-1,0,1,ruled,0.04756,0.12277")

    heading(d, "9.7 元数据 JSON", H2_SZ)
    para(d, "文件名 meta.json，包含模式标识、输入文件列表、每面每段的误差摘要：")
    code(d, '{\n  "mode": "ruled",\n  "files": ["Blade-raw1.STEP", "Blade-raw2.STEP"],\n'
             '  "surfaces": [{"name":"Blade-1","segments":[{"index":0,"maxErr":0.064,...}]}]\n}')
    para(d, "注意：可视化 UI 仅加载 .obj 文件用于 3D 渲染。.txt、.json 和 .csv 文件不会被 UI 读取，仅供外部数据分析使用。")

    # ═══ §10 ═══
    heading(d, "10. 变更记录", H1_SZ)
    tbl(d, ["版本", "日期", "变更"],
        [["v2.0", "2026-08-09", "新增4个全自动函数：pressure_ruled_fitting、pressure_plane_fitting、suction_ruled_fitting、suction_plane_fitting。单文件输入自动完成识别→分割→拟合。"],
         ["v1.0", "2026-07", "初始版本：ruled_surface_fitting、plane_surface_fitting 双文件手动模式。"]])

    # ── Save ──
    d.save(OUT_PATH)
    print(f"Written: {OUT_PATH}")


if __name__ == "__main__":
    build()
