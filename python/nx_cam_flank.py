#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
nx_cam_flank.py — 在 NX 里自动给拟合直纹面的每个格子建一个「侧刃铣」工序并生成刀轨。

依据你录制的 journal.py 精简而来（刀轴 = SwarfDrive = 侧刃驱动体）。

前提：
  - 已打开 blade2_fitted.step（或 blade1_fitted.step，缝合后的片体零件）。
  - 已进过加工（加工环境已建好 MCS_MILL / WORKPIECE / NC_PROGRAM）。
  - 若还没进加工，脚本会自动 ApplicationSwitchImmediate("UG_APP_MANUFACTURING")。

运行：Developer → Journal → 选本脚本（或 File → Execute → NX Open）。

参数（改 main 顶部变量）：
  - TOOL_DIAMETER：刀具直径（默认 3.0 mm）
  - MAX_FACES：最多建几个工序（0 = 全部面，建议先设 3 试跑）
  - SKIP_GENERATE：True = 只建工序不生成刀轨（先试建，确认无误再 False）
"""

import NXOpen
import NXOpen.CAM

TOOL_DIAMETER = 3.0
MAX_FACES = 3          # 先试 3 个面；确认 OK 改成 0 跑全部
SKIP_GENERATE = True   # 先 True 试建；确认 OK 改成 False

# 调试日志文件（NX 里 print 进信息窗口不好找，直接写文件）
DEBUG_LOG = r"D:\Projects\lpSurface\Simple\output\nx_cam_debug.log"


def log(msg):
    with open(DEBUG_LOG, "a", encoding="utf-8") as f:
        f.write(msg + "\n")


def get_sheet_body_and_faces(workPart):
    """返回 (片体 body, 该片体的所有面)。"""
    for body in workPart.Bodies:
        try:
            if not body.IsSolidBody:      # 片体
                return body, list(body.GetFaces())
        except Exception:
            pass
    # 兜底：用 journal 里的名字
    body = workPart.Bodies.FindObject("UNPARAMETERIZED_FEATURE(1)")
    return body, list(body.GetFaces())


def create_tool(workPart, camSetup, diameter):
    machine = camSetup.CAMGroupCollection.FindObject("GENERIC_MACHINE")
    tool = camSetup.CAMGroupCollection.CreateToolWithUserName(
        machine, "mill_multi-axis", "MILL",
        NXOpen.CAM.NCGroupCollection.UseDefaultName.TrueValue, "MILL", "Mill")
    tb = camSetup.CAMGroupCollection.CreateMillToolBuilder(tool)
    tb.TlDiameterBuilder.Value = diameter
    tb.Commit()
    tb.Destroy()
    return tool


def create_flank_op(workPart, camSetup, program, method, geom_group, workgeom,
                    body, face, tool, name):
    op = camSetup.CAMOperationCollection.CreateWithUserName(
        program, method, geom_group, workgeom,
        "mill_multi-axis", "VARIABLE_CONTOUR",
        NXOpen.CAM.OperationCollection.UseDefaultName.FalseValue,
        name, name)

    # 一个 builder 会话里全部设完（驱动方法 + 驱动面 + 刀轴），最后提交一次
    b = camSetup.CAMOperationCollection.CreateSurfaceContourBuilder(op)
    b.DriveMethod = NXOpen.CAM.SurfaceContourBuilder.DriveMethodTypes.SurfaceArea

    driveGeo = b.DmSurfBuilder.DriveGeometry
    driveSet = driveGeo.GeometryList.FindItem(0)
    driveSet.Surface = face
    driveGeo.Validate()
    driveGeo.Commit()

    tav = b.ToolAxisVariable
    log(f"[{name}] 设置前 ToolAxisType = {tav.ToolAxisType}")
    tav.ToolAxisType = NXOpen.CAM.ToolAxisVariable.Types.SwarfDrive
    log(f"[{name}] 设置后 ToolAxisType = {tav.ToolAxisType}")
    op = b.Commit()
    b.Destroy()

    # 验证：Commit 后重新打开读回
    b4 = camSetup.CAMOperationCollection.CreateSurfaceContourBuilder(op)
    log(f"[{name}] Commit后读回 ToolAxisType = {b4.ToolAxisVariable.ToolAxisType}")
    b4.Destroy()

    # ---- 挂刀具 + 生成刀轨 ----
    workPart.CAMSetup.MoveObjects(NXOpen.CAM.CAMSetup.View.MachineTool,
                                  [op], tool, NXOpen.CAM.CAMSetup.Paste.Inside)
    if not SKIP_GENERATE:
        workPart.CAMSetup.GenerateToolPath([op])
    return op


def main():
    # 清空上次的调试日志
    try:
        open(DEBUG_LOG, "w", encoding="utf-8").close()
    except Exception:
        pass

    session = NXOpen.Session.GetSession()
    workPart = session.Parts.Work
    session.ApplicationSwitchImmediate("UG_APP_MANUFACTURING")

    camSetup = workPart.CAMSetup
    program = camSetup.CAMGroupCollection.FindObject("NC_PROGRAM")
    method = camSetup.CAMGroupCollection.FindObject("METHOD")
    geom_group = camSetup.CAMGroupCollection.FindObject("NONE")
    workgeom = camSetup.CAMGroupCollection.FindObject("WORKPIECE")

    body, faces = get_sheet_body_and_faces(workPart)
    if not faces:
        raise RuntimeError("没有找到片体面，请确认已打开拟合面零件")
    print(f"[nx] 找到 {len(faces)} 个面")

    tool = create_tool(workPart, camSetup, TOOL_DIAMETER)

    n = len(faces) if MAX_FACES == 0 else min(MAX_FACES, len(faces))
    for i in range(n):
        name = f"FLANK_{i:03d}"
        try:
            create_flank_op(workPart, camSetup, program, method, geom_group,
                            workgeom, body, faces[i], tool, name)
            print(f"[nx] 面 {i} 工序已建: {name}")
        except Exception as e:
            print(f"[nx] 面 {i} 失败: {e}")

    print(f"[nx] done: {n} 个侧刃铣工序")


if __name__ == "__main__":
    main()
