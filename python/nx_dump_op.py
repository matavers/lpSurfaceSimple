#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
nx_dump_op.py — 把指定工序的关键参数 dump 到文件，用于定位"侧刃驱动体"缺失的隐藏参数。

用法：在 NX 里 Developer→Journal 运行，然后看 output/nx_dump.log
"""

import NXOpen
import NXOpen.CAM

DUMP = r"D:\Projects\lpSurface\Simple\output\nx_dump.log"
OP_NAME = "VARIABLE_CONTOUR"   # 改成你手动建的正确工序名


def log(msg):
    with open(DUMP, "a", encoding="utf-8") as f:
        f.write(str(msg) + "\n")


def main():
    try:
        open(DUMP, "w", encoding="utf-8").close()
    except Exception:
        pass

    session = NXOpen.Session.GetSession()
    workPart = session.Parts.Work

    op = workPart.CAMSetup.CAMOperationCollection.FindObject(OP_NAME)
    if op is None:
        log(f"没找到工序 {OP_NAME}")
        return
    log(f"=== 工序 {OP_NAME} ===")

    builder = workPart.CAMSetup.CAMOperationCollection.CreateSurfaceContourBuilder(op)

    # 刀轴
    tav = builder.ToolAxisVariable
    log(f"ToolAxisType = {tav.ToolAxisType}")
    for name in ["RulingType", "SwarfTiltAngle", "LeadAngle", "TiltAngle",
                 "ApplySmoothing", "FanDistance", "RotationAngle"]:
        try:
            v = getattr(tav, name)
            log(f"ToolAxisVariable.{name} = {v}")
        except Exception as e:
            log(f"ToolAxisVariable.{name} = <err {e}>")

    # 驱动方法
    try:
        log(f"DriveMethod = {builder.DriveMethod}")
    except Exception as e:
        log(f"DriveMethod = <err {e}>")

    # 投影矢量
    try:
        pv = builder.ProjectionVector
        log(f"ProjectionVector = {pv}")
        for name in ["Method", "Vector", "Plane"]:
            try:
                log(f"ProjectionVector.{name} = {getattr(pv, name)}")
            except Exception as e:
                log(f"ProjectionVector.{name} = <err {e}>")
    except Exception as e:
        log(f"ProjectionVector = <err {e}>")

    # 曲面区域驱动方法 (DmSurfBuilder)
    try:
        dm = builder.DmSurfBuilder
        log("=== DmSurfBuilder ===")
        for name in ["CutDirection", "CutPattern", "Stepover", "PatternType",
                     "CutLevel", "MaterialSide", "AreaMethod", "CutDirectionType"]:
            try:
                v = getattr(dm, name)
                log(f"DmSurfBuilder.{name} = {v}")
            except Exception:
                pass
        # 驱动几何
        try:
            dg = dm.DriveGeometry
            gl = dg.GeometryList
            log(f"DriveGeometry items = {gl.Length}")
        except Exception as e:
            log(f"DriveGeometry = <err {e}>")
    except Exception as e:
        log(f"DmSurfBuilder = <err {e}>")

    builder.Destroy()
    log("=== done ===")


if __name__ == "__main__":
    main()
