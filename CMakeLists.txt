cmake_minimum_required(VERSION 3.16)
project(simple VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# OCCT static build root (relative to Simple, pointing to Simple/3rdPart/OCCT-build-static)
set(OCCT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/3rdPart/OCCT-build-static"
    CACHE PATH "Path to OCCT static build root")
set(OCCT_LIB_DIR "${OCCT_ROOT}/win64/vc14/lib")
# Headers from the combined package (OCCT 8.0.0, complete, includes deprecated
# collection aliases like TColgp_Array1OfPnt.hxx absent from the build-tree inc).
set(OCCT_INC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../OCCT/opencascade-8.0.0-vc14-64-combined/opencascade-8.0.0-vc14-64/inc")

set(OCCT_3RDPARTY "${CMAKE_CURRENT_SOURCE_DIR}/../OCCT/opencascade-8.0.0-vc14-64-combined/3rdparty-vc14-64")

set(OCCT_MODULES
    TKernel TKMath TKG2d TKG3d TKGeomBase TKBRep
    TKGeomAlgo TKTopAlgo TKPrim TKBO TKShHealing TKBool
    TKHLR TKHelix TKFillet TKOffset TKFeat TKMesh TKXMesh TKExpress
    TKService TKV3d TKOpenGl TKMeshVS
    TKCDF TKLCAF TKCAF TKBinL TKXmlL TKBin TKXml TKStdL TKStd
    TKTObj TKBinTObj TKXmlTObj TKVCAF
    TKDE TKXSBase TKDESTEP TKXCAF TKDEIGES TKDESTL TKDEVRML TKRWMesh
    TKDECascade TKBinXCAF TKXmlXCAF TKDEOBJ TKDEPLY
)

foreach(mod ${OCCT_MODULES})
    find_library(OCCT_${mod}_LIB ${mod}
        PATHS "${OCCT_LIB_DIR}"
        NO_DEFAULT_PATH)
    if(NOT OCCT_${mod}_LIB)
        message(FATAL_ERROR "OCCT library ${mod} not found in ${OCCT_LIB_DIR}")
    endif()
endforeach()

set(OCCT_LIBRARIES)
foreach(mod ${OCCT_MODULES})
    list(APPEND OCCT_LIBRARIES ${OCCT_${mod}_LIB})
endforeach()

add_library(OCCT INTERFACE)
target_include_directories(OCCT INTERFACE "${OCCT_INC_DIR}")
target_compile_definitions(OCCT INTERFACE OCCT_STATIC_BUILD OCCT_NO_PLUGINS)
target_link_libraries(OCCT INTERFACE ${OCCT_LIBRARIES})
if(MSVC)
    # Windows system libraries required by OCCT static libraries.
    target_link_libraries(OCCT INTERFACE wsock32 advapi32 gdi32 user32 shell32 psapi winmm opengl32)
endif()

# Eigen (bundled in VTK)
set(EIGEN3_ROOT "${OCCT_3RDPARTY}/vtk-9.4.1-x64/include/vtk-9.4/vtkeigen")
if(NOT EXISTS "${EIGEN3_ROOT}/Eigen/Core")
    message(FATAL_ERROR "Eigen not found at ${EIGEN3_ROOT}")
endif()
add_library(Eigen3::Eigen INTERFACE IMPORTED)
set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_ROOT}"
)

# Sources
file(GLOB CORE_SOURCES "src/simple/*.cpp")
file(GLOB API_SOURCES "src/api/*.cpp")

set(SIMPLE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# ── Shared library (DLL) ─────────────────────────
add_library(simple_dll SHARED ${CORE_SOURCES} ${API_SOURCES})
target_include_directories(simple_dll PRIVATE ${SIMPLE_INCLUDE_DIR})
target_include_directories(simple_dll PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/api)
target_link_libraries(simple_dll PRIVATE OCCT Eigen3::Eigen)
set_target_properties(simple_dll PROPERTIES
    OUTPUT_NAME ruledSurfaceFitting
    PREFIX ""
    WINDOWS_EXPORT_ALL_SYMBOLS ON
)

if(MSVC)
    target_compile_options(simple_dll PRIVATE /W3 /bigobj /utf-8 /DNOMINMAX /DWIN32_LEAN_AND_MEAN)
    target_compile_definitions(simple_dll PRIVATE RULED_DLL_EXPORTS)
endif()

# ── Executable (CLI) ─────────────────────────────
add_executable(simple ${CORE_SOURCES})
target_include_directories(simple PRIVATE ${SIMPLE_INCLUDE_DIR})
target_link_libraries(simple PRIVATE OCCT Eigen3::Eigen)

if(MSVC)
    target_compile_options(simple PRIVATE /W3 /bigobj /utf-8 /DNOMINMAX /DWIN32_LEAN_AND_MEAN)
    set_target_properties(simple PROPERTIES
        VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:simple>"
    )
endif()

# OCCT static libraries are self-contained (no runtime DLL deployment needed).
