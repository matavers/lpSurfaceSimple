cmake_minimum_required(VERSION 3.16)
project(simple VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# OCCT root (relative to Simple, pointing to ../OCCT/opencascade-8.0.0-vc14-64-combined)
set(OCCT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../OCCT/opencascade-8.0.0-vc14-64-combined/opencascade-8.0.0-vc14-64"
    CACHE PATH "Path to OCCT root")
set(OCCT_BIN_DIR "${OCCT_ROOT}/win64/vc14/bin")
set(OCCT_LIB_DIR "${OCCT_ROOT}/win64/vc14/lib")
set(OCCT_INC_DIR "${OCCT_ROOT}/inc")

set(OCCT_3RDPARTY "${CMAKE_CURRENT_SOURCE_DIR}/../OCCT/opencascade-8.0.0-vc14-64-combined/3rdparty-vc14-64")

set(OCCT_MODULES
    TKernel TKMath TKGeomBase TKGeomAlgo
    TKG2d TKG3d TKMesh TKBRep TKTopAlgo
    TKBO TKService TKXSBase TKDESTEP TKShHealing
    TKDEIGES TKXSDRAW
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
target_link_libraries(OCCT INTERFACE ${OCCT_LIBRARIES})

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

add_custom_command(TARGET simple_dll POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Copying OCCT runtime DLLs to DLL output..."
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${OCCT_BIN_DIR}"
        "$<TARGET_FILE_DIR:simple_dll>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${OCCT_3RDPARTY}/freetype-2.13.3-x64/bin/freetype.dll"
        "$<TARGET_FILE_DIR:simple_dll>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${OCCT_3RDPARTY}/freeimage-3.18.0-x64/bin/FreeImage.dll"
        "$<TARGET_FILE_DIR:simple_dll>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${OCCT_3RDPARTY}/tbb-2021.13.0-x64/bin/tbb12.dll"
        "$<TARGET_FILE_DIR:simple_dll>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${OCCT_3RDPARTY}/tbb-2021.13.0-x64/bin/tbbmalloc.dll"
        "$<TARGET_FILE_DIR:simple_dll>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${OCCT_3RDPARTY}/jemalloc-vc14-64/bin/jemalloc.dll"
        "$<TARGET_FILE_DIR:simple_dll>"
    COMMENT "Deploying runtime dependencies for DLL"
)

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

add_custom_command(TARGET simple POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Copying OCCT runtime DLLs to EXE output..."
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${OCCT_BIN_DIR}"
        "$<TARGET_FILE_DIR:simple>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${OCCT_3RDPARTY}/freetype-2.13.3-x64/bin/freetype.dll"
        "$<TARGET_FILE_DIR:simple>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${OCCT_3RDPARTY}/freeimage-3.18.0-x64/bin/FreeImage.dll"
        "$<TARGET_FILE_DIR:simple>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${OCCT_3RDPARTY}/tbb-2021.13.0-x64/bin/tbb12.dll"
        "$<TARGET_FILE_DIR:simple>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${OCCT_3RDPARTY}/tbb-2021.13.0-x64/bin/tbbmalloc.dll"
        "$<TARGET_FILE_DIR:simple>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${OCCT_3RDPARTY}/jemalloc-vc14-64/bin/jemalloc.dll"
        "$<TARGET_FILE_DIR:simple>"
    COMMENT "Deploying runtime dependencies for EXE"
)
