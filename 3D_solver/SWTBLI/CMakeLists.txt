cmake_minimum_required(VERSION 3.15)

project(OUxSBLI LANGUAGES Fortran)

# OUxSBLI source tree: -DOUXSBLI_ROOT=<checkout>, else ../.. (a case inside the tree),
# else the OUXSBLI_ROOT environment variable (a case directory anywhere else).
if(NOT DEFINED OUXSBLI_ROOT)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../src/prec.fypp")
    get_filename_component(OUXSBLI_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
  else()
    set(OUXSBLI_ROOT "$ENV{OUXSBLI_ROOT}")
  endif()
endif()
if(NOT EXISTS "${OUXSBLI_ROOT}/src/prec.fypp")
  message(FATAL_ERROR "OUxSBLI source tree not found (OUXSBLI_ROOT='${OUXSBLI_ROOT}'): "
                      "pass -DOUXSBLI_ROOT=<checkout> or set the OUXSBLI_ROOT environment variable")
endif()
include("${OUXSBLI_ROOT}/3D_solver/CMakeLists.txt")

set(CASE_EXTRA_SOURCES
    "${OUXSBLI_ROOT}/3D_solver/src/set_bc_tbl_sbli.f90"
    "${OUXSBLI_ROOT}/src/set_compressible_bl.f90"
    "${OUXSBLI_ROOT}/3D_solver/src/set_init_common.f90")
set(CASE_MAIN_SOURCE "${OUXSBLI_ROOT}/src/sbli.f90")
add_ouxsbli_case()
