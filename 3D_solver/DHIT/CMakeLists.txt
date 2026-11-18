cmake_minimum_required(VERSION 3.15)

project(OUxSBLI LANGUAGES Fortran)


set(CASE_EXTRA_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/set_init_dhit.f90")
set(CASE_PHYSICAL_QUANTITIES_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/calc_physical_quantities_dhit.f90")
set(CASE_LINK_OPTIONS "-cudalib=cufft")
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
add_ouxsbli_case()
