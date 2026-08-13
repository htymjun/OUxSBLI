cmake_minimum_required(VERSION 3.15)

project(OUxSBLI LANGUAGES Fortran)

include("${CMAKE_CURRENT_SOURCE_DIR}/../CMakeLists.txt")

set(CASE_EXTRA_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/../src/set_bc_tbl_sbli.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../src/set_compressible_bl.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src/set_init_common.f90")
set(CASE_MAIN_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../../src/sbli.f90")
add_ouxsbli_case()
