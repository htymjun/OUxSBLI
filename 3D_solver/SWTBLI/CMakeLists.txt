cmake_minimum_required(VERSION 3.15)

project(OUxSBLI LANGUAGES Fortran)

include("${CMAKE_CURRENT_SOURCE_DIR}/../CMakeLists.txt")

# two coupled domains (TBL generator + SBLI) driven by src/sbli.f90; 4 MPI ranks
set(CASE_EXTRA_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/../../src/set_compressible_bl.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src/set_init_common.f90")
set(CASE_MAIN_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../../src/sbli.f90")
add_ouxsbli_case()
