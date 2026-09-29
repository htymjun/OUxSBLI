cmake_minimum_required(VERSION 3.15)

project(OUxSBLI LANGUAGES Fortran)

# The production GPU is an NVIDIA A100 (compute capability 8.0).  The build
# node may not expose a GPU, so ccnative would otherwise fall back to sm_50.
set(CASE_GPU_CC "80")

include("${CMAKE_CURRENT_SOURCE_DIR}/../CMakeLists.txt")

set(CASE_EXTRA_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/grid_partitioned.f90")
add_ouxsbli_case()
