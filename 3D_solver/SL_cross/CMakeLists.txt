cmake_minimum_required(VERSION 3.15)

project(OUxSBLI LANGUAGES Fortran)

# The production GPU is an NVIDIA A100 (compute capability 8.0).  The build
# node may not expose a GPU, so ccnative would otherwise fall back to sm_50.
set(CASE_GPU_CC "80")

include("${CMAKE_CURRENT_SOURCE_DIR}/../CMakeLists.txt")

set(CASE_EXTRA_SOURCES "")
add_ouxsbli_case()

# Enable separable z-volume correction only for this case.
target_compile_definitions(SL_cross PRIVATE OUXSBLI_STRETCHED_Z)
