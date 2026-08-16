cmake_minimum_required(VERSION 3.18)

if(NOT DEFINED CMAKE_Fortran_COMPILER)
  set(CMAKE_Fortran_COMPILER nvfortran)
endif()

project(weno_microbenchmark LANGUAGES Fortran)

if(NOT DEFINED CASE_GPU_CC)
  set(CASE_GPU_CC "native")
endif()

if(CASE_GPU_CC STREQUAL "native")
  set(_CASE_GPU_CC_FLAG "ccnative")
else()
  set(_CASE_GPU_CC_FLAG "cc${CASE_GPU_CC}")
endif()

set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} -cuda -fast -Mfma -Mpreprocess -Mchkptr -Minfo=accel -gpu=ptxinfo,rdc,lto,${_CASE_GPU_CC_FLAG},lineinfo,maxregcount:96")

add_executable(weno_micro weno_micro.f90)

target_include_directories(weno_micro PRIVATE ../src)
set_target_properties(weno_micro PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
