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

# Echo the target architecture. A stale cached CASE_GPU_CC that does not match
# the GPU in the machine builds fine and then fails at launch with
# cudaErrorInvalidPtx (218) -- the `lto` device code cannot be JIT'd across
# architectures. That failure was silent until weno_micro.f90 grew a
# cudaGetLastError check; this message makes the cause visible in build.log.
message(STATUS "weno_micro target architecture: ${_CASE_GPU_CC_FLAG} (CASE_GPU_CC=${CASE_GPU_CC})")


# maxregcount raised 96 -> 128 for the WENO9-Z kernels. At 96 they pinned to
# exactly the cap and spilled 64-84 bytes, which injects LDL/STL traffic into
# the very instruction stream this benchmark measures. At 128 they settle at
# their natural 114 registers with zero spill, and the WENO5/7 kernels are
# unaffected (61/64/60 registers either way, verified) so prior WENO5 timings
# remain comparable. The solver keeps 96; this is a microbenchmark-only change.
set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} -cuda -fast -Mfma -Mpreprocess -Mchkptr -Minfo=accel -gpu=ptxinfo,rdc,lto,${_CASE_GPU_CC_FLAG},lineinfo,maxregcount:128")

# ../src/fltflt.f90 is the double-float (DF) emulation module used by the
# family G modes. Only this one file is compiled -- fltflt_operator.f90,
# fltflt_subroutines.f90 and their two interface files are textually `include`d
# by weno_micro.f90 so nvfortran can inline the bodies without crossing a
# module boundary (see the rationale in src/fltflt.f90). That is why the
# include directory below is load-bearing rather than vestigial.
add_executable(weno_micro ../src/fltflt.f90 weno_micro.f90)

target_include_directories(weno_micro PRIVATE ../src)
set_target_properties(weno_micro PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
