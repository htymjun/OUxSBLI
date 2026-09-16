cmake_minimum_required(VERSION 3.15)

if(NOT DEFINED CMAKE_Fortran_COMPILER)
  set(CMAKE_Fortran_COMPILER mpif90)
endif()

if(NOT DEFINED CASE_GPU_CC)
  set(CASE_GPU_CC "native")
endif()

if(NOT DEFINED CASE_MAIN_SOURCE)
  set(CASE_MAIN_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../src/main.f90")
endif()

set(SUBROUTINES "name:calc_tau_straight,calc_tau_cross,flux4")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cuda -acc -fast -Mpreprocess -Mchkptr -Minline=reshape,${SUBROUTINES} -Minfo=accel,inline -gpu=ptxinfo,rdc,lto,cc${CASE_GPU_CC},lineinfo -lnvhpcwrapnvtx")

find_package(MPI REQUIRED Fortran)
find_program(FYPP_EXE fypp REQUIRED)

# Macro: call from each case's CMakeLists.txt after setting CASE_EXTRA_SOURCES.
# CMAKE_CURRENT_SOURCE_DIR / BINARY_DIR resolve to the caller's case directory.
macro(add_ouxsbli_case)
  if(DEFINED CASE_DIR)
    set(_CASE_DIR "${CASE_DIR}")
  else()
    set(_CASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()
  set(FYPP_FLAGS "-I${_CASE_DIR}" "-I${CMAKE_CURRENT_SOURCE_DIR}/../src" "-DDIM=3")

  set(_PARENT_FYPP "mod_constant" "calc_muscl" "print")
  set(_SHARED_FYPP
      "calc_visc2"
      "calc_visc_high"
      "calc_visc_high_internal"
      "load_smem_visc_cent"
      "calc_keep_kernel"
      "calc_keep_kernel_internal"
      "calc_keep_visc_kernel"
      "calc_slau_kernel"
      "calc_slau_kernel_internal"
      "calc_slau_visc_kernel"
      "calc_hybrid_kernel"
      "calc_hybrid_kernel_internal"
      "calc_flux_base"
      "calc_time_dev"
      "calc_para"
      "preprocess")

  # calc_scheme_math.f90.fypp is include'd as raw text (no `module` wrapper) by
  # calc_keep_kernel/calc_slau_kernel/calc_hybrid_kernel(+_internal) below, not
  # compiled as its own translation unit -- generate it via a custom target so
  # it exists before those files compile, without adding it to _GEN/add_executable.
  set(_SCHEME_MATH_IN  "${CMAKE_CURRENT_SOURCE_DIR}/../../src/calc_scheme_math.f90.fypp")
  set(_SCHEME_MATH_OUT "${CMAKE_CURRENT_BINARY_DIR}/calc_scheme_math.f90")
  add_custom_command(OUTPUT ${_SCHEME_MATH_OUT}
    COMMAND ${FYPP_EXE} ${FYPP_FLAGS} ${_SCHEME_MATH_IN} ${_SCHEME_MATH_OUT}
    DEPENDS ${_SCHEME_MATH_IN} "${_CASE_DIR}/config.fypp"
    COMMENT "fypp ../../src: calc_scheme_math (include-only)")

  set(_GEN "")
  foreach(_f ${_PARENT_FYPP})
    set(_in  "${CMAKE_CURRENT_SOURCE_DIR}/../../src/${_f}.f90.fypp")
    set(_out "${CMAKE_CURRENT_BINARY_DIR}/${_f}.f90")
    add_custom_command(OUTPUT ${_out}
      COMMAND ${FYPP_EXE} ${FYPP_FLAGS} ${_in} ${_out}
      DEPENDS ${_in} "${_CASE_DIR}/config.fypp"
      COMMENT "fypp ../../src: ${_f}")
    list(APPEND _GEN ${_out})
  endforeach()

  foreach(_f ${_SHARED_FYPP})
    set(_in  "${CMAKE_CURRENT_SOURCE_DIR}/../src/${_f}.f90.fypp")
    set(_out "${CMAKE_CURRENT_BINARY_DIR}/${_f}.f90")
    add_custom_command(OUTPUT ${_out}
      COMMAND ${FYPP_EXE} ${FYPP_FLAGS} ${_in} ${_out}
      DEPENDS ${_in} "${_CASE_DIR}/config.fypp"
      COMMENT "fypp ../src: ${_f}")
    list(APPEND _GEN ${_out})
  endforeach()

  if(NOT DEFINED CASE_PHYSICAL_QUANTITIES_SOURCE)
    set(CASE_PHYSICAL_QUANTITIES_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../../src/calc_physical_quantities.f90")
  endif()

  set(_BASE
      "${_CASE_DIR}/mod_globals.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/set.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/cpu_gpu_mpi.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/set_coordinate.f90"
      "${CASE_PHYSICAL_QUANTITIES_SOURCE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/load_smem_visc2.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/load_smem_visc_me4_base.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_div.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_les.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_hybrid.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_steps.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_rescale.f90"
      #"${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_para.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/set_bc_common.f90")

  get_filename_component(_CASE "${CMAKE_CURRENT_SOURCE_DIR}" NAME)
  add_custom_target(${_CASE}_scheme_math DEPENDS ${_SCHEME_MATH_OUT})
  add_executable(${_CASE} ${_BASE} ${CASE_MAIN_SOURCE} ${CASE_EXTRA_SOURCES} ${_GEN})
  add_dependencies(${_CASE} ${_CASE}_scheme_math)
  target_link_libraries(${_CASE} PRIVATE MPI::MPI_Fortran)
  # CUDA Lib
  if (DEFINED CASE_LINK_OPTIONS)
    target_link_options(${_CASE} PRIVATE ${CASE_LINK_OPTIONS})
  endif()
  target_include_directories(${_CASE} PRIVATE
      "${CMAKE_CURRENT_BINARY_DIR}"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src")
  set_target_properties(${_CASE} PROPERTIES
      OUTPUT_NAME "a.out"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endmacro()
