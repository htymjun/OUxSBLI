cmake_minimum_required(VERSION 3.15)

# This file is <root>/3D_solver_curv/CMakeLists.txt; a case's CMakeLists.txt normally sets
# OUXSBLI_ROOT before including it (see any case directory).
if(NOT DEFINED OUXSBLI_ROOT)
  get_filename_component(OUXSBLI_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
endif()

if(NOT DEFINED CMAKE_Fortran_COMPILER)
  set(CMAKE_Fortran_COMPILER mpif90)
endif()

if(NOT DEFINED CASE_GPU_CC)
  set(CASE_GPU_CC "native")
endif()

set(SUBROUTINES "name:calc_tau_straight,calc_tau_cross,flux4")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cuda -acc -fast -Mpreprocess -Mchkptr -Minline=reshape,${SUBROUTINES} -Minfo=accel,inline -gpu=ptxinfo,rdc,lto,cc${CASE_GPU_CC},lineinfo -lnvhpcwrapnvtx")

find_package(MPI REQUIRED Fortran)
find_program(FYPP_EXE fypp REQUIRED)

# Macro: call from each case's CMakeLists.txt after setting CASE_EXTRA_SOURCES.
# CMAKE_CURRENT_SOURCE_DIR / BINARY_DIR resolve to the caller's case directory, which may be
# anywhere; the shared sources are taken from OUXSBLI_ROOT.
macro(add_ouxsbli_curv_case)
  if(DEFINED CASE_DIR)
    set(_CASE_DIR "${CASE_DIR}")
  else()
    set(_CASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()
  set(FYPP_FLAGS "-I${_CASE_DIR}" "-DDIM=3")

  set(_PARENT_FYPP "mod_constant" "mod_prec" "calc_muscl" "print")
  set(_SHARED_FYPP
      "calc_flux_base_curv"
      "preprocess_curv"
      "calc_time_dev_curv"
      "main_curv")

  # calc_scheme_math.f90.fypp is include'd as raw text (no `module` wrapper) by
  # calc_keep_kernel_curv/calc_slau_kernel_curv/calc_hybrid_kernel_curv, not
  # compiled as its own translation unit -- generate it via a custom target so
  # it exists before those files compile, without adding it to _GEN/add_executable.
  set(_SCHEME_MATH_IN  "${OUXSBLI_ROOT}/src/calc_scheme_math.f90.fypp")
  set(_SCHEME_MATH_OUT "${CMAKE_CURRENT_BINARY_DIR}/calc_scheme_math.f90")
  add_custom_command(OUTPUT ${_SCHEME_MATH_OUT}
    COMMAND ${FYPP_EXE} ${FYPP_FLAGS} ${_SCHEME_MATH_IN} ${_SCHEME_MATH_OUT}
    DEPENDS ${_SCHEME_MATH_IN} "${_CASE_DIR}/config.fypp"
    COMMENT "fypp ../../src: calc_scheme_math (include-only)")

  set(_GEN "")
  foreach(_f ${_PARENT_FYPP})
    set(_in  "${OUXSBLI_ROOT}/src/${_f}.f90.fypp")
    set(_out "${CMAKE_CURRENT_BINARY_DIR}/${_f}.f90")
    add_custom_command(OUTPUT ${_out}
      COMMAND ${FYPP_EXE} ${FYPP_FLAGS} ${_in} ${_out}
      DEPENDS ${_in} "${_CASE_DIR}/config.fypp"
      COMMENT "fypp ../../src: ${_f}")
    list(APPEND _GEN ${_out})
  endforeach()

  foreach(_f ${_SHARED_FYPP})
    set(_in  "${OUXSBLI_ROOT}/3D_solver_curv/src/${_f}.f90.fypp")
    set(_out "${CMAKE_CURRENT_BINARY_DIR}/${_f}.f90")
    add_custom_command(OUTPUT ${_out}
      COMMAND ${FYPP_EXE} ${FYPP_FLAGS} ${_in} ${_out}
      DEPENDS ${_in} "${_CASE_DIR}/config.fypp"
      COMMENT "fypp ../src: ${_f}")
    list(APPEND _GEN ${_out})
  endforeach()

  set(_BASE
      "${_CASE_DIR}/mod_globals.f90"
      "${_CASE_DIR}/set.f90"
      "${OUXSBLI_ROOT}/src/cpu_gpu_mpi.f90"
      "${OUXSBLI_ROOT}/src/set_coordinate.f90"
      "${OUXSBLI_ROOT}/src/calc_physical_quantities.f90"
      "${OUXSBLI_ROOT}/src/print_curv.f90"
      "${OUXSBLI_ROOT}/3D_solver/src/calc_hybrid.f90"
      "${OUXSBLI_ROOT}/3D_solver/src/set_bc_common.f90"
      "${OUXSBLI_ROOT}/3D_solver_curv/src/calc_hybrid_curv.f90"
      "${OUXSBLI_ROOT}/3D_solver_curv/src/calc_keep_kernel_curv.f90"
      "${OUXSBLI_ROOT}/3D_solver_curv/src/calc_slau_kernel_curv.f90"
      "${OUXSBLI_ROOT}/3D_solver_curv/src/calc_hybrid_kernel_curv.f90"
      "${OUXSBLI_ROOT}/3D_solver_curv/src/calc_steps_curv.f90"
      "${OUXSBLI_ROOT}/3D_solver_curv/src/load_smem_visc2_curv.f90"
      "${OUXSBLI_ROOT}/3D_solver_curv/src/calc_visc2_curv.f90"
      "${OUXSBLI_ROOT}/3D_solver_curv/src/calc_para_curv.f90")

  get_filename_component(_CASE "${CMAKE_CURRENT_SOURCE_DIR}" NAME)
  string(MAKE_C_IDENTIFIER "${_CASE}" _CASE)
  add_custom_target(${_CASE}_scheme_math DEPENDS ${_SCHEME_MATH_OUT})
  add_executable(${_CASE} ${_BASE} ${CASE_EXTRA_SOURCES} ${_GEN})
  add_dependencies(${_CASE} ${_CASE}_scheme_math)
  target_link_libraries(${_CASE} PRIVATE MPI::MPI_Fortran)
  if(DEFINED CASE_LINK_OPTIONS)
    target_link_options(${_CASE} PRIVATE ${CASE_LINK_OPTIONS})
  endif()
  target_include_directories(${_CASE} PRIVATE
      "${CMAKE_CURRENT_BINARY_DIR}"
      "${OUXSBLI_ROOT}/3D_solver_curv/src"
      "${OUXSBLI_ROOT}/src"
      "${OUXSBLI_ROOT}/3D_solver/src")
  set_target_properties(${_CASE} PROPERTIES
      OUTPUT_NAME "a.out"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endmacro()
