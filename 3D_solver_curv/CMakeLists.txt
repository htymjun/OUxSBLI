cmake_minimum_required(VERSION 3.15)

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
# CMAKE_CURRENT_SOURCE_DIR / BINARY_DIR resolve to the caller's case directory.
macro(add_ouxsbli_curv_case)
  if(DEFINED CASE_DIR)
    set(_CASE_DIR "${CASE_DIR}")
  else()
    set(_CASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()
  set(FYPP_FLAGS "-I${_CASE_DIR}")

  set(_PARENT_FYPP "mod_constant" "calc_muscl")
  set(_SHARED_FYPP
      "calc_flux_base_curv"
      "preprocess_curv"
      "calc_time_dev_curv"
      "main_curv")

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

  set(_BASE
      "${_CASE_DIR}/mod_globals.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/set.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/cpu_gpu_mpi.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/set_coordinate.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/calc_physical_quantities.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/print.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/print_curv.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../3D_solver/src/calc_hybrid.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../3D_solver/src/set_bc_common.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_hybrid_curv.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_keep_kernel_curv.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_slau_kernel_curv.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_hybrid_kernel_curv.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_steps_curv.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/load_smem_visc2_curv.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_visc2_curv.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/calc_para_curv.f90")

  get_filename_component(_CASE "${CMAKE_CURRENT_SOURCE_DIR}" NAME)
  add_executable(${_CASE} ${_BASE} ${CASE_EXTRA_SOURCES} ${_GEN})
  target_link_libraries(${_CASE} PRIVATE MPI::MPI_Fortran)
  if(DEFINED CASE_LINK_OPTIONS)
    target_link_options(${_CASE} PRIVATE ${CASE_LINK_OPTIONS})
  endif()
  target_include_directories(${_CASE} PRIVATE
      "${CMAKE_CURRENT_BINARY_DIR}"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../3D_solver/src")
  set_target_properties(${_CASE} PROPERTIES
      OUTPUT_NAME "a.out"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endmacro()
