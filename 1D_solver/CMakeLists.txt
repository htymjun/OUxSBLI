cmake_minimum_required(VERSION 3.15)

if(NOT DEFINED CMAKE_Fortran_COMPILER)
  set(CMAKE_Fortran_COMPILER mpif90)
endif()

if(NOT DEFINED CASE_GPU_CC)
  set(CASE_GPU_CC "native")
endif()

if(CASE_GPU_CC STREQUAL "native")
  set(_CASE_GPU_CC_FLAG "ccnative")
else()
  set(_CASE_GPU_CC_FLAG "cc${CASE_GPU_CC}")
endif()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cuda -acc -fast -Mfma -Mpreprocess -Mchkptr -Minfo=accel -gpu=ptxinfo,rdc,lto,${_CASE_GPU_CC_FLAG},lineinfo,maxregcount:96")

find_package(MPI REQUIRED Fortran)
find_program(FYPP_EXE fypp REQUIRED)

# Macro: call from each case's CMakeLists.txt after setting CASE_EXTRA_SOURCES.
# CMAKE_CURRENT_SOURCE_DIR / BINARY_DIR resolve to the caller's case directory.
macro(add_ouxsbli_1d_case)
  if(DEFINED CASE_DIR)
    set(_CASE_DIR "${CASE_DIR}")
  else()
    set(_CASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()
  set(FYPP_FLAGS "-I${_CASE_DIR}")

  # mod_constant and calc_muscl come from the shared repo-root src/
  set(_PARENT_FYPP "mod_constant" "calc_muscl")
  set(_SHARED_FYPP
      "calc_flux_base"
      "calc_time_dev"
      "preprocess"
      "calc_keep_kernel"
      "calc_slau_kernel"
      "calc_visc"
      "calc_seq_kernel"
      "calc_fused_kernel"
      "calc_warp_kernel"
      "calc_steps")

  # Fortran `include` fragments (device-function bodies shared by several
  # modules). These are fypp-generated but never compiled on their own, so they
  # are kept out of _GEN and instead made a DEPENDS of every generated module.
  set(_FRAGMENT_FYPP "calc_keep_1d" "calc_keep_1d_df" "calc_visc_1d" "calc_slau_1d")

  set(_FRAG "")
  foreach(_f ${_FRAGMENT_FYPP})
    set(_in  "${CMAKE_CURRENT_SOURCE_DIR}/../src/${_f}.f90.fypp")
    set(_out "${CMAKE_CURRENT_BINARY_DIR}/${_f}.f90")
    add_custom_command(OUTPUT ${_out}
      COMMAND ${FYPP_EXE} ${FYPP_FLAGS} ${_in} ${_out}
      DEPENDS ${_in} "${_CASE_DIR}/config.fypp"
      COMMENT "fypp ../src (fragment): ${_f}")
    list(APPEND _FRAG ${_out})
  endforeach()

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
      DEPENDS ${_in} "${_CASE_DIR}/config.fypp" ${_FRAG}
      COMMENT "fypp ../src: ${_f}")
    list(APPEND _GEN ${_out})
  endforeach()

  set(_BASE
      "${_CASE_DIR}/mod_globals.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/fltflt.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/set.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/calc_physical_quantities.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src/calc_weno.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/print_1d.f90"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/main.f90")

  get_filename_component(_CASE "${CMAKE_CURRENT_SOURCE_DIR}" NAME)
  add_executable(${_CASE} ${_BASE} ${CASE_EXTRA_SOURCES} ${_GEN})
  target_link_libraries(${_CASE} PRIVATE MPI::MPI_Fortran)
  target_include_directories(${_CASE} PRIVATE
      "${CMAKE_CURRENT_BINARY_DIR}"
      "${CMAKE_CURRENT_SOURCE_DIR}/../src"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../src")
  set_target_properties(${_CASE} PROPERTIES
      OUTPUT_NAME "a.out"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endmacro()
