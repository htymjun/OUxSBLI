cmake_minimum_required(VERSION 3.15)

# This file is <root>/2D_solver/CMakeLists.txt; a case's CMakeLists.txt normally sets
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
macro(add_ouxsbli_2d_case)
  if(DEFINED CASE_DIR)
    set(_CASE_DIR "${CASE_DIR}")
  else()
    set(_CASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()
  # -m re: the FP32 literal filter of src/prec.fypp
  set(FYPP_FLAGS "-I${_CASE_DIR}" "-I${OUXSBLI_ROOT}/2D_solver/src" "-I${OUXSBLI_ROOT}/src" "-DDIM=2" "-m" "re")

  # templates of the shared src/
  set(_PARENT_FYPP "mod_constant" "mod_prec" "calc_weno" "calc_muscl" "print")
  # templates of 2D_solver/src/
  set(_SHARED_FYPP
      "calc_physical_quantities"
      "calc_steps"
      "calc_hybrid"
      "calc_div"
      "load_smem_visc2"
      "load_smem_visc4"
      "set_bc_common"
      "calc_flux_base"
      "calc_time_dev"
      "preprocess"
      "calc_keep_kernel"
      "calc_keep_kernel_internal"
      "calc_visc2"
      "calc_visc4"
      "calc_visc4_internal"
      "calc_slau_kernel"
      "calc_slau_kernel_internal"
      "calc_hybrid_kernel"
      "calc_hybrid_kernel_internal")

  # calc_scheme_math (shared src/) and calc_visc_me4_base (2D_solver/src/) are include'd as raw
  # text (no `module` wrapper) by the kernel templates, not compiled as their own translation
  # units -- generate them via a custom target so they exist before those files compile,
  # without adding them to _GEN/add_executable.
  set(_PREC_FYPP "${OUXSBLI_ROOT}/src/prec.fypp")
  set(_SCHEME_MATH_IN  "${OUXSBLI_ROOT}/src/calc_scheme_math.f90.fypp")
  set(_SCHEME_MATH_OUT "${CMAKE_CURRENT_BINARY_DIR}/calc_scheme_math.f90")
  add_custom_command(OUTPUT ${_SCHEME_MATH_OUT}
    COMMAND ${FYPP_EXE} ${FYPP_FLAGS} ${_SCHEME_MATH_IN} ${_SCHEME_MATH_OUT}
    DEPENDS ${_SCHEME_MATH_IN} ${_PREC_FYPP} "${_CASE_DIR}/config.fypp"
    COMMENT "fypp src: calc_scheme_math (include-only)")
  set(_VISC_BASE_IN  "${OUXSBLI_ROOT}/2D_solver/src/calc_visc_me4_base.f90.fypp")
  set(_VISC_BASE_OUT "${CMAKE_CURRENT_BINARY_DIR}/calc_visc_me4_base.f90")
  add_custom_command(OUTPUT ${_VISC_BASE_OUT}
    COMMAND ${FYPP_EXE} ${FYPP_FLAGS} ${_VISC_BASE_IN} ${_VISC_BASE_OUT}
    DEPENDS ${_VISC_BASE_IN} ${_PREC_FYPP} "${_CASE_DIR}/config.fypp"
    COMMENT "fypp 2D_solver/src: calc_visc_me4_base (include-only)")

  set(_GEN "")
  foreach(_f ${_PARENT_FYPP})
    set(_in  "${OUXSBLI_ROOT}/src/${_f}.f90.fypp")
    set(_out "${CMAKE_CURRENT_BINARY_DIR}/${_f}.f90")
    add_custom_command(OUTPUT ${_out}
      COMMAND ${FYPP_EXE} ${FYPP_FLAGS} ${_in} ${_out}
      DEPENDS ${_in} ${_PREC_FYPP} "${_CASE_DIR}/config.fypp"
      COMMENT "fypp src: ${_f}")
    list(APPEND _GEN ${_out})
  endforeach()

  foreach(_f ${_SHARED_FYPP})
    set(_in  "${OUXSBLI_ROOT}/2D_solver/src/${_f}.f90.fypp")
    set(_out "${CMAKE_CURRENT_BINARY_DIR}/${_f}.f90")
    add_custom_command(OUTPUT ${_out}
      COMMAND ${FYPP_EXE} ${FYPP_FLAGS} ${_in} ${_out}
      DEPENDS ${_in} ${_PREC_FYPP} "${_CASE_DIR}/config.fypp"
      COMMENT "fypp 2D_solver/src: ${_f}")
    list(APPEND _GEN ${_out})
  endforeach()

  set(_BASE
      "${_CASE_DIR}/mod_globals.f90"
      "${_CASE_DIR}/set.f90"
      "${OUXSBLI_ROOT}/src/cpu_gpu_mpi.f90"
      "${OUXSBLI_ROOT}/src/set_coordinate.f90"
      "${OUXSBLI_ROOT}/2D_solver/src/main.f90")

  get_filename_component(_CASE "${CMAKE_CURRENT_SOURCE_DIR}" NAME)
  string(MAKE_C_IDENTIFIER "${_CASE}" _CASE)
  add_custom_target(${_CASE}_scheme_math DEPENDS ${_SCHEME_MATH_OUT} ${_VISC_BASE_OUT})
  add_executable(${_CASE} ${_BASE} ${CASE_EXTRA_SOURCES} ${_GEN})
  add_dependencies(${_CASE} ${_CASE}_scheme_math)
  # PREC (config.fypp) selects the state type of set.f90 / mod_globals.f90 through cpp:
  #   REAL_T state (QJ_*), REAL_W metrics (Jacobian); PREC_DF / PREC_FP32 for the rest.
  # PREC = 'DF' adds FmaX fltflt and builds at -O1 -Mflushz: nvfortran's -fast gave wrong
  # values in fltflt kernels (see "Precision modes" in 3D_solver/CLAUDE.md).
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_CASE_DIR}/config.fypp")
  file(STRINGS "${_CASE_DIR}/config.fypp" _prec REGEX "^#:set +PREC *=")
  if(_prec MATCHES "'DF'")
    target_sources(${_CASE} PRIVATE "${OUXSBLI_ROOT}/extern/FmaX/FmaX/fltflt.f90")
    target_compile_definitions(${_CASE} PRIVATE PREC_DF "REAL_T=type(fltflt)" "REAL_W=real(8)")
    string(REPLACE "-fast" "-O1 -Mflushz" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  elseif(_prec MATCHES "'FP32'")
    target_compile_definitions(${_CASE} PRIVATE PREC_FP32 "REAL_T=real(4)" "REAL_W=real(4)")
  else()
    target_compile_definitions(${_CASE} PRIVATE "REAL_T=real(8)" "REAL_W=real(8)")
  endif()
  target_link_libraries(${_CASE} PRIVATE MPI::MPI_Fortran)
  if(DEFINED CASE_LINK_OPTIONS)
    target_link_options(${_CASE} PRIVATE ${CASE_LINK_OPTIONS})
  endif()
  target_include_directories(${_CASE} PRIVATE
      "${CMAKE_CURRENT_BINARY_DIR}"
      "${OUXSBLI_ROOT}/2D_solver/src"
      "${OUXSBLI_ROOT}/src")
  set_target_properties(${_CASE} PROPERTIES
      OUTPUT_NAME "a.out"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endmacro()
