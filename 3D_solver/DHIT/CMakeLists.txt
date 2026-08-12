cmake_minimum_required(VERSION 3.15)

project(OUxSBLI LANGUAGES Fortran)


set(CASE_EXTRA_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/set_init_dhit.f90")
set(CASE_PHYSICAL_QUANTITIES_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/calc_physical_quantities_dhit.f90")
set(CASE_LINK_OPTIONS "-cudalib=cufft")
include("${CMAKE_CURRENT_SOURCE_DIR}/../CMakeLists.txt")
add_ouxsbli_case()
