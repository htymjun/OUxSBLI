cmake_minimum_required(VERSION 3.15)

project(OUxSBLI LANGUAGES Fortran)

include("${CMAKE_CURRENT_SOURCE_DIR}/../CMakeLists.txt")

# tracer particles (TRACER = True in config.fypp)
set(CASE_EXTRA_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/../src/tracer.f90")
add_ouxsbli_case()
