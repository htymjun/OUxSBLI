cmake_minimum_required(VERSION 3.15)

project(OUxSBLI LANGUAGES Fortran)

include("${CMAKE_CURRENT_SOURCE_DIR}/../CMakeLists.txt")

set(CASE_EXTRA_SOURCES "")
add_ouxsbli_case()
