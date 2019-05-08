cmake_minimum_required (VERSION 2.6)
project (Cddlib_Test6D)
# The version number.
set (Tutorial_VERSION_MAJOR 1)
set (Tutorial_VERSION_MINOR 0)
#set ( libipopt.a )

# configure a header file to pass some of the CMake settings
# to the source code
# configure_file (
#  "${PROJECT_SOURCE_DIR}/MyNLP.hpp"
#  "${PROJECT_BINARY_DIR}/MyNLP.hpp"
#  )

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h

find_package(PkgConfig)
pkg_check_modules(EIGEN3 eigen3 REQUIRED)
#message("${Eigen3_INCLUDE_DIRS}")
#include_directories(${Eigen3_INCLUDE_DIRS})

add_executable(Cddlib_Test6D minkowski.c)

target_link_libraries (Cddlib_Test6D cdd)
