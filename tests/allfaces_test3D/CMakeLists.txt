cmake_minimum_required (VERSION 2.6)
project (allfaces_Test3D)
# The version number.
set (Tutorial_VERSION_MAJOR 1)
set (Tutorial_VERSION_MINOR 0)
#set ( libipopt.a )

find_package(PkgConfig)
pkg_check_modules(EIGEN3 eigen3 REQUIRED)

# configure a header file to pass some of the CMake settings
# to the source code
# configure_file (
#  "${PROJECT_SOURCE_DIR}/MyNLP.hpp"
#  "${PROJECT_BINARY_DIR}/MyNLP.hpp"
#  )

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
# include_directories("${PROJECT_BINARY_DIR}" /usr/local/include/coin/)

add_executable(allfaces_Test3D allfaces.c)

target_link_libraries (allfaces_Test3D cdd)
