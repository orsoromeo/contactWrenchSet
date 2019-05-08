cmake_minimum_required (VERSION 2.6)
project (Cws_Thread CXX)

SET(CMAKE_CXX_FLAGS "-std=c++0x")


# The version number.
set (Tutorial_VERSION_MAJOR 1)
set (Tutorial_VERSION_MINOR 0)
#set ( libipopt.a )


# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h

find_package(PkgConfig)
pkg_check_modules(EIGEN3 eigen3 REQUIRED)
#message("${Eigen3_INCLUDE_DIRS}")
#include_directories(${Eigen3_INCLUDE_DIRS})

add_executable(Cws_Thread cws_thread.c)
target_compile_features(Cws_Thread PRIVATE cxx_range_for)
