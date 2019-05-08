cmake_minimum_required (VERSION 2.6)
project (minksum_test)
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

#set(MINKSUM_DIRECTORY CACHE STRING)
message("---------------- $ENV{MINKSUM_DIRECTORY}")
# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}" $ENV{MINKSUM_DIRECTORY}/lib-src 
											$ENV{MINKSUM_DIRECTORY}/external/include 
											$ENV{MINKSUM_DIRECTORY}/external/lib
											/usr/local/include
											/usr/include)

add_executable(minksum_test minkSum.cc)

target_link_libraries (minksum_test MINKSUM gmp gmpxx cdd)
