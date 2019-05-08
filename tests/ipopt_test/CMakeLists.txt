cmake_minimum_required (VERSION 2.6)
project (Cpp_Example)
# The version number.
set (Tutorial_VERSION_MAJOR 1)
set (Tutorial_VERSION_MINOR 0)
set ( libipopt.a )

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/MyNLP.hpp"
  "${PROJECT_BINARY_DIR}/MyNLP.hpp"
  )

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}" /usr/local/include/coin/
											 ~/catkin_ws/src/dls/common_utilities/contact_wrench_set_margin/tests/ipopt_test/)


add_executable(Cpp_Example cpp_example.cpp
			   			    MyNLP.cpp
			   			    mink_glpk.cpp
)
#add_directories(/usr/local/lib/)
#add_directories(/usr/local/include/coin/)
target_link_libraries (Cpp_Example ipopt glpk)
