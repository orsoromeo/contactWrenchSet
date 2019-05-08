cmake_minimum_required(VERSION 2.8.3)
project(contact_wrench_set)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")

find_package(catkin REQUIRED
                         dwl
                         roscpp
  			                 roslib
  			                 math_utils
  			                 dwl_rviz_plugin
  			                 std_msgs)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME}
  CATKIN_DEPENDS roscpp std_msgs dwl dwl_rviz_plugin pluginlib
)
## Adding unnecessary libraries for debugging
find_package(PkgConfig)
pkg_check_modules(EIGEN3 eigen3 REQUIRED)
include_directories(${LCM_INCLUDE_DIRS})
pkg_check_modules (LCM lcm)
#include_directories(/home/rorsolino/git/CoinIpopt/Ipopt/src/Interfaces)
#include_directories(/home/rorsolino/git/CoinIpopt/Ipopt/src/Common)
#include_directories(/home/rorsolino/git/CoinIpopt/Ipopt/src/LinAlg)
#include_directories(/home/rorsolino/git/CoinIpopt/Ipopt/src/Apps)
#include_directories(/home/rorsolino/git/CoinIpopt/Ipopt/src/Algorithm)
#include_directories(/home/rorsolino/git/CoinIpopt/Ipopt/src/contrib)

link_directories(${LCM_LIBRARY_DIRS})
link_directories(/usr/local/lib)
set(GENCODELIB iitgenhyq2max iitgenhyqreal iitgenhyq)
#this is to find robot.h for liblocomotionutils
include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(/usr/local/include/iit/commons)
include_directories(/usr/local/include/iit/robots/hyq2max)
include_directories(/usr/local/include/iit/robots/hyq)
include_directories(/usr/local/include/iit/robots/hyqreal)
set(LOC_LIB_NAME "locomotionutils${ROBOT_LOCATION}")
add_definitions(-DLOC_LIB_NAME=${LOC_LIB_NAME})
add_definitions(-DROBOT_LOCATION=${ROBOT_LOCATION})
add_definitions(-DROBOT=${ROBOT_NAME})
if(ROBOT_NAME STREQUAL "HyQ2Max" )
  add_definitions(-DUSEHYQ2MAX=1)
  list(APPEND COMMON_LIBS iithyq2maxcommons)
  list(APPEND COMMON_LIBS iitgenhyq2max)
  endif()
if( ROBOT_NAME STREQUAL "HyQ" )
  add_definitions(-DUSEHYQ=1)
  list(APPEND COMMON_LIBS iithyqcommons)
  list(APPEND COMMON_LIBS iitgenhyq)
endif()
if( ROBOT_NAME STREQUAL "HyQReal" )
  add_definitions(-DUSEHYQREAL=1)
  list(APPEND COMMON_LIBS iithyqrealcommons)
  list(APPEND COMMON_LIBS iitgenhyqreal)
endif()

## Include directories ipopt
if(POSTURE_OPTIM)
    list(APPEND CMAKE_MODULE_PATH ${SL_ROOT}/${ROBOT}User/cmake)
    find_package(IPOPT REQUIRED)
    include_directories(${IPOPT_INCLUDE_DIRS})#must be after find_package
    add_definitions(-DPOSTURE_OPTIM=1)
    set(LOCOMOTIONUTILS_LIB ${LOC_LIB_NAME} log4cxx ipoptEigWrap ${IPOPT_LIBRARIES})
else(POSTURE_OPTIM)
    set(LOCOMOTIONUTILS_LIB ${LOC_LIB_NAME} log4cxx)
endif(POSTURE_OPTIM)

catkin_package(
  INCLUDE_DIRS include ${LCM_INCLUDE_DIRS}
  LIBRARIES ${PROJECT_NAME}
  CATKIN_DEPENDS roscpp pluginlib)

include_directories(include ${catkin_INCLUDE_DIRS} polito)

LINK_DIRECTORIES(~/dls_ws/src/dls-distro/dls_utils/common_utilities/contact_wrench_set/thirdparty/politopix/bin/)

## Declare a cpp library
add_library(${PROJECT_NAME} 	  src/FeasibleWrenchPolytope_API.cpp
                                src/mink_glpk.cpp
                                src/VertexBasedTO.cpp
                                src/get_weights_lp.cpp
                                src/chebyshev_center.cpp
                                src/OptimizationProb.cpp
                                src/testOptim.cpp
                                src/LP/glpk_wrapper.cpp
                                src/LP/test_glpk_wrapper.cpp)

target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES} polito glpk) 

## this is in case you want to generate an executable (in this case a main() function will be necessary)
#add_executable(test_ipopt tests/test_ipopt.cpp)
#target_link_libraries(test_ipopt ${PROJECT_NAME})
#set_target_properties(test_ipopt PROPERTIES COMPILE_FLAGS "-std=c++11 -o3" RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

#add_executable(test_glpk_wrapper src/LP/test_glpk_wrapper.cpp)
#target_link_libraries(test_glpk_wrapper ${PROJECT_NAME})
#set_target_properties(test_glpk_wrapper PROPERTIES COMPILE_FLAGS "-std=c++11 -o3" RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")           
