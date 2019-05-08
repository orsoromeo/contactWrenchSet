cmake_minimum_required(VERSION 2.8.3)
project(FWSOptim)

find_package(dwl REQUIRED)
find_package(octomap REQUIRED)

# Setting flags for optimization
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")

include_directories(include ${dwl_INCLUDE_DIRS} ${OCTOMAP_INCLUDE_DIRS}) 

#add_executable(testOptim testOptim.cpp src/OptimizationProblem.cpp)
#target_link_libraries(testOptim ${dwl_LIBRARIES} ${OCTOMAP_LIBRARY_DIRS})