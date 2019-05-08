set(politopix_VERSION_MAJOR 4)
set(politopix_VERSION_MINOR 0)
set(politopix_VERSION_PATCH 0)
# To export your lib and build the dll from CMake Version 2.8.6 add the following line.
include(GenerateExportHeader)
# configure a header file to pass some of the CMake settings to the source code
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Config.h.in" "${CMAKE_CURRENT_SOURCE_DIR}/Config.h")
#configure_file("${PROJECT_SOURCE_DIR}/Config.h.in" "${PROJECT_BINARY_DIR}/Config.h")
# add the binary tree to the search path for include files so that we will find Config.h
#include_directories("${PROJECT_BINARY_DIR}")
find_package(Boost)
if(Boost_FOUND)
  message(STATUS "Boost found by CMake")
  include_directories(${Boost_INCLUDE_DIRS}) 
  #target_link_libraries(progname ${Boost_LIBRARIES})
elseif(BOOST_PATH)
  message(STATUS "Use of variable BOOST_PATH")
  include_directories(${BOOST_PATH}) 
elseif(MSVC)
  message(STATUS "Use of predefined variable")
  set(Boost_INCLUDE_DIR D:/boost/boost_1_55_0)
  include_directories(${Boost_INCLUDE_DIR})
endif()

set(SOURCES 
  DoubleDescription_Rn.h  GeometricObjectIterator_Rn.h  IO_Polytope.h    NormalFan_Rn.h  politopixAPI.h       Polytope_Rn.h  Tracking.h
  Generator_Rn.h          HalfSpace_Rn.h                Neighbours_Rn.h  Point_Rn.h      PolyhedralCone_Rn.h  Rn.h           VolumeOfPolytopes_Rn.h
  Generator_Rn.cpp  HalfSpace_Rn.cpp  IO_Polytope.cpp NormalFan_Rn.cpp  Point_Rn.cpp  PolyhedralCone_Rn.cpp  Polytope_Rn.cpp  VolumeOfPolytopes_Rn.cpp
  Rn.cpp  politopixAPI.cpp	PolyhedralAlgorithms_Rn.h		PolyhedralAlgorithms_Rn.cpp	     Voronoi_Rn.cpp  Voronoi_Rn.h 
)

add_library(polito SHARED ${SOURCES})
# To have the .lib created and to generate the symbols in the .dll
generate_export_header(polito
             BASE_NAME polito
             EXPORT_MACRO_NAME polito_EXPORT
             EXPORT_FILE_NAME ${CMAKE_CURRENT_SOURCE_DIR}/polito_Export.h
             STATIC_DEFINE polito_BUILT_AS_STATIC
)
# Make sure the compiler can find include files for our libpolito library
# when other libraries or executables link to libpolito
target_include_directories(polito PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

#add_executable(politopix main.cpp ${SOURCES})
add_executable(politopix main.cpp)
# Add executable called "test_politopix" that is built from the source files "test_politopix".
add_executable(test_politopix test_politopix.cpp)
target_link_libraries(politopix LINK_PUBLIC polito)
# Link the executable to the polito library. Since the politopix library has
# public include directories we will use those link directories when building test_politopix
target_link_libraries(test_politopix LINK_PUBLIC polito)
