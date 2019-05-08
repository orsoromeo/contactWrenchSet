# cd build
# cmake -DCMAKE_BUILD_TYPE=Release .. (or)
# cmake -DCMAKE_BUILD_TYPE=Debug  ..
# make clean && make -j VERBOSE=1     (or)
# make VERBOSE=1
cmake_minimum_required(VERSION 2.8.6)

project(politopix)

# Send all libraries and executables in the bin directory
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

message(STATUS "### begin compilers ###")
if (MSVC)
  message(STATUS "MSVC!!!!!!!!")
  if (CMAKE_VERSION VERSION_GREATER 2.8.12 OR CMAKE_VERSION VERSION_EQUAL 2.8.12)
    message(STATUS "CMake version >= 2.8.12")
    #CSF_CXX_SysOptions
    add_compile_options(/nologo)
    add_compile_options(/arch:SSE2)
    add_compile_options(/W4)
    #CSF_CXX_Options
    add_compile_options(/DWNT)
    add_compile_options(/D_WINDOWS)
    add_compile_options(/D_CRT_SECURE_NO_DEPRECATE)
    if (CMAKE_BUILD_TYPE MATCHES Release)
      #CSF_CXX_ModeOpt
      add_compile_options(/DCSFDB)
      add_compile_options(/MD)
      add_compile_options(/EHa)
      add_compile_options(/O2)
      add_compile_options(/Gy)
      add_compile_options(/DNDEBUG)
      add_compile_options(/DNo_Exception)
      message(STATUS "Release build with CMake version >= 2.8.12")
    elseif (CMAKE_BUILD_TYPE MATCHES Debug)
      #CSF_CXX_ModeDeb
      add_compile_options(/MDd)
      add_compile_options(/EHa)
      add_compile_options(/Od)
      add_compile_options(/Gm)
      add_compile_options(/Zi)
      add_compile_options(/DDEB)
      add_compile_options(/D_DEBUG)
      message(STATUS "Debug build with CMake version >= 2.8.12")
    endif()
  else()
    # V < 2.8.12
    message(STATUS "CMake version < 2.8.12")
    #CSF_CXX_SysOptions
    add_definitions(/nologo)
    add_definitions(/arch:SSE2)
    add_definitions(/W4)
    #CSF_CXX_Options
    add_definitions(/DWNT)
    add_definitions(/D_WINDOWS)
    add_definitions(/D_CRT_SECURE_NO_DEPRECATE)
    if (CMAKE_BUILD_TYPE MATCHES Debug)
      #CSF_CXX_ModeDeb
      add_definitions(/MDd)
      add_definitions(/EHa)
      add_definitions(/Od)
      add_definitions(/Gm)
      add_definitions(/Zi)
      add_definitions(/DDEB)
      add_definitions(/D_DEBUG)
      add_definitions(/DNDEBUG)
      add_definitions(/DBOOST_UBLAS_NDEBUG)
      message(STATUS "Debug build with CMake version < 2.8.12")
    else()
      #if (CMAKE_BUILD_TYPE MATCHES Release)
      #CSF_CXX_ModeOpt
      add_definitions(/DCSFDB)
      add_definitions(/MD)
      add_definitions(/EHa)
      add_definitions(/O2)
      add_definitions(/Gy)
      add_definitions(/DNDEBUG)
      add_definitions(/DBOOST_UBLAS_NDEBUG)
      add_definitions(/DNo_Exception)
      message(STATUS "Release build with CMake version < 2.8.12")
    endif()
  endif()
  #  @set %CSF_LINKER_SysShared = "-subsystem:windows -dll ";
  #  @set %CSF_LINKER_SysExe  = "CSF_LINKER_SysExeC";
  #  @set %CSF_LINKER_SysExeC = "-subsystem:CONSOLE";
  #  @set %CSF_LINKER_SysExeW = "-subsystem:WINDOWS";
  #  @string %CSF_LINKER_SysOpt = "-nologo -incremental:no" %CSF_VC_MACHINE " -debug ";
  set(CMAKE_SHARED_LINKER_FLAGS "-subsystem:windows -dll -subsystem:CONSOLE -subsystem:WINDOWS -nologo -incremental:no")
elseif (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wswitch -Wformat -Wchar-subscripts -Wparentheses -Wmultichar -Wtrigraphs -Wpointer-arith -Wcast-align -Wreturn-type -Wshadow -Wundef -Woverloaded-virtual -Wno-unused-function -fPIC -msse -msse2 -mfpmath=sse")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -DDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -fno-strict-aliasing -DNDEBUG -DBOOST_UBLAS_NDEBUG")
endif()
message(STATUS "###  end  compilers ###")

# Recurse into the "trunk" and "test" subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory(${CMAKE_SOURCE_DIR}/trunk)
#add_subdirectory(${CMAKE_SOURCE_DIR}/test)