cmake_minimum_required (VERSION 2.6)
if (COMMAND cmake_policy)
  cmake_policy (SET CMP0003 NEW)
  cmake_policy (SET CMP0005 NEW)
  if (POLICY CMP0011)
    cmake_policy (SET CMP0011 NEW)
  endif (POLICY CMP0011)
endif (COMMAND cmake_policy)

project (kinematic_elastic)

set (CMAKE_VERBOSE_MAKEFILE ON)
if (NOT BUILD_SHARED_LIBS)
  set (BUILD_SHARED_LIBS True)
endif (NOT BUILD_SHARED_LIBS)
if (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Debug)
endif (NOT CMAKE_BUILD_TYPE)

include (CheckCXXCompilerFlag)

check_cxx_compiler_flag (-pipe CXX_FLAG_pipe)
if (CXX_FLAG_pipe)
  add_definitions (-pipe)
endif (CXX_FLAG_pipe)
check_cxx_compiler_flag (-Wall CXX_FLAG_Wall)
if (CXX_FLAG_Wall)
  add_definitions (-Wall)
endif (CXX_FLAG_Wall)

if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  check_cxx_compiler_flag (-O0 CXX_FLAG_O0)
  if (CXX_FLAG_O0)
    add_definitions (-O0)
  endif (CXX_FLAG_O0)
endif (${CMAKE_BUILD_TYPE} STREQUAL "Debug")

if (APPLE)
  include_directories (
    /opt/local/include
    /opt/local/include/eigen3
    )
  list (APPEND CMAKE_REQUIRED_INCLUDES
    /opt/local/include
    /opt/local/include/eigen3
    )
endif (APPLE)

include_directories (
  /usr/include/eigen3
  )
list (APPEND
  CMAKE_REQUIRED_INCLUDES
  /usr/include/eigen3
  )

### if (EIGEN3_DIR)
###   list (APPEND CMAKE_REQUIRED_INCLUDES
###     ${EIGEN3_DIR}/include
###     ${EIGEN3_DIR}
###     ${EIGEN3_DIR}/include/eigen3
###     ${EIGEN3_DIR}/eigen3
###     )
###   include_directories (
###     ${EIGEN3_DIR}/include
###     ${EIGEN3_DIR}
###     ${EIGEN3_DIR}/include/eigen3
###     ${EIGEN3_DIR}/
###     eigen3
###     )
### endif (EIGEN3_DIR)

include (CheckIncludeFileCXX)

check_include_file_cxx (Eigen/Core HAVE_EIGEN3)
if (HAVE_EIGEN3)
  message ("*** found Eigen3 (required)")
else (HAVE_EIGEN3)
  message (FATAL_ERROR "Please install Eigen3.")
endif (HAVE_EIGEN3)

find_package(GTK2 REQUIRED gtk)
if (GTK2_FOUND)
  message ("*** found GTK2 (required)")
  include_directories (${GTK2_INCLUDE_DIRS})
else (GTK2_FOUND)
  message (FATAL_ERROR "Please install GTK2 development files (or tweak build sys to skip GTK-dependent stuff).")
endif (GTK2_FOUND)

list (APPEND SRCS
  kinematic_elastic.cpp
  task.cpp
  model.cpp
  algorithm.cpp
  mistry_algorithm.cpp
  pseudo_inverse.cpp
  print.cpp
  )

add_library (kinematic_elastic SHARED ${SRCS})

add_executable (elastic_gtk elastic_gtk.cpp)
target_link_libraries (elastic_gtk kinematic_elastic ${GTK2_LIBRARIES})
