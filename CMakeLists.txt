cmake_minimum_required (VERSION 2.6)
if (COMMAND cmake_policy)
  cmake_policy (SET CMP0003 NEW)
  cmake_policy (SET CMP0005 NEW)
  if (POLICY CMP0011)
    cmake_policy (SET CMP0011 NEW)
  endif (POLICY CMP0011)
endif (COMMAND cmake_policy)

project (kinematic_objectives)

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
  include
  )
list (APPEND
  CMAKE_REQUIRED_INCLUDES
  /usr/include/eigen3
  )

if (EIGEN3_DIR)
  list (APPEND CMAKE_REQUIRED_INCLUDES
    ${EIGEN3_DIR}/include
    ${EIGEN3_DIR}
    ${EIGEN3_DIR}/include/eigen3
    ${EIGEN3_DIR}/eigen3
    )
  include_directories (
    ${EIGEN3_DIR}/include
    ${EIGEN3_DIR}
    ${EIGEN3_DIR}/include/eigen3
    ${EIGEN3_DIR}/
    eigen3
    )
endif (EIGEN3_DIR)

include (CheckIncludeFileCXX)

check_include_file_cxx (Eigen/Core HAVE_EIGEN3)
if (HAVE_EIGEN3)
  message ("*** found Eigen3 (required)")
else (HAVE_EIGEN3)
  message (FATAL_ERROR "Please install Eigen3 or set the EIGEN3_DIR variable to tell me where to find it.")
endif (HAVE_EIGEN3)

list (APPEND SRCS
  src/util.cpp
  src/objective.cpp
  src/compound_objective.cpp
  src/constraint_teleporting_blender.cpp
  src/unconstrained_blender.cpp
  src/constraint_bouncing_blender.cpp
  src/prioritization.cpp
  src/joint_limit_objective.cpp
  src/obstacle_objective.cpp
  src/point_attraction_objective.cpp
  src/point_repulsion_objective.cpp
  src/joint_position_objective.cpp
  src/pseudo_inverse.cpp
  src/print.cpp
  src/achievability.cpp
  src/integrator.cpp
  src/blender.cpp
  )

add_library (kinematic_objectives SHARED ${SRCS})

add_subdirectory (demo)
