CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

if(COMMAND cmake_policy)
 cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

Project( Tweezers )

SET(CMAKE_VERBOSE_MAKEFILE ON)
SET( CMAKE_CXX_FLAGS "-g")
SET( CMAKE_C_FLAGS "" )
SET( CMAKE_BUILD_TYPE "Release" )
#SET( CUDA_BUILD_TYPE "Emulation")#Device" )
#SET( CUDA_BUILD_EMULATION ON)
#SET( CUDA_INSTALL_PREFIX "C:/CUDA" )

# Load cuda macros
INCLUDE(${CMAKE_SOURCE_DIR}/CMake/cuda/FindCUDA.cmake)
SET( CUDA_NVCC_FLAGS "--ptxas-options=-v;-use_fast_math;-arch;sm_13")#;-code;compute_11")#;-code;sm_13")


MESSAGE( STATUS "CUDA_NVCC_FLAGS=" ${CUDA_NVCC_FLAGS} )
SET( BOOST_ROOT $ENV{BOOST_ROOT} )
#FIND_PACKAGE( REQUIRED)
find_package( Boost 1.40 COMPONENTS program_options )
FIND_PACKAGE(GLUT REQUIRED)

SET( PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR} )

SET( EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin )

execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
  ${CMAKE_BINARY_DIR}/../font
  ${EXECUTABLE_OUTPUT_PATH}/font)

execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
  ${CMAKE_BINARY_DIR}/../textures
  ${EXECUTABLE_OUTPUT_PATH}/textures)

execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
  ${CMAKE_BINARY_DIR}/../force_5_microns_silica.txt
  ${EXECUTABLE_OUTPUT_PATH}/force_5_microns_silica.txt)

ADD_SUBDIRECTORY( src )

