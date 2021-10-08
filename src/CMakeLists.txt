SET( ARG_PARSER_SRCS arg_parser.cpp )

SET( SRCS tweezer_render_applet.cpp 
  physobj.cpp particle.cpp constants.cpp 
  position_log.cpp laser.cpp io.cpp 
  particle_list.cpp cuda_driver.cpp 
  driver_kernels.cu )

SET( EXPERIMENT_SRCS 
  physobj.cpp particle.cpp constants.cpp 
  position_log.cpp laser.cpp io.cpp 
  particle_list.cpp exp_cuda_driver.cpp 
  experiment_driver.cu random.cu gpu.cu )

SET( CPU_SRCS physobj.cpp particle.cpp 
  constants.cpp laser.cpp io.cpp
  particle_list.cpp 
  cpu_driver.cpp )

CUDA_INCLUDE_DIRECTORIES(
  ${SOURCE_DIR}
  ${Boost_INCLUDE_DIRS}
  ${PROJECT_ROOT}/include

  # On the cluster
  /fs/gvilsw-usersupported/data/cuda-3.0-SDK/C/common/inc
  #/opt/cudasdk_2.3/C/common/inc

  # On generic machines
  /usr/local/cuda/include

  # On our machines
  /home/rob/SoftwareStaging/CUDA-3.0-SDK/C/common/inc
  /home/spook/NVIDIA_GPU_Computing_SDK/C/common/inc
)

LINK_DIRECTORIES(
  ${Boost_LIBRARY_DIRS}
  /fs/gvilsw-usersupported/data/gcc-4.2/lib
  /fs/gvilsw-usersupported/data/gcc-4.2/lib64

  # On the cluster
  /fs/gvilsw-usersupported/data/cuda-3.0-SDK/C/lib

  # On generic machines
  /usr/local/cuda/lib64

  # On our machines
  /home/rob/SoftwareStaging/CUDA-3.0-SDK/C/lib
  /home/spook/NVIDIA_GPU_Computing_SDK/C/lib
)

INCLUDE_DIRECTORIES(
  /usr/local/include
  ${GLUT_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${PROJECT_ROOT}/include
  ${PROJECT_ROOT}/include/

  # On the cluster
  /opt/cudasdk_2.3/C/common/inc

  # On generic machines
  /usr/local/cuda/include

  # On our machines
  /home/rob/SoftwareStaging/CUDA-3.0-SDK/C/common/inc
  /home/spook/NVIDIA_GPU_Computing_SDK/C/common/inc
  )

ADD_LIBRARY( arg_parser ${ARG_PARSER_SRCS} )

LINK_LIBRARIES( 
  m 
  pthread 
  cutil_x86_64 cudart ${Boost_LIBRARIES} )

#CUDA_ADD_EXECUTABLE( tweezers ${SRCS} )
CUDA_ADD_EXECUTABLE( exp_tweezers ${EXPERIMENT_SRCS} )
ADD_EXECUTABLE( cpu_tweezers ${CPU_SRCS} )

TARGET_LINK_LIBRARIES( cpu_tweezers arg_parser )
TARGET_LINK_LIBRARIES( exp_tweezers arg_parser )
