project( OTRegion1 )

set(CMAKE_CXX_FLAGS "-g")

set(SOURCES 
  cpuSphericalParticleSim.cpp
  main.cpp
  RNG.cpp
  velocity.cpp
  CPUCompute.hpp
  cpuSphericalParticleSim.hpp
  PhysicalConstants.hpp
  RNG.hpp
  SphericalParticle.hpp
  io.cpp
)

INCLUDE_DIRECTORIES(
  ${PROJECT_ROOT}/include
  /opt/cudasdk_2.3/C/common/inc
  /usr/local/cuda/include
  /usr/local/stow/cudasdk-2.3/C/common/inc
  /usr/local/stow/NVIDIA-CUDA-SDK-2.3/C/common/inc
  )

add_definitions( 
  -D_CONSOLE
  -D_DEBUG
)

add_executable( OTRegion1 ${SOURCES} )

