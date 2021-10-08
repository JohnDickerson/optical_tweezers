#ifndef __TWEEZERS_DEFINES_H__
#define __TWEEZERS_DEFINES_H__

// VisualLibrary hack
#ifndef GL_FRAMEBUFFER
	#define GL_FRAMEBUFFER 0
#endif
#ifndef GLEW_ARB_draw_instanced
	#define GLEW_ARB_draw_instanced 0
#endif


// Eventually set RANDOM_COUNT to maximum constant memory
#define RANDOM_COUNT 2048


// NUM_DOUBLE_CONSTANTS refers to the single float constants our CUDA kernels must reference
#define NUM_FLOAT_CONSTANTS 11
// These are used to index the __constant__ CUDA array of our float constants
#define TIME_STEP 0
#define LANGEVIN_XI_SQUARED 1
#define LANGEVIN_GAMMA 2
#define INV_SPHERE_MASS 3
#define GRID_Y_LOW 4
#define GRID_Y_HIGH 5
#define GRID_Z_LOW 6
#define GRID_Z_HIGH 7
#define GRID_SPACING 8
#define CONST1 9
#define CONST2 10

// NUM_DOUBLE3_CONSTANTS refers to the float3 constants our CUDA kernels must reference
#define NUM_FLOAT3_CONSTANTS 2
// These are used to index the __constant__ CUDA array of our float3 constants
#define LANGEVIN_F_G 0
#define LANGEVIN_F_B 1

// NUM_UINT2_CONSTANTS refers to the uint2 constants our CUDA
// kernels must reference
#define NUM_UINT2_CONSTANTS 2
// These are used to index the __constant__ CUDA array of our uint2 constants
#define DEV_A 0
#define DEV_C 1

#endif // defines.h

