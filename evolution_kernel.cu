#ifndef _EVOLUTION_KERNEL_H_
#define _EVOLUTION_KERNEL_H_

#include <cutil_inline.h>
#include <cutil_math.h>

#include <stdio.h>

#include "constants.h"
#include "particle.h"
#include "io.h"
#include "laser.h"

// DO NOT USE RIGHT NOW DO NOT USE RIGHT NOW DO NOT USE RIGHT NOW

// Store our force data in a texture
// [y, z] -> {fy, fz}
texture<float2, 2, cudaReadModeNormalizedFloat> force_tex;

__global__ void
evolve(particle* particles, laser* lasers, float3* random_numbers) 
{
	//float normalized_y = 0.0f;
	//float normalized_z = 0.0f;
	//float2 interp_force = tex2D(force_tex, normalized_y, normalized_z);
}

#endif

