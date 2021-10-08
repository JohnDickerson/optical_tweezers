#ifndef DRIVER_KERNELS_CU
#define DRIVER_KERNELS_CU

#include "constants.h"
#include "defines.h"
#include "particle.h"
#include "laser.h"
#include "particle_list.h"
#include <cutil_inline.h>

//////  //  //  //  //  //  //  //////
//  REMEMBER TO STANDARDIZE UNITS!  //
//////  //  //  //  //  //  //  //////

using namespace tweezers;
using namespace std;

// Texture that stores the force data (laser)


extern "C" void evolve_particle(const int& num_particles, particle_list particles, laser* lasers, int step);
extern "C" void setup_cuda_constants();
extern "C" void to_force_texture(float2* h_force_array, cudaArray* d_force_array,   cudaChannelFormatDesc& channelDesc, const int& force_tex_width, const int& force_tex_height);
extern "C" void copy_random_numbers(float3* h_random_numbers, size_t nbytes);

// Texture that stores the force data (laser)
texture<float2, 2, cudaReadModeElementType> force_tex;
// Precompute random numbers on host, store on device
__constant__ float3 cuda_random_numbers[RANDOM_COUNT];

// Store physical constants on device
__constant__ float cuda_float_constants[NUM_FLOAT_CONSTANTS];
__constant__ float3 cuda_float3_constants[NUM_FLOAT3_CONSTANTS];


// float3 vector add
__device__ float3
operator+(const float3 &a, const float3 &b)
{
	return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);
}

// float3 scalar multiply
__device__ float3
operator*(const float3 &a, const float &c)
{
	return make_float3(a.x*c, a.y*c, a.z*c);
}

// evolves the laser (position, velocity, accel)
__global__ void
evolve_laser(laser* l)
{

}


__device__ float3
apply_f_ext(float3 p_pos, laser* l)
{
	// XY forces are symmetric, distance from particle -> laser
	float rel_y_pos = sqrtf( (p_pos.x - l->pos.x)*(p_pos.x - l->pos.x) + (p_pos.y - l->pos.y)*(p_pos.y - l->pos.y));
	// Z forces are not symmetric
	float rel_z_pos = p_pos.z - l->pos.z;  
	
	//printf("ry: %f, rz: %f\n", rel_y_pos*1e6, rel_z_pos*1e6);
	
	
	// We are told to ignore laser forces outside of the grid,
	// so just update the particle with gravity and buoyancy
	if(rel_y_pos < cuda_float_constants[GRID_Y_LOW] || rel_y_pos + cuda_float_constants[GRID_SPACING] > cuda_float_constants[GRID_Y_HIGH] || 
	rel_z_pos < cuda_float_constants[GRID_Z_LOW] || rel_z_pos + cuda_float_constants[GRID_SPACING] > cuda_float_constants[GRID_Z_HIGH])
	{
		return cuda_float3_constants[LANGEVIN_F_G] + cuda_float3_constants[LANGEVIN_F_B];
	}
	
	// Index into normalized texture:  x in [low, high] -> norm_x in [0,1]
	float normalized_rel_y = __fdividef((rel_y_pos - cuda_float_constants[GRID_Y_LOW]), (cuda_float_constants[GRID_Y_HIGH] - cuda_float_constants[GRID_Y_LOW]));
	float normalized_rel_z = __fdividef((rel_z_pos - cuda_float_constants[GRID_Z_LOW]), (cuda_float_constants[GRID_Z_HIGH] - cuda_float_constants[GRID_Z_LOW]));

	
	// Textures automatically bilinearly interpolate!
	float2 interp_force = tex2D(force_tex, normalized_rel_y, normalized_rel_z);
	//printf("fy: %f, fz: %f\n", interp_force.x*1e12, interp_force.y*1e12);

	// Need to move XY forces from local (laser) coordinate system to world coordinate system
	// The calculated Y force is in the direction of the laser
	float xDif = l->pos.x - p_pos.x;
	float yDif = l->pos.y - p_pos.y;
	float mag = sqrtf(xDif*xDif + yDif*yDif);

	xDif = __fdividef(xDif, mag);
	yDif = __fdividef(yDif, mag);

	float3 fInterpVec = make_float3(-xDif*interp_force.x, -yDif*interp_force.x, interp_force.y);
	
	
	// return the F_ext force vector in the paper -- this is the sum
	// of gravity (in +Z direction), buoyancy forces (in +Z direction), and
	// the trapping force estimated above
	return cuda_float3_constants[LANGEVIN_F_G] + cuda_float3_constants[LANGEVIN_F_B] + fInterpVec;
}

__global__ void
evolve_particle_kernel(particle_list particles, laser* lasers, int step,int num_particles) 
{
	int p_idx = blockDim.x * blockIdx.x + threadIdx.x;

	if(p_idx>=num_particles)
		return;
	
	// Hack for single lasers, very few particles
	float3 p_pos = particles.pos[p_idx];
	float3 p_vel = particles.vel[p_idx];
	float3 p_acc = particles.acc[p_idx];

	laser* l = lasers;	

	// Get the new acceleration for the single particle
	// -gamma/m    * V(t)
	float c1 = -cuda_float_constants[LANGEVIN_GAMMA] * cuda_float_constants[INV_SPHERE_MASS];
	float3 first_term = p_vel * c1;

	// 1/m * sqrt(xi / timestep)    * N[0,1]
	float c2 = cuda_float_constants[INV_SPHERE_MASS] * sqrt( __fdividef(cuda_float_constants[LANGEVIN_XI_SQUARED], cuda_float_constants[TIME_STEP]) );
	float3 randf3 = cuda_random_numbers[step + p_idx]; 
	float3 second_term = randf3 * c2 * 10.0f;
	
	// sum of buoyancy, gravity, optical force
	float3 third_term = apply_f_ext(p_pos, l);
	third_term = third_term * cuda_float_constants[INV_SPHERE_MASS];
	

	// a(t+dt) = v(t+dt) - v(t) / dt = c1+c2+c3
	float3 accel = first_term + second_term + third_term;

	// Evolve particle's physics using second order Verlet velocity integrator
	p_pos = p_pos + (p_vel*cuda_float_constants[TIME_STEP]) + (p_acc*cuda_float_constants[TIME_STEP]*cuda_float_constants[TIME_STEP]*0.5f);
	p_vel = p_vel + (p_acc + accel)*0.5f*cuda_float_constants[TIME_STEP];
	p_acc = accel;

	particles.pos[p_idx] = p_pos;
	particles.vel[p_idx] = p_vel;
	particles.acc[p_idx] = p_acc;
}

__global__ void 
test_force_tex()
{
	
	float rel_y_pos = cuda_float_constants[GRID_Y_LOW];
	float rel_z_pos = cuda_float_constants[GRID_Z_LOW];
	for(; rel_y_pos < cuda_float_constants[GRID_Y_HIGH]; rel_y_pos += cuda_float_constants[GRID_SPACING])
	{
		for(; rel_z_pos < cuda_float_constants[GRID_Z_HIGH]; rel_z_pos += cuda_float_constants[GRID_SPACING])
		{
			float normalized_rel_y = __fdividef((rel_y_pos - cuda_float_constants[GRID_Y_LOW]), (cuda_float_constants[GRID_Y_HIGH] - cuda_float_constants[GRID_Y_LOW]));
		float normalized_rel_z = __fdividef((rel_z_pos - cuda_float_constants[GRID_Z_LOW]), (cuda_float_constants[GRID_Z_HIGH] - cuda_float_constants[GRID_Z_LOW]));

			float2 interp_force = tex2D(force_tex, normalized_rel_y, normalized_rel_z);

			//printf("y: %f, z: %f\n", rel_y_pos*1e6, rel_z_pos*1e6);
			//printf("fy: %f, fz: %f\n", interp_force.x*1e12, interp_force.y*1e12);
		}
	}

}

void copy_random_numbers(float3* h_random_numbers, size_t nbytes) {
  // Transfer our random numbers over to constant memory
  cudaMemcpyToSymbol( cuda_random_numbers, h_random_numbers, nbytes );
}

void to_force_texture(float2* h_force_array, cudaArray* d_force_array,   cudaChannelFormatDesc& channelDesc, const int& force_tex_width, const int& force_tex_height) {
  int force_size = force_tex_width * force_tex_height * sizeof(float2);
  cutilSafeCall( cudaMallocArray( &d_force_array, &channelDesc, force_tex_width, force_tex_height )); 
  cutilSafeCall( cudaMemcpyToArray( d_force_array, 0, 0, h_force_array, force_size, cudaMemcpyHostToDevice));

  force_tex.addressMode[0] = cudaAddressModeClamp;
  force_tex.addressMode[1] = cudaAddressModeClamp;
  force_tex.filterMode = cudaFilterModeLinear;
  force_tex.normalized = true;    // texture access normalized to [0,1]

  cutilSafeCall( cudaBindTextureToArray(force_tex, d_force_array, channelDesc));
}


void evolve_particle(const int& num_particles, particle_list particles, laser* lasers, int step) {
  // Determine block/grid sizes based on num particles
  // Hacked for very few particles right now
  dim3 dimBlock(   1, 1, 1  );
  dim3 dimGrid(   num_particles, 1,  1 );

  evolve_particle_kernel<<<dimGrid, dimBlock>>>(particles, lasers, step,num_particles);
}


void setup_cuda_constants()
{
	float host_float_constants[NUM_FLOAT_CONSTANTS];
	host_float_constants[TIME_STEP] = constants::time_step;
	host_float_constants[LANGEVIN_XI_SQUARED] = constants::langevin_xi_squared;
	host_float_constants[LANGEVIN_GAMMA] = constants::langevin_gamma;
	host_float_constants[INV_SPHERE_MASS] = constants::inv_sphere_mass;
	host_float_constants[GRID_Y_LOW] = constants::grid_Y_extent_low;
	host_float_constants[GRID_Y_HIGH] = constants::grid_Y_extent_high;
	host_float_constants[GRID_Z_LOW] = constants::grid_Z_extent_low;
	host_float_constants[GRID_Z_HIGH] = constants::grid_Z_extent_high;
	host_float_constants[GRID_SPACING] = constants::grid_spacing;

	cudaMemcpyToSymbol(cuda_float_constants, host_float_constants, sizeof(host_float_constants) );

	float3 host_float3_constants[NUM_FLOAT3_CONSTANTS];
	host_float3_constants[LANGEVIN_F_G] = constants::langevin_F_g;
	host_float3_constants[LANGEVIN_F_B] = constants::langevin_F_b;
	
	cudaMemcpyToSymbol(cuda_float3_constants, host_float3_constants, sizeof(host_float3_constants) );
}

#endif // DRIVER_KERNELS_CU
