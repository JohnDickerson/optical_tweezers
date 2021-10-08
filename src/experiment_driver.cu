#ifndef EXPERIMENT_DRIVER_CU
#define EXPERIMENT_DRIVER_CU

#include <cutil.h>
#include <cutil_inline_runtime.h>

#include "random.hpp"

#include "constants.h"
#include "defines.h"
#include "particle.h"
#include "laser.h"
#include "particle_list.h"
#include "laser_list.hpp"


// Thrust
//#include <thrust/random/linear_congruential_engine.h>
//#include <thrust/random/normal_distribution.h>


//////  //  //  //  //  //  //  //////
//  REMEMBER TO STANDARDIZE UNITS!  //
//////  //  //  //  //  //  //  //////

using namespace tweezers;
using namespace std;

// Texture that stores the force data (laser)

extern "C" void evolve_particle(dim3& dimBlock, dim3& dimGrid,
                                uint2* states,
                                const int& num_particles, 
                                particle_list particles, 
                                laser_list lasers, int total_steps);

extern "C" void evolve_particle_compute_bounds(dim3& dimBlock, dim3& dimGrid,
                                               uint2* states,
                                               const int& num_particles, 
                                               particle_list particles, 
                                laser_list lasers, int total_steps);

extern "C" void setup_cuda_constants(uint2 A, uint2 C);
extern "C" void to_force_texture(float2* h_force_array, cudaArray* d_force_array,   
                                 cudaChannelFormatDesc& channelDesc, const int& force_tex_width, 
                                 const int& force_tex_height);

//extern "C" void copy_random_numbers(float3* h_random_numbers, size_t nbytes);

// Texture that stores the force data (laser)
texture<float2, 2, cudaReadModeElementType> force_tex;
// Precompute random numbers on host, store on device
//__constant__ float3 cuda_random_numbers[RANDOM_COUNT];

// Store physical constants on device
__constant__ float cuda_float_constants[NUM_FLOAT_CONSTANTS];
__constant__ float3 cuda_float3_constants[NUM_FLOAT3_CONSTANTS];
__constant__ uint2 cuda_uint2_constants[NUM_UINT2_CONSTANTS];


// float3 vector add
__device__ inline float3
operator+(const float3 &a, const float3 &b)
{
  return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);
}

// float3 vector subtract
__device__ inline float3
operator-(const float3 &a, const float3 &b)
{
  return make_float3(a.x-b.x, a.y-b.y, a.z-b.z);
}

// float3 scalar multiply
__device__ inline float3
operator*(const float3 &a, const float &c)
{
  return make_float3(a.x*c, a.y*c, a.z*c);
}

// evolves the laser (position, velocity, accel)
__device__ inline void
evolve_laser(float3& l_pos, const float3& l_vel)
{
  l_pos = l_pos + l_vel;
}


__device__ float3
apply_f_ext(float3 p_pos, float3 l_pos, bool& active)
{
	// XY forces are symmetric, distance from particle -> laser
	//float rel_y_pos = sqrtf( (p_pos.x - l_pos.x)*(p_pos.x - l_pos.x) + (p_pos.y - l_pos.y)*(p_pos.y - l_pos.y));
	float tempA=(p_pos.x - l_pos.x);
	float tempB=(p_pos.y - l_pos.y);
	float rel_y_pos = sqrtf(tempA*tempA+tempB*tempB);
	
	// Z forces are not symmetric
	float rel_z_pos = p_pos.z - l_pos.z;  
	
	// <HACK> Ashis does this
	if (rel_z_pos == 0.0f) 
	{
		rel_z_pos = 5e-4 * 1e-6;
	}
	// </HACK>


  // We are told to ignore laser forces outside of the grid,
  // so just update the particle with gravity and buoyancy

  if(	rel_y_pos < cuda_float_constants[GRID_Y_LOW] || rel_y_pos > cuda_float_constants[GRID_Y_HIGH] || 
        rel_z_pos < cuda_float_constants[GRID_Z_LOW] || rel_z_pos > cuda_float_constants[GRID_Z_HIGH] ||
        !active ) {
      active = false;
      return make_float3(0.0f, 0.0f, 
                         cuda_float3_constants[LANGEVIN_F_G].z - cuda_float3_constants[LANGEVIN_F_B].z  * 9.81);	
	}
	
  // Index into normalized texture:  x in [low, high] -> norm_x in [0,1]
  float normalized_rel_y = __fdividef((rel_y_pos - cuda_float_constants[GRID_Y_LOW]), (cuda_float_constants[GRID_Y_HIGH] - cuda_float_constants[GRID_Y_LOW]));
  float normalized_rel_z = __fdividef((rel_z_pos - cuda_float_constants[GRID_Z_LOW]), (cuda_float_constants[GRID_Z_HIGH] - cuda_float_constants[GRID_Z_LOW]));
  
  /*
  float2 tl = tex2D(force_tex, 0.0, 0.0);
  float2 bl = tex2D(force_tex, 0.0, 1.0);
  float2 tr = tex2D(force_tex, 1.0, 0.0);
  float2 br = tex2D(force_tex, 1.0, 1.0);
  */
  //printf("top left: %f %f , bottom left: %f %f, top right: %f %f, bottom right: %f %f\n",
  //      tl.x, tl.y, bl.x, bl.y, tr.x, tr.y, br.x, br.y);

  // Textures automatically bilinearly interpolate!
  float2 interp_force = tex2D(force_tex, normalized_rel_y, normalized_rel_z);
  float2 scale_interp_force = make_float2(interp_force.x * 1e-12, interp_force.y * 1e-12); 

  //printf("interp force: %.20e %.20e %.20e %.20e\n\n", 

  /*
  printf("normalized_rel_y: %f, normalized_rel_z: %f\n",
         normalized_rel_y, normalized_rel_z);
  printf("fy: %f, fz: %f\n", dinterp_force.x*1e12, dinterp_force.y*1e12);
  */

  /*
  // Need to move XY forces from local (laser) coordinate system to world coordinate system
  // The calculated Y force is in the direction of the laser
  float mag = sqrtf(xDif*xDif + yDif*yDif);

  if (rel_y_pos > 1e-20) { 
    float inv_mag = __fdividef(1.0f, mag);
    xDif *= inv_mag;
    yDif *= inv_mag;
  } else {
    xDif = 0.0f;
    yDif = 0.0f;
  }

  float3 fInterpVec = make_float3(xDif*scale_interp_force.x, yDif*scale_interp_force.x, scale_interp_force.y);
  */
	
  // From Ashis' code
  float xDif = p_pos.x - l_pos.x;
  float yDif = p_pos.y - l_pos.y;
  float angle = atan2f(yDif, xDif);
  float xForce = scale_interp_force.x * cosf(angle);
  float yForce = scale_interp_force.x * sinf(angle);
  float zForce = scale_interp_force.y;

  if( fabsf(sqrtf(xDif*xDif + yDif*yDif)) < 1e-18 )
  {	xForce = 0; yForce = 0; }
  if( fabsf(scale_interp_force.y) < 1e-17 )
  { zForce = 0; }
  float3 fInterpVec = make_float3(xForce, yForce, zForce);


  // return the F_ext force vector in the paper -- this is the sum
  // of gravity (in +Z direction), buoyancy forces (in -Z direction), and
  // the trapping force estimated above
  active = true;
  fInterpVec.z = fInterpVec.z + (cuda_float3_constants[LANGEVIN_F_G].z - cuda_float3_constants[LANGEVIN_F_B].z) * 9.81;	
  return fInterpVec;
}


/** propagate an rand48 RNG one iteration.
    @param Xn  the current RNG state, in 2x 24-bit formg
    @param A,C the magic constants for the RNG. For striding,
    this constants have to be adapted, see the constructor
    @result    the new RNG state X(n+1)
*/
__device__ uint2 RNG_rand48_iterate_single(uint2 Xn)//, uint2 A, uint2 C)
{
  // results and Xn are 2x 24bit to handle overflows optimally, i.e.
  // in one operation.

  // the multiplication commands however give the low and hi 32 bit,
  // which have to be converted as follows:
  // 48bit in bytes = ABCD EF (space marks 32bit boundary)
  // R0             = ABC
  // R1             =    D EF

  unsigned int R0, R1;
  uint2 A = cuda_uint2_constants[DEV_A];
  uint2 C = cuda_uint2_constants[DEV_C];

  // low 24-bit multiplication
  const unsigned int lo00 = __umul24(Xn.x, A.x);
  const unsigned int hi00 = __umulhi(Xn.x, A.x);

  // 24bit distribution of 32bit multiplication results
  R0 = (lo00 & 0xFFFFFF);
  R1 = (lo00 >> 24) | (hi00 << 8);

  R0 += C.x; R1 += C.y;

  // transfer overflows
  R1 += (R0 >> 24);
  R0 &= 0xFFFFFF;

  // cross-terms, low/hi 24-bit multiplication
  R1 += __umul24(Xn.y, A.x);
  R1 += __umul24(Xn.x, A.y);

  R1 &= 0xFFFFFF;

  return make_uint2(R0, R1);
}

/************************************************
 * Box-Muller transform from Uniform [-1,1] to Normal [0,1]
 ************************************************/
__device__ void BoxMuller(float& u1, float& u2)
{
  u1 = max(u1, 1e-20);
  u2 = max(u2, 1e-20);
  float   r = sqrtf(-2.0f * logf(u1));
  float phi = 2 * M_PI * u2;
  u1 = r * __cosf(phi);
  u2 = r * __sinf(phi);
}


/** create a set of random numbers. The random numbers are generated in blocks.
    In each block, a thread calculates one random number, the first thread the
    first one, the second the second and so on.
    @param state      the current states of the RNGS, one per thread.
    @param res        where to put the generated numbers
    @param num_blocks how many random numbers each thread generates.
    The total number of random numbers is therefore
    num_blocks*nThreads.
    @param A,C        the magic constants for the iteration. They need
    to be chosen as to advance the RNG by nThreads iterations
    at once, see the constructor.
*/
__device__ void RNG_rand48_get_float(uint2& lstate, float4& t) {
  //const int nThreads = blockDim.x*gridDim.x;
  const float inv_rand_max = __fdividef(1.0f, 2147483647.0f);
  // Assumes positive even num_blocks
  t.x = float(( lstate.x >> 17 ) | ( lstate.y << 7));
  t.x *= inv_rand_max;
  lstate = RNG_rand48_iterate_single(lstate);

  t.y = float(( lstate.x >> 17 ) | ( lstate.y << 7));
  t.y *= inv_rand_max;
  lstate = RNG_rand48_iterate_single(lstate);

  t.z = float(( lstate.x >> 17 ) | ( lstate.y << 7));
  t.z *= inv_rand_max;
  lstate = RNG_rand48_iterate_single(lstate);

  t.w = float(( lstate.x >> 17 ) | ( lstate.y << 7));
  t.w *= inv_rand_max;
  lstate = RNG_rand48_iterate_single(lstate);

  BoxMuller(t.x,t.y);
  BoxMuller(t.z,t.w);
  //  printf("rand = %f, %f, %f, %f\n", t.x, t.y, t.z, t.w);
}

__global__ __device__ void
evolve_particle_kernel(uint2* states, 
                       particle_list particles, 
                       laser_list lasers, 
                       int total_steps) {
  /*
	int p_idx = blockDim.x * blockIdx.x + threadIdx.x;
  */
  int p_idx =
    ((gridDim.x * blockIdx.y) + blockIdx.x) * (blockDim.x*blockDim.y)  +
    (blockDim.x * threadIdx.y) + threadIdx.x ;
  
  if (p_idx < *(particles.num_particles)) {

    float3 p_pos = particles.pos[p_idx];
    float3 p_vel = particles.vel[p_idx];
    float3 p_acc = particles.acc[p_idx];
    float3 l_pos = lasers.pos[p_idx];
    float3 l_vel = lasers.vel[p_idx];
    uint2 lstate = states[p_idx];
    bool active = true;


    //float c1 = -cuda_float_constants[LANGEVIN_GAMMA] * cuda_float_constants[INV_SPHERE_MASS];
    //float c2 = cuda_float_constants[INV_SPHERE_MASS] * 
    //  sqrt( cuda_float_constants[LANGEVIN_XI_SQUARED] / cuda_float_constants[TIME_STEP] );

	float const_CONST1=cuda_float_constants[CONST1];
	float const_CONST2=cuda_float_constants[CONST2];
	float const_INV_SPHERE_MASS=cuda_float_constants[INV_SPHERE_MASS];
	float const_TIME_STEP=cuda_float_constants[TIME_STEP];
	float const_TIME_STEP_HALF=const_TIME_STEP*.5;
	float const_TIME_STEP_SQ_HALF=const_TIME_STEP*const_TIME_STEP*.5;
	

    int step = 0;
    while( step < total_steps && active ) 
	{
		//for (int step = 0; step < total_steps; ++step) {      
		// -gamma/m    * V(t)
		float3 first_term = p_vel * const_CONST1;

		// Generate the next 4 random numbers into res
		float4 res;
		RNG_rand48_get_float(lstate, res);         

		// 1/m * sqrt(xi / timestep)    * N[0,1]
		float3 second_term = make_float3(res.x, res.y, res.z) * const_CONST2;

		// sum of buoyancy, gravity, optical force
		float3 third_term = apply_f_ext(p_pos, l_pos, active);
		third_term = third_term * const_INV_SPHERE_MASS;

		// a(t+dt) = v(t+dt) - v(t) / dt = c1+c2+c3
		float3 next_acc = first_term + second_term + third_term;

		// Evolve particle's physics using second order Verlet velocity integrator
		/*
		particles.pos[p_idx] = particles.pos[p_idx] +
		( particles.vel[p_idx] * cuda_float_constants[TIME_STEP] ) + 
		( particles.acc[p_idx] * cuda_float_constants[TIME_STEP] * cuda_float_constants[TIME_STEP] * 0.5f );

		particles.vel[p_idx] = particles.vel[p_idx] + 
		((particles.acc[p_idx] + next_acc) * ( 0.5f * cuda_float_constants[TIME_STEP] ));

		particles.acc[p_idx] = next_acc;
		*/

	  /*p_pos = p_pos + 
        ( p_vel * cuda_float_constants[TIME_STEP] ) + 
        ( p_acc * cuda_float_constants[TIME_STEP] * cuda_float_constants[TIME_STEP] * 0.5f );
      p_vel = p_vel + ((p_acc + next_acc) * ( 0.5f * cuda_float_constants[TIME_STEP] ));*/

		p_pos = p_pos + ( p_vel * const_TIME_STEP ) + ( p_acc * const_TIME_STEP_SQ_HALF);
		p_vel = p_vel + ((p_acc + next_acc) *const_TIME_STEP_HALF);
		p_acc = next_acc;

		evolve_laser(l_pos, l_vel);
		++step;
    }
    
    particles.pos[p_idx] = p_pos;
    particles.vel[p_idx] = p_vel;
    particles.acc[p_idx] = p_acc;
    particles.active[p_idx] = active;

    lasers.pos[p_idx] = l_pos;
    //lasers.vel[p_idx] = l_vel;    
  }
}


__global__ void
evolve_particle_compute_bounds_kernel(uint2* states, 
                                      particle_list particles, 
                                      laser_list lasers, 
                                      int total_steps) {
  /*
	int p_idx = blockDim.x * blockIdx.x + threadIdx.x;
  */
  int p_idx =
    ((gridDim.x * blockIdx.y) + blockIdx.x) * (blockDim.x*blockDim.y)  +
    (blockDim.x * threadIdx.y) + threadIdx.x ;
  
  if (p_idx < *(particles.num_particles)) {

    float3 p_pos = particles.pos[p_idx];
    float3 p_vel = particles.vel[p_idx];
    float3 p_acc = particles.acc[p_idx];
    uint2 lstate = states[p_idx];
    bool active = true;

    // if testing for trapping conditions, store max trapping distances
    float2 max_dists = make_float2(0.0f,0.0f);

    // Hack for single lasers
    float3 l_pos = lasers.pos[p_idx];
    float3 l_vel = lasers.vel[p_idx];
    float c1 = -cuda_float_constants[LANGEVIN_GAMMA] * cuda_float_constants[INV_SPHERE_MASS];
    float c2 = cuda_float_constants[INV_SPHERE_MASS] * 
      sqrt( cuda_float_constants[LANGEVIN_XI_SQUARED] / cuda_float_constants[TIME_STEP] );

    for (int step = 0; step < total_steps; ++step)
      {      
        // -gamma/m    * V(t)
        float3 first_term = p_vel * c1;

        // Generate the next 4 random numbers into res
        float4 res;
        RNG_rand48_get_float(lstate, res);         

        // 1/m * sqrt(xi / timestep)    * N[0,1]
        float3 second_term = make_float3(res.x, res.y, res.z) * c2;

        // sum of buoyancy, gravity, optical force
        float3 third_term = apply_f_ext(p_pos, l_pos, active);
        third_term = third_term * cuda_float_constants[INV_SPHERE_MASS];

        // a(t+dt) = v(t+dt) - v(t) / dt = c1+c2+c3
        float3 next_acc = first_term + second_term + third_term;

        // Evolve particle's physics using second order Verlet velocity integrator
        p_pos = p_pos + 
          ( p_vel * cuda_float_constants[TIME_STEP] ) + 
          ( p_acc * cuda_float_constants[TIME_STEP] * cuda_float_constants[TIME_STEP] * 0.5f );
        p_vel = p_vel + ((p_acc + next_acc) * ( 0.5f * cuda_float_constants[TIME_STEP] ));
        p_acc = next_acc;

        // See if our trapped particle has strayed farther from the optical
        // trap than it has before.  If so, store radial and z delta
        float delta_y = sqrtf( (p_pos.x - l_pos.x)*(p_pos.x - l_pos.x) + (p_pos.y - l_pos.y)*(p_pos.y - l_pos.y));
        float delta_z = p_pos.z - l_pos.z;  

        evolve_laser(l_pos, l_vel);

        // delta_y is guaranteed >= 0
        if(delta_y > max_dists.x)
          max_dists.x = delta_y;
        // we only care about z below the laser; however, z-axis
        // is flipped so we keep the greater than
        if(delta_z > max_dists.y)
          max_dists.y = delta_z;
        if (!active) break;		
      }



    particles.pos[p_idx] = p_pos;
    particles.vel[p_idx] = p_vel;
    particles.acc[p_idx] = p_acc;
    
    if(!active) {
      particles.active[p_idx] = false;
    } 

    // If we are generating trapping bounds, overwrite the particle's
    // calculated position with the max_radial and max_z distances
    particles.pos[p_idx].x = 0.0f;			
    particles.pos[p_idx].y = max_dists.x;
    particles.pos[p_idx].z = max_dists.y;
  }
}



__global__ void
find_trapping_bounds(uint2* states, 
                       int num_particles, 
                       particle_list particles, 
                       laser_list lasers, 
                       int total_steps) {
  /*
	int p_idx = blockDim.x * blockIdx.x + threadIdx.x;
  */
  int p_idx =
    ((gridDim.x * blockIdx.y) + blockIdx.x) * (blockDim.x*blockDim.y)  +
    (blockDim.x * threadIdx.y) + threadIdx.x ;
  
  if (p_idx < num_particles) {

    float3 p_pos = particles.pos[p_idx];
    float3 p_vel = particles.vel[p_idx];
    float3 p_acc = particles.acc[p_idx];
    uint2 lstate = states[p_idx];
    bool active = true;

    // Hack for single lasers
    float3 l_pos = lasers.pos[p_idx];
    float3 l_vel = lasers.vel[p_idx];

    float c1 = -cuda_float_constants[LANGEVIN_GAMMA] * cuda_float_constants[INV_SPHERE_MASS];
    float c2 = cuda_float_constants[INV_SPHERE_MASS] * 
      sqrt( cuda_float_constants[LANGEVIN_XI_SQUARED] / cuda_float_constants[TIME_STEP] );

    for (int step = 0; step < total_steps; ++step)
		{      
      // -gamma/m    * V(t)
      float3 first_term = p_vel * c1;

			// Generate the next 4 random numbers into res
      float4 res;
      RNG_rand48_get_float(lstate, res);         

      // 1/m * sqrt(xi / timestep)    * N[0,1]
      float3 second_term = make_float3(res.x, res.y, res.z) * c2;

      // sum of buoyancy, gravity, optical force
      float3 third_term = apply_f_ext(p_pos, l_pos, active);
      third_term = third_term * cuda_float_constants[INV_SPHERE_MASS];

      // a(t+dt) = v(t+dt) - v(t) / dt = c1+c2+c3
      float3 next_acc = first_term + second_term + third_term;

      // Evolve particle's physics using second order Verlet velocity integrator
      p_pos = p_pos + 
        ( p_vel * cuda_float_constants[TIME_STEP] ) + 
        ( p_acc * cuda_float_constants[TIME_STEP] * cuda_float_constants[TIME_STEP] * 0.5f );
      p_vel = p_vel + ((p_acc + next_acc) * ( 0.5f * cuda_float_constants[TIME_STEP] ));
			p_acc = next_acc;

      evolve_laser(l_pos, l_vel);
      if (!active) break;		
    }
    
    particles.pos[p_idx] = p_pos;
		particles.vel[p_idx] = p_vel;
		particles.acc[p_idx] = p_acc;

    if (!active) {
      particles.active[p_idx] = false;
    } 

  }
}




__global__ void 
test_force_tex()
{
  for(float rel_y_pos = cuda_float_constants[GRID_Y_LOW];// + 0.125e-6; 
			rel_y_pos < cuda_float_constants[GRID_Y_HIGH]; 
			rel_y_pos += cuda_float_constants[GRID_SPACING])
	{
      for(float rel_z_pos = cuda_float_constants[GRID_Z_LOW];// +0.125e-6; 
					rel_z_pos < cuda_float_constants[GRID_Z_HIGH]; 
					rel_z_pos += cuda_float_constants[GRID_SPACING])
		{
           float normalized_rel_y = __fdividef((rel_y_pos - cuda_float_constants[GRID_Y_LOW]), (cuda_float_constants[GRID_Y_HIGH] - cuda_float_constants[GRID_Y_LOW]));
  float normalized_rel_z = __fdividef((rel_z_pos - cuda_float_constants[GRID_Z_LOW]), (cuda_float_constants[GRID_Z_HIGH] - cuda_float_constants[GRID_Z_LOW]));

          float2 interp_force = tex2D(force_tex, normalized_rel_y, normalized_rel_z);

					//printf("%.10f %.10f %.10f %.10f\n", rel_y_pos * 1e6, rel_z_pos * 1e6, interp_force.x, interp_force.y);
		
		}
	}

}

/*
  void copy_random_numbers(float3* h_random_numbers, size_t nbytes) {
  // Transfer our random numbers over to constant memory
  cudaMemcpyToSymbol( cuda_random_numbers, h_random_numbers, nbytes );
  }
*/
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

void evolve_particle(dim3& dimBlock,
                     dim3& dimGrid,
                     uint2* states,
                     const int& num_particles, 
                     particle_list particles, 
                     laser_list lasers, 
                     int total_steps) {

  evolve_particle_kernel<<<dimGrid, dimBlock>>>(states, particles, lasers, total_steps);
  cudaError_t error = cudaGetLastError();
  if(error!=cudaSuccess) {
    cerr << cudaGetErrorString(error) << "\n";
    exit(-1);
  }         
	//dim3 ones(1,1,1);
	//test_force_tex<<<ones,ones>>>();
}


void evolve_particle_compute_bounds(dim3& dimBlock,
                                    dim3& dimGrid,
                                    uint2* states,
                                    const int& num_particles, 
                                    particle_list particles, 
                                    laser_list lasers, 
                                    int total_steps) {
        
  evolve_particle_compute_bounds_kernel<<<dimGrid, dimBlock>>>(states,particles, lasers, total_steps);
  cudaError_t error = cudaGetLastError();
  if(error!=cudaSuccess) {
    cerr << cudaGetErrorString(error) << "\n";
    exit(-1);
  }         
  //dim3 ones(1,1,1);
  //test_force_tex<<<ones,ones>>>();
}

void setup_cuda_constants(uint2 A, uint2 C)
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
  host_float_constants[CONST1] = constants::c1;
  host_float_constants[CONST2] = constants::c2;
  cerr << "C2 = " << host_float_constants[CONST2] << "\n";
  cudaMemcpyToSymbol(cuda_float_constants, host_float_constants, sizeof(host_float_constants) );

  float3 host_float3_constants[NUM_FLOAT3_CONSTANTS];
  host_float3_constants[LANGEVIN_F_G] = constants::langevin_F_g;
  host_float3_constants[LANGEVIN_F_B] = constants::langevin_F_b;
	
  cudaMemcpyToSymbol(cuda_float3_constants, host_float3_constants, sizeof(host_float3_constants) );

  uint2 host_uint2_constants[NUM_UINT2_CONSTANTS];
  host_uint2_constants[DEV_A] = A;
  host_uint2_constants[DEV_C] = C;
  cudaMemcpyToSymbol(cuda_uint2_constants, host_uint2_constants, sizeof(host_uint2_constants) );
}

#endif // EXPERIMENT_DRIVER_CU
