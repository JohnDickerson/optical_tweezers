#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <map>

//#define _WINDOWS_64_
#ifdef _WINDOWS_64_
	#include <fstream>
	#include <iomanip>
	#include <conio.h>
#endif

#include <boost/random/normal_distribution.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/cstdint.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>

#include "defines.h"

#include "constants.h"
#include "particle.h"
#include "io.h"
#include "laser.h"
#include "particle_list.h"
#include "laser_list.hpp"
#include "position_log.hpp"

#include <cutil_inline.h>
#include "cutil.h"
#include "cutil_math.h"
#include "cuda_runtime.h"

// Random Number Generation
#include "gpu.hpp"
#include "arg_parser.hpp"
#include "random.hpp"
#include "util.hpp"


//////  //  //  //  //  //  //  //////
//  REMEMBER TO STANDARDIZE UNITS!  //
//////  //  //  //  //  //  //  //////

using namespace tweezers;
using namespace std;
using boost::lexical_cast;
using boost::is_any_of;
using boost::split;

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
extern "C" void to_force_texture(float2* h_force_array, cudaArray* d_force_array,   cudaChannelFormatDesc& channelDesc, const int& force_tex_width, const int& force_tex_height);
//extern "C" void copy_random_numbers(float3* h_random_numbers, size_t nbytes);

// RNG
boost::normal_distribution<float> norm_dist(0.0, 1.0);
boost::lagged_fibonacci19937 engine;
/** Namespace declarations */
namespace po = boost::program_options;

// Only works for Ashis' three examples
particle_list setup_particles(int num_particles)
{
	particle_list particles;
	particles.pos = new float3[num_particles];
	particles.vel = new float3[num_particles];
	particles.acc = new float3[num_particles];
	particles.active = new bool[num_particles];
	particles.num_particles = new int;

	return particles;
}

laser_list setup_lasers(int num_particles)
{
	laser_list lasers;
	lasers.pos = new float3[num_particles];
	lasers.vel = new float3[num_particles];
	lasers.num_lasers = new int;
	*lasers.num_lasers = num_particles;

	return lasers;
}

// Only works for a single laser right now
laser** setup_laser(int num_lasers)
{
	laser** lasers = new laser*[num_lasers];

	lasers[0] = new laser(0,0,0);

	return lasers;
}

void generateRandomNumbers(float3* h_random_numbers, int amount)
{
	for(int idx=0; idx<amount; idx++)
	{
		h_random_numbers[idx].x = norm_dist.operator () <boost::lagged_fibonacci19937>((engine));
		h_random_numbers[idx].y = norm_dist.operator () <boost::lagged_fibonacci19937>((engine));
		h_random_numbers[idx].z = norm_dist.operator () <boost::lagged_fibonacci19937>((engine));
	}
}

// Squishes our loaded force data structure into a sorted form that can be
// pushed into texture memory as a cudaArray
float2* form_force_array(io::forcemap* force_data, int width, int height)
{
	float2* force_array = new float2[width*height];
	int w_idx=0, h_idx=0;

	size_t idx = 0;
	io::forcemap::iterator it = force_data->begin();
	for(; it != force_data->end(); it++, ++idx) {
		force_array[idx].x = it->second.x;
		force_array[idx].y = it->second.y;
	}
	std::cerr << "idx = " << idx << ", size = " << width*height << "\n";

	return force_array;
}

// If the particle is within 1 micron XY and 1 micron Z of the laser trap, then
// it was successfully trapped -- rule of thumb?  Ashis' thesis is unclear on #s	
bool isTrapped(const float3& particle, const float3& laser, const float& radial_trap_limit, const float& axial_trap_limit)
{
	//cout << particle.x << " " << particle.y << " " << particle.z << endl;
	if ( abs(particle.x - laser.x) < 2.0*radial_trap_limit && 
		abs(particle.y - laser.y) < 2.0*radial_trap_limit 
		&& (particle.z - laser.z) >= 0.0 && 
		(particle.z - laser.z) < 2.0*axial_trap_limit ) {
			return true;
	} else {
		return false;
	}

}

int main(int argc, char** argv)
{
	/*  if(argc < 4)
	{
	cerr << "Usage: " << argv[0] << " NUM_PARTICLES Y_LOW Y_HIGH" << endl;
	return -1;
	}
	*/
	
	po::options_description desc("Allowed Options");
	desc.add_options()
		("help,h", "produce this message")
		("y_start,l", po::value<float>(), 
		"pdb file")
		("y_end,r", po::value<float>(), 
		"chain file")
		("z_start,t", po::value<float>(), 
		"configuration file")
		("z_end,b", po::value<float>(), 
		"curvature output file")
		("duration,d", po::value<float>(), 
		"duration of the simulation")
		("num_particles,n", po::value<int>(), 
		"# of particles")
		("moving_laser,m", po::value<bool>()->default_value(false), 
		"is the laser moving?")
		("laser_velocity,v", po::value<string>()->default_value("0.0:0.0:0.0"), 
		"What's the velocity of the meter (in m/s)")
		;

	try
	{

		#ifndef	_WINDOWS_64_
			ArgParser ap( desc );
			ap.parse_args( argc, argv );
			float total_time = ap["duration"].as<float>();
			float y_start = ap["y_start"].as<float>();
			float y_end = ap["y_end"].as<float>();
			float z_start = ap["z_start"].as<float>();
			float z_end = ap["z_end"].as<float>();
			int num_particles = ap["num_particles"].as<int>();//atoi(argv[1]);//1521;	

			string laser_velocity_s;
			bool moving_laser = ap["moving_laser"].as<bool>();
			laser_velocity_s = ap["laser_velocity"].as<string>();    

			typedef vector< string > split_vector_t;
			split_vector_t split_vec;
			split( split_vec, laser_velocity_s, is_any_of(":") ); 
			assert(split_vec.size() == 3);
			
			float3 laser_velocity = make_float3(lexical_cast<float>(split_vec[0]),
				lexical_cast<float>(split_vec[1]),
				lexical_cast<float>(split_vec[2]));
			laser_velocity *= constants::time_step;
			cerr << "laser velocity: " << util::to_s(laser_velocity) << "\n";

			// Load pre-computed force data, [Y Z ForceY ForceZ]
			#ifndef	WIN32
				string forces_path = "/chimerahomes/sujal/tweezers/force_5_microns_silica.txt";	
			#else
				string forces_path = "..\\..\\force_5_microns_silica.txt";	
			#endif
			//init device
			//cudaSetDevice( cutGetMaxGflopsDeviceId());
			cudaSetDevice(1);
		#else
			// Load pre-computed force data, [Y Z ForceY ForceZ]
			string forces_path = "..\\..\\force_5_microns_silica.txt";	

			//init device
			cudaSetDevice( cutGetMaxGflopsDeviceId());
		#endif

		io::forcemap* force_data = io::load_force_data(forces_path,constants::dist_data_scale,constants::force_data_scale);
		//io::print_loaded_force_data(force_data);

		const int seed_value = 2345;
		srand(seed_value);
		cerr << "# Seed: " << seed_value << endl;

		// Create timer for later
		unsigned int timer = 0;
		cutilCheckError( cutCreateTimer( &timer));

		#ifdef	_WINDOWS_64_

		ofstream out("..\\outputGPU.txt");
		out<<setprecision(9);
		for(unsigned y=0;y<21;y++)
		{
			float total_time = 1.0f;
			float y_start = y*1.0e-6f;
			float y_end = y_start+.99e-6f;
			float z_start = -20.0e-6f;
			float z_end = 8.1e-6f;
			size_t num_particles = 100;	

			bool moving_laser=false;
			float3 laser_velocity=make_float3(0,0,0);
			
			out<<"\n#"<<y<<"\t";
			out<<"-l "<<y_start<<" ";
			out<<"-r "<<y_end<<" ";
			out<<"-t "<<z_start<<" ";
			out<<"-b "<<z_end<<" ";
			out<<"-d "<<total_time<<" ";
			out<<"-n "<<num_particles<<"\n";

		#endif


			//double y_start = static_cast<float>(atoi(argv[2])) * 1e-6;//constants::grid_Y_extent_low;
			//double y_end = static_cast<float>(atoi(argv[3])) * 1e-6;//constants::grid_Y_extent_high;
			//double z_start = -20e-6; // -40 to 10 microns on Z
			//double z_end = 8e-6;

			//double total_time = 1.0; // in seconds

			double print_time = 0.25e-3;  // print every .25ms
			int total_iterations = (total_time / constants::time_step) + 1;
			float steps_per_micron = constants::grid_spacing;
			double test_step = constants::grid_spacing;  // Every 1.0 microns

			size_t num_y_steps = ceil(static_cast<float>(std::abs(y_end - y_start)) / steps_per_micron);
			size_t num_z_steps = ceil(static_cast<float>(std::abs(z_end - z_start)) / steps_per_micron);

			cout << "# Max time: " << total_time << " seconds" << endl;
			cout << "# Timestep: " << constants::time_step * 1e6 << " microsec" << endl;
			cout << "# Num particles " << num_particles << endl;
			cout << "# y steps: " << num_y_steps << endl;
			cout << "# z steps: " << num_z_steps << endl;			
			cout << "# Iterations: " << total_iterations << " per particle per run" << endl;
			cout << "# y: [" << y_start << ", " << y_end << "]" << endl;
			cout << "# z: [" << z_start << ", " << z_end << "]" << endl;
			cout << "# Test step: " << test_step * 1e6 << " microns." << endl;


			//---------------- CPU-side Setup Starts Here ----------------//
			
			// Initiate some random particles
			//int num_particles = atoi(argv[1]);//1521;	
			particle_list h_particles = setup_particles(num_particles);

			// Initialize laser with beam focus at origin
			int num_lasers = num_particles;
			laser_list h_lasers = setup_lasers(num_lasers);
			//laser** h_lasers = setup_lasers(num_lasers);
			//laser_list h_lasers;

			// CPU-side create random number array for Brownian motion
			//int random_amount = (int)(total_time / constants::time_step) + 1;
			int random_amount = RANDOM_COUNT;
			float3* h_random_numbers = new float3[random_amount];
			generateRandomNumbers(h_random_numbers, random_amount);

			//---------------- Interactions with the GPU start here ----------------//
			
			// Malloc our initial particles - transfer occurs in testing loop below
			particle_list d_particles;
			cutilSafeCall( cudaMalloc( (void**) &d_particles.pos, sizeof(float3) * num_particles ));
			cutilSafeCall( cudaMalloc( (void**) &d_particles.vel, sizeof(float3) * num_particles ));
			cutilSafeCall( cudaMalloc( (void**) &d_particles.acc, sizeof(float3) * num_particles ));
			cutilSafeCall( cudaMalloc( (void**) &d_particles.active, sizeof(bool) * num_particles ));
			cutilSafeCall( cudaMalloc( (void**) &d_particles.num_particles, sizeof(int)) );

			// Malloc our initial lasers - transfer occurs in testing loop below
			laser_list d_lasers;
			cutilSafeCall( cudaMalloc( (void**) &d_lasers.pos, sizeof(float3) * num_lasers ));
			cutilSafeCall( cudaMalloc( (void**) &d_lasers.vel, sizeof(float3) * num_lasers ));
			cutilSafeCall( cudaMalloc( (void**) &d_lasers.num_lasers, sizeof(int)) );

			// Transfer our initial laser location over as a 1D array
			//laser* d_lasers = NULL;
			//cutilSafeCall( cudaMalloc( (void**) &d_lasers, sizeof(laser) * num_lasers) );
			//cutilSafeCall( cudaMemcpy( d_lasers, h_lasers, sizeof(laser) * num_lasers, cudaMemcpyHostToDevice) );

			// Transfer our force data over as a texture
			cudaArray* d_force_array=NULL;
			cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float2>();
			/*int force_tex_width = 
			round(((constants::grid_Y_extent_high - constants::grid_Y_extent_low) /  
			constants::grid_spacing)) + 1;
			int force_tex_height =
			round(((constants::grid_Z_extent_high - constants::grid_Z_extent_low) /  
			constants::grid_spacing)) + 1;*/
			int force_tex_width = 
				floor(((constants::grid_Y_extent_high - constants::grid_Y_extent_low) /  
				constants::grid_spacing)+.5f) + 1;
			int force_tex_height =
				floor(((constants::grid_Z_extent_high - constants::grid_Z_extent_low) /  
				constants::grid_spacing)+.5f) + 1;

			float2* h_force_array = form_force_array(force_data, force_tex_width, force_tex_height);

			to_force_texture(h_force_array, d_force_array, channelDesc, force_tex_width, force_tex_height);

			//  copy_random_numbers(h_random_numbers, sizeof(float3)*random_amount);

			//---------------- Evolution on the GPU Starts Here ----------------//
			//cerr << "# Setup complete, beginning evolution." << endl;

			int total_y_idx, total_z_idx, y_idx=0, z_idx=0;
			

			//double y_start = static_cast<float>(atoi(argv[2])) * 1e-6;//constants::grid_Y_extent_low;
			//double y_end = static_cast<float>(atoi(argv[3])) * 1e-6;//constants::grid_Y_extent_high;
			//double z_start = -20e-6; // -40 to 10 microns on Z
			//double z_end = 8e-6;

			total_y_idx = ceil( (y_end - y_start) / test_step );
			total_z_idx = ceil( (z_end - z_start) / test_step );	//sujal: changed to get exact timestep as CPU
			//total_z_idx = ceil( (z_end - z_start) / test_step ) + 1;
			
			unsigned int old_timer = 0;

			GPU_init();
			{
				
				RNG_rand48 r48;
				// Determine block/grid sizes based on num particles
				size_t s = min( static_cast<size_t>(ceil( sqrt((double)num_particles) )), 
					static_cast<size_t>(8));
				size_t gdx = ceil( num_particles / static_cast<float>((s*s)) );

				dim3 dimBlock(s, s, 1);
				dim3 dimGrid(gdx, 1, 1);

				const int nThreads = dimBlock.x * dimBlock.y * dimBlock.z * dimGrid.x * dimGrid.y; 
				//std::cerr << "NTHREADS : " << nThreads << "\n";
				//std::cerr << "dimBlock = " << s << " " << s << " " << 1 << "\n";
				//std::cerr << "dimGrid = " << gdx << " " << 1 << " " << 1 << "\n";

				unsigned int sA0, sA1, sC0, sC1;
				r48.init(seed_value, nThreads, sA0, sA1, sC0, sC1);

				uint2* states =  r48.get_state_ptr();
				// Transfer all the random physical constants and grid spacing stuff to constant memory
				uint2 A = make_uint2(sA0, sA1);
				uint2 C = make_uint2(sC0, sC1);

				setup_cuda_constants(A, C);


				// Setup for the trapping bounds
				//cerr << "#### CALC LIMITS #######" << endl;
				for (int i=0; i<num_particles; i++) {
					h_particles.pos[i] = make_float3(0,0,0);
					h_particles.vel[i] = make_float3(0,0,0);
					h_particles.acc[i] = make_float3(0,0,0);
					h_particles.active[i] = true;
				}
				*(h_particles.num_particles) = num_particles;

				for(int i=0; i<num_lasers; i++) {
					h_lasers.pos[i] = make_float3(0.0, 0.0, 0.0);
					h_lasers.vel[i] = laser_velocity;
				}
				*(h_lasers.num_lasers) = num_lasers;

				// Transfer particles to the device
				cutilSafeCall( cudaMemcpy( d_particles.pos, h_particles.pos, 
					sizeof(float3) * num_particles, cudaMemcpyHostToDevice) );
				cutilSafeCall( cudaMemcpy( d_particles.vel, h_particles.vel, 
					sizeof(float3) * num_particles, cudaMemcpyHostToDevice) );
				cutilSafeCall( cudaMemcpy( d_particles.acc, h_particles.acc, 
					sizeof(float3) * num_particles, cudaMemcpyHostToDevice) );
				cutilSafeCall( cudaMemcpy( d_particles.active, h_particles.active, 
					sizeof(bool) * num_particles, cudaMemcpyHostToDevice) );
				cutilSafeCall( cudaMemcpy( d_particles.num_particles, h_particles.num_particles, 
					sizeof(int), cudaMemcpyHostToDevice) );

				// Transfer lasers to the device
				cutilSafeCall( cudaMemcpy( d_lasers.pos, h_lasers.pos, 
					sizeof(float3) * num_lasers, cudaMemcpyHostToDevice) );
				cutilSafeCall( cudaMemcpy( d_lasers.vel, h_lasers.vel, 
					sizeof(float3) * num_lasers, cudaMemcpyHostToDevice) );
				cutilSafeCall( cudaMemcpy( d_lasers.num_lasers, 
					h_lasers.num_lasers, 
					sizeof(int), cudaMemcpyHostToDevice) );			

				double limits_world_time = 0.1;  // in seconds
				int limits_iterations = (limits_world_time / constants::time_step) + 1;
				evolve_particle_compute_bounds(dimBlock, dimGrid, 
					states, num_particles, 
					d_particles, d_lasers, limits_iterations);

				cutilSafeCall( cudaMemcpy( h_particles.pos, d_particles.pos, sizeof(float3) * num_particles, cudaMemcpyDeviceToHost) );
				cutilSafeCall( cudaMemcpy( h_particles.active, d_particles.active, sizeof(bool) * num_particles, cudaMemcpyDeviceToHost) );

				// We have a list of the radial_max and z_max for each of the particles
				// Now find the maximums over all radial_max, z_max

				float radial_trap_limit = 0.0f;
				float axial_trap_limit = 0.0f;
				for(size_t idx=0; idx<num_particles; idx++)
				{
					float3 pos = h_particles.pos[idx];
					if(pos.y > radial_trap_limit)
						radial_trap_limit = pos.y;
					if(pos.z > axial_trap_limit)
						axial_trap_limit = pos.z;
				}

				//radial_trap_limit = 5e-7f;
				//axial_trap_limit = 5e-7f;
				//cerr << "#Radial Trap Limit: " << radial_trap_limit << "\n";
				//cerr << "#Axial Trap Limit: " << axial_trap_limit << "\n";



				// Now that we have experimentally justified trapping bounds, we can start
				// doing the real evolution, checking whether particles are trapped or not
				//cerr << "#### EVOLVE START ######" << endl;
				for(y_idx=0; y_idx < total_y_idx; y_idx++)
				{
					for(z_idx=0; z_idx < total_z_idx; z_idx++)
					{
						old_timer = timer;

						double y_dist = y_start + ((double) y_idx * test_step);
						double z_dist = z_start + ((double) z_idx * test_step);

						//y_dist =  0.000010;
						//z_dist = -0.000040;

						// Setup the particles for this specific Y, Z
						for(int i=0; i<num_particles; i++) {
							h_particles.pos[i] = make_float3(0, y_dist, z_dist);
							h_particles.vel[i] = make_float3(0,0,0);
							h_particles.acc[i] = make_float3(0,0,0);
							h_particles.active[i] = true;
						}

						for(int i=0; i<num_lasers; i++) {
							h_lasers.pos[i] = make_float3(0.0, 0.0, 0.0);
							h_lasers.vel[i] = laser_velocity;
						}

						// Move particle positions to the device
						cutilSafeCall( cudaMemcpy( d_particles.pos, h_particles.pos, 
							sizeof(float3) * num_particles, cudaMemcpyHostToDevice) );
						cutilSafeCall( cudaMemcpy( d_particles.vel, h_particles.vel, 
							sizeof(float3) * num_particles, cudaMemcpyHostToDevice) );
						cutilSafeCall( cudaMemcpy( d_particles.acc, h_particles.acc, 
							sizeof(float3) * num_particles, cudaMemcpyHostToDevice) );
						cutilSafeCall( cudaMemcpy( d_particles.active, h_particles.active, 
							sizeof(bool) * num_particles, cudaMemcpyHostToDevice) );
						cutilSafeCall( cudaMemcpy( d_particles.num_particles, 
							h_particles.num_particles, 
							sizeof(int), cudaMemcpyHostToDevice) );

						// Move lasers to the device
						cutilSafeCall( cudaMemcpy( d_lasers.pos, h_lasers.pos, 
							sizeof(float3) * num_lasers, cudaMemcpyHostToDevice) );
						cutilSafeCall( cudaMemcpy( d_lasers.vel, h_lasers.vel, 
							sizeof(float3) * num_lasers, cudaMemcpyHostToDevice) );
						cutilSafeCall( cudaMemcpy( d_lasers.num_lasers, 
							h_lasers.num_lasers, 
							sizeof(int), cudaMemcpyHostToDevice) );			

						// Evolve the system
						cutilCheckError( cutResetTimer(timer));
						cutilCheckError( cutStartTimer(timer));	

						// Evolve the system
						double print_counter = print_time;
						int random_index = 0;

						// evolve_laser<<<>>>() when we care about movement

						// Invoke particle evolution kernel
						evolve_particle(dimBlock, dimGrid, states, num_particles, d_particles, d_lasers, total_iterations);

						// We need the final particle positions and we need to
						// know if they are active
						cutilSafeCall( cudaMemcpy( h_particles.pos, d_particles.pos, 
							sizeof(float3) * num_particles, cudaMemcpyDeviceToHost) );
						cutilSafeCall( cudaMemcpy( h_particles.active, d_particles.active, 
							sizeof(bool) * num_particles, cudaMemcpyDeviceToHost) );

						// We need the evolved laser positions as well
						cutilSafeCall( cudaMemcpy( h_lasers.pos, d_lasers.pos, 
							sizeof(float3) * num_lasers, cudaMemcpyDeviceToHost) );
						/*
						//---------------- How long did this take? -----------------//
						cerr << "# Iteration time: " << cutGetTimerValue(timer) << " ms" << endl;

						//---------------- For Ashis, Avg distance travelled ---------------//

						cerr << "# resting_x resting_y resting_z dist_x dist_y dist_z dist_travelled" << endl;
						float3 starting_point = make_float3(0,0,-4e-5);
						float4 tot_dist_travelled = make_float4(0,0,0,0);
						size_t num_nans = 0;
						for(int p_idx=0; p_idx<num_particles; p_idx++)
						{
						cerr << h_particles.pos[p_idx].x << " " << h_particles.pos[p_idx].y << " " << h_particles.pos[p_idx].z << " ";

						double dist_x = abs( h_particles.pos[p_idx].x - starting_point.x );
						double dist_y = abs( h_particles.pos[p_idx].y - starting_point.y );
						double dist_z = abs( h_particles.pos[p_idx].z - starting_point.z );
						cerr << dist_x << " " << dist_y << " " << dist_z << " ";

						double dist_travelled = sqrt( 
						pow( h_particles.pos[p_idx].x - starting_point.x, 2) + 
						pow( h_particles.pos[p_idx].y - starting_point.y, 2) + 
						pow( h_particles.pos[p_idx].z - starting_point.z, 2) );
						cerr << dist_travelled << endl; 

						if (!isnan(dist_x)) {
						tot_dist_travelled.x += dist_x;
						} else { ++num_nans; }
						if (!isnan(dist_y)) {
						tot_dist_travelled.y += dist_y;
						} else { ++num_nans; }
						if (!isnan(dist_z)) {
						tot_dist_travelled.z += dist_z;
						} else { ++num_nans; }
						if (!isnan(dist_travelled)) {
						tot_dist_travelled.w += dist_travelled;
						} else { ++num_nans; }
						}
						cerr << "# NaNs " << num_nans << "\n";
						cerr << "# avg_x avg_y avg_z" << endl;
						cerr << (tot_dist_travelled.x / (float) num_particles) << " "
						<< (tot_dist_travelled.y / (float) num_particles) << " "
						<< (tot_dist_travelled.z / (float) num_particles) << " "
						<< (tot_dist_travelled.w / (float) num_particles) << endl;


						*/

						//  How long did the whole run take?
						cutStopTimer(timer);
						//cerr << "# Total evolution time: " << cutGetTimerValue(timer) - old_timer << " ms" << endl;
						cout << "# Elapsed time: " << cutGetTimerValue(timer) - old_timer << " ms" << endl;

						//---------------- Check to see which particles were trapped ----------------//	
						int trapped_count = 0;
						//              float3 a_pos = make_float3(0.0, 0.0, 0.0);
						//              float scale = 1.0f / num_particles;
						for(int p_idx=0; p_idx<num_particles; p_idx++)
						{
							float3 p_pos = h_particles.pos[p_idx];
							float3 l_pos = h_lasers.pos[p_idx];
							//  a_pos += p_pos * scale;
							if (h_particles.active[p_idx]) {                    
								//cerr << p_pos.x << " " << p_pos.y << " " << p_pos.z  << "\n";
								if( isTrapped(p_pos, l_pos, radial_trap_limit, axial_trap_limit) )
								{
									trapped_count++;
								}
							}
						}

						double trapped_percentage = (double) trapped_count / (double) num_particles;

						// y_dist z_dist %trapped
						#ifdef	_WINDOWS_64_
							out << y_dist << " " << z_dist << " " << trapped_percentage << endl;
						#else
							cout << y_dist << " " << z_dist << " " << trapped_percentage << endl;
						#endif
						//              cout << "a_pos = " << a_pos.x << ", " << a_pos.y << ", " << a_pos.z << "\n";
						//cout << flush;

					} // end of z_idx loop
				}  // end of y_idx loop
				


				//---------------- Cleanup and Quitting Starts Here ----------------//
				cerr << "# Evolution complete, performing cleanup." << endl;

				//---------------- Kill the CUDA stuff we no longer need -----------//
				// Cleanup particles
				delete[] h_particles.pos;
				delete[] h_particles.vel;
				delete[] h_particles.acc;
				delete[] h_particles.active;
				delete h_particles.num_particles;

				delete[] h_lasers.pos;
				delete[] h_lasers.vel;
				delete h_lasers.num_lasers;

				//delete h_particles;
				cutilSafeCall(cudaFree(d_particles.pos));
				cutilSafeCall(cudaFree(d_particles.vel));
				cutilSafeCall(cudaFree(d_particles.acc));
				cutilSafeCall(cudaFree(d_particles.active));
				//cutilSafeCall(cudaFree(d_particles));

				cutilSafeCall(cudaFree(d_lasers.pos));
				cutilSafeCall(cudaFree(d_lasers.vel));
				cutilSafeCall(cudaFree(d_lasers.num_lasers));

				// Cleanup force data
				delete[] h_force_array;
				cutilSafeCall(cudaFreeArray(d_force_array));

				// Cleanup random numbers
				delete[] h_random_numbers;
				//cutilSafeCall(cudaFree(d_random_numbers));
			}

		#ifdef	_WINDOWS_64_
		}
		#endif

		io::delete_loaded_force_data(force_data);


		cutStopTimer(timer);
		cutDeleteTimer(timer);		

		#ifdef	_WINDOWS_64_
			out.close();
		#endif
		// Kill CUDA and return
		cudaThreadExit();

	}catch( exception& e ){
		cerr << "Caught Exception: [" << e.what() << "]" << endl;
		abort();
	}

	#ifdef	_WINDOWS_64_
	//  cutilExit(argc, argv);
	cout<<"Press any key to continue..."<<endl;
	getch();
	#endif

	return 0;
}

