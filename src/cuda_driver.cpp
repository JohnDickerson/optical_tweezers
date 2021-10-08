

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <map>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/lagged_fibonacci.hpp>

#include "defines.h"

#include <GL/glew.h>
#include <vl/VisualizationLibrary.hpp>
#include <vlGLUT/GLUTWindow.hpp>

//#include "evolution_kernel.cu"
#include "constants.h"
#include "particle.h"
#include "io.h"
#include "laser.h"
#include "particle_list.h"
#include "position_log.hpp"
#include "tweezer_render_applet.hpp"

#include <cutil_inline.h>
#include "cutil.h"
#include "cutil_math.h"
#include "cuda_runtime.h"


//////  //  //  //  //  //  //  //////
//  REMEMBER TO STANDARDIZE UNITS!  //
//////  //  //  //  //  //  //  //////

using namespace tweezers;
using namespace std;

extern "C" void evolve_particle(const int& num_particles, particle_list particles, laser* lasers, int step);
extern "C" void setup_cuda_constants();
extern "C" void to_force_texture(float2* h_force_array, cudaArray* d_force_array,   cudaChannelFormatDesc& channelDesc, const int& force_tex_width, const int& force_tex_height);
extern "C" void copy_random_numbers(float3* h_random_numbers, size_t nbytes);

// RNG
boost::normal_distribution<float> norm_dist(0.0, 1.0);
boost::lagged_fibonacci19937 engine;

// Only works for Ashis' three examples
particle_list setup_particles(int num_particles)
{
  particle_list particles;
  particles.pos = new float3[num_particles];
  particles.vel = new float3[num_particles];
  particles.acc = new float3[num_particles];

  particles.pos[0] = make_float3(0, 0.000003, -0.000005);  // Fig 3.4
  particles.pos[1] = make_float3(0, 0.000003, 0.00000275); // Fig 3.5
  particles.pos[2] = make_float3(0, 0.000009, -0.000007);  // Fig 3.6
	
  particles.vel[0] = make_float3(0,0,0);
  particles.vel[1] = make_float3(0,0,0);
  particles.vel[2] = make_float3(0,0,0);

  particles.acc[0] = make_float3(0,0,0);
  particles.acc[1] = make_float3(0,0,0);
  particles.acc[2] = make_float3(0,0,0);

  return particles;
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

  // maps are internally sorted lowest -> highest
  io::forcemap::iterator y_it = force_data->begin();
  for(; y_it != force_data->end(); y_it++)
    {
      map<int, io::forcepair>* y_row = y_it->second;
      map<int, io::forcepair>::iterator z_it = y_row->begin();
      for(; z_it != y_row->end(); z_it++)
	{
	  io::forcepair fp = z_it->second;
	  force_array[(w_idx*width)+h_idx].x = fp.fy;
	  force_array[(w_idx*width)+h_idx].y = fp.fz;
	  //cout << z_it->first << endl;
	  //cout << fp.fy*1e12 << " " << fp.fz*1e12 << endl;
	  w_idx++;
			
	}
      h_idx++;
      w_idx = 0;
    }

  return force_array;
}

int main(int argc, char** argv)
{
  double total_time = 0.050;  //23ms for the Figure 3.4 runs
  double print_time = 0.00025;  //print every .25ms
  int total_iterations = (total_time / constants::time_step) + 1;

  //---------------- CPU-side Setup Starts Here ----------------//
  unsigned seed_value = 2345;
  srand(seed_value);
  cout << "# Seed: " << seed_value << endl;

  // Load pre-computed force data, [Y Z ForceY ForceZ]
  string forces_path = "force_5_microns_silica.txt";
  io::forcemap* force_data = io::load_force_data(forces_path, constants::dist_data_scale, constants::force_data_scale);
  //io::print_loaded_force_data(force_data);

	
  // Initiate some random particles
  int num_particles = 3;	
  particle_list h_particles = setup_particles(num_particles);
  position_log plog(num_particles);
  

  // Initialize laser with beam focus at origin
  int num_lasers = 1;
  //laser** h_lasers = setup_lasers(num_lasers);
  laser* h_lasers = new laser(0,0,0);

  // CPU-side create random number array for Brownian motion
  //int random_amount = (int)(total_time / constants::time_step) + 1;
  int random_amount = RANDOM_COUNT;
  float3* h_random_numbers = new float3[random_amount];
  generateRandomNumbers(h_random_numbers, random_amount);
	
  //---------------- Interactions with the GPU start here ----------------//
  cudaSetDevice( cutGetMaxGflopsDeviceId() );

  int deviceCount;
  cudaGetDeviceCount(&deviceCount);

  // This function call returns 0 if there are no CUDA capable devices.
  if (deviceCount == 0) {
    std::cerr << "There is no device supporting CUDA\n";
    std::abort();
  } else {
    cudaDeviceProp deviceProp;
    int dev;
    cudaGetDeviceProperties(&deviceProp, dev);

    if (dev == 0) {
      // This function call returns 9999 for both major & minor fields, if no CUDA capable devices are present
      if (deviceProp.major == 9999 && deviceProp.minor == 9999)
	std::cerr << "There is no device supporting CUDA.\n";
      else if (deviceCount == 1)
	std::cerr << "There is 1 device supporting CUDA\n";
      else
	std::cerr << "There are " << deviceCount << "devices supporting CUDA\n";
    }

  }
 
  // Transfer our initial particles
  particle_list d_particles;
  cutilSafeCall( cudaMalloc( (void**) &d_particles.pos, sizeof(float3) * num_particles ));
  cutilSafeCall( cudaMalloc( (void**) &d_particles.vel, sizeof(float3) * num_particles ));
  cutilSafeCall( cudaMalloc( (void**) &d_particles.acc, sizeof(float3) * num_particles ));
  cutilSafeCall( cudaMemcpy( d_particles.pos, h_particles.pos, sizeof(float3) * num_particles, cudaMemcpyHostToDevice) );
  cutilSafeCall( cudaMemcpy( d_particles.vel, h_particles.vel, sizeof(float3) * num_particles, cudaMemcpyHostToDevice) );
  cutilSafeCall( cudaMemcpy( d_particles.acc, h_particles.acc, sizeof(float3) * num_particles, cudaMemcpyHostToDevice) );

  // Transfer our initial laser location over as a 1D array
  laser* d_lasers = NULL;
  cutilSafeCall( cudaMalloc( (void**) &d_lasers, sizeof(laser) * num_lasers) );
  cutilSafeCall( cudaMemcpy( d_lasers, h_lasers, sizeof(laser) * num_lasers, cudaMemcpyHostToDevice) );


  // Transfer our force data over as a texture
  cudaArray* d_force_array;
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float2>();
  int force_tex_width = force_data->size();
  int force_tex_height = force_data->begin()->second->size();

  float2* h_force_array = form_force_array(force_data, force_tex_width, force_tex_height);

  to_force_texture(h_force_array, d_force_array, channelDesc, force_tex_width, force_tex_height);
	
  copy_random_numbers(h_random_numbers, sizeof(float3)*random_amount);


  // Transfer all the random physical constants and grid spacing stuff to constant memory
  setup_cuda_constants();


	
  //---------------- Evolution on the GPU Starts Here ----------------//
  cout << "# Setup complete, beginning evolution." << endl;

  // Create timer for later
  unsigned int timer = 0;
  cutilCheckError( cutCreateTimer( &timer));


  // Evolve the system
  cutilCheckError( cutStartTimer( timer));	
  double print_counter = print_time;
  int random_index = 0;
  for(int iter_count=0; iter_count < total_iterations; iter_count++)
    {
      // evolve_laser<<<>>>() when we care about movement

      // Invoke particle evolution kernel
      evolve_particle(num_particles, d_particles, d_lasers, random_index);

      //cudaThreadSynchronize();
		
      random_index+=num_particles;
      // If we've used all our random numbers in constant GPU memory, make some new ones
      if(random_index > RANDOM_COUNT-num_particles)
	{
	  generateRandomNumbers(h_random_numbers, random_amount);
	  copy_random_numbers(h_random_numbers, sizeof(float3)*random_amount);
	  random_index = 0;
	}

      // Every print_time seconds, update our output
      print_counter+=constants::time_step;

      if( print_counter >= print_time )
        {
          // only copy back the updated positions for now
          cutilSafeCall( cudaMemcpy( h_particles.pos, d_particles.pos, sizeof(float3) * num_particles, cudaMemcpyDeviceToHost) );
          plog.append(0, h_particles.pos[0]);
          plog.append(1, h_particles.pos[1]);
          plog.append(2, h_particles.pos[2]);

          cout << h_particles.pos[0].x * 1e6 << " " << h_particles.pos[0].y * 1e6 << " " << h_particles.pos[0].z * 1e6 << " "; 
          cout << h_particles.pos[1].x * 1e6 << " " << h_particles.pos[1].y * 1e6 << " " << h_particles.pos[1].z * 1e6  << " ";
          cout << h_particles.pos[2].x * 1e6 << " " << h_particles.pos[2].y * 1e6 << " " << h_particles.pos[2].z * 1e6 << endl;
          print_counter = 0.0;
        }
    }

  cout << "# Particle data as stored in the plog\n";
  cout << plog;

  //---------------- Cleanup and Quitting Starts Here ----------------//
  cout << "# Evolution complete, performing cleanup." << endl;

  //---------------- Kill the CUDA stuff we no longer need -----------//
  // Cleanup particles
  delete[] h_particles.pos;
  delete[] h_particles.vel;
  delete[] h_particles.acc;
  //delete h_particles;
  cutilSafeCall(cudaFree(d_particles.pos));
  cutilSafeCall(cudaFree(d_particles.vel));
  cutilSafeCall(cudaFree(d_particles.acc));
  //cutilSafeCall(cudaFree(d_particles));

  // Cleanup laser
  delete h_lasers;
  cutilSafeCall(cudaFree(d_lasers));

  // Cleanup force data
  io::delete_loaded_force_data(force_data);
  delete[] h_force_array;
  cutilSafeCall(cudaFreeArray(d_force_array));
	
  // Cleanup random numbers
  delete[] h_random_numbers;
  //cutilSafeCall(cudaFree(d_random_numbers));

  // Kill CUDA and return
  cudaThreadExit();
  

  //---------------------- Visualization Code -----------------------//
  // Visualization code goes here until it's properly tied in with
  // the simulation
  vl::fvec2 grid_range_x(0.0, 20.0); 
  vl::fvec2 grid_range_y(0.0, 20.0); 
  vl::fvec2 grid_range_z(-15, 8.0); 

  // Initalize GLUT and the Visualization Library
  glutInit(&argc, argv);
  vl::VisualizationLibrary::init();
  // Install Visualization Library shutdown function 
  atexit( vlGLUT::atexit_visualization_library_shutdown );
  
  // Setup the OpenGL context format
  vl::OpenGLContextFormat format;
  format.setDoubleBuffer(true);
  format.setRGBABits(8,8,8,8);  
  format.setDepthBufferBits(24);
  format.setStencilBufferBits(8);
  format.setFullscreen(false);
  
  // Create our applet
  vl::ref< tweezer_render_applet > tweezer_applet = 
    new tweezer_render_applet(grid_range_x, grid_range_y, grid_range_z, plog);

  tweezer_applet->initialize();
  
  // Create the hosting GLUT window
  vl::ref<vlGLUT::GLUTWindow> glut_window = new vlGLUT::GLUTWindow;
  // Bind the applet to the GLUT window to receive GUI events
  glut_window->addEventListener(tweezer_applet.get());
  
  // Set the render target
  vl::VisualizationLibrary::rendering()->as<vl::Rendering>()->setRenderTarget(glut_window->renderTarget());
  // Set the background 
  vl::fvec4 bg_color(0.8f, 0.8f, 0.8f, 1.0f);
  vl::VisualizationLibrary::rendering()->as<vl::Rendering>()->camera()->viewport()->setClearColor(bg_color);

  // Define the camera position and orientation
  float3 pos = plog.get_position(0,0);
  vl::vec3 eye((1e6 * pos.x) + 10.0f, -1e6 * pos.z, (1e6 * pos.y) + 50.0f);
  vl::vec3 center(1e6 * pos.x, -1e6 * pos.z, 1e6 * pos.y);
  vl::vec3 up(0.0f, 1.0f, 0.0f);
  vl::mat4 view_mat = vl::mat4::lookAt(eye, center, up).inverse();
  vl::VisualizationLibrary::rendering()->as<vl::Rendering>()->camera()->setViewMatrix( view_mat );

  int x = 0;
  int y = 0;
  int width = 1024;
  int height = 1024;
  glut_window->initGLUTWindow("Nano-Tweezers Visualization", format, x, y, width, height);

  glutMainLoop();
  //----------- End of visualization code --------------------------//

  //  cutilExit(argc, argv);
  return 0;
}

