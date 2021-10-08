// C++ headers
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <map>
#include <ctime>

#define _WINDOWS_64_
#ifdef _WINDOWS_64_
	#include <fstream>
	#include <iomanip>
	#include <conio.h>
#endif

// Boost headers
#include <boost/random/normal_distribution.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/iterator/counting_iterator.hpp>

// C/C++ for CUDA headers
#include <vector_types.h>
#include <cutil_math.h>

// Our headers
#include "constants.h"
#include "particle_list.h"
#include "io.h"
#include "laser.h"
#include "interpolation.hpp"
#include "arg_parser.hpp"

/**
 *  REMEMBER TO STANDARDIZE UNITS!
 **/
using namespace tweezers;
using namespace std;
using namespace blitz;

/** Namespace declarations */
namespace po = boost::program_options;

boost::normal_distribution<double> norm_dist(0.0, 1.0);
boost::lagged_fibonacci19937 engine;

typedef Array< TinyVector<float, 2>, 2 > force_array_t;

//runs simulation several times to compute the trapping probability
void checkTrappingProbability(unsigned int simulationCount);

// Evolve the laser's position / velocity / acceleration / etc
void do_laser_time_step(float3& laserbeam)
{
  // Right now, laser is stationary
  // Fix calculateF before moving laser
}


// Calculates the effect of gravity, buoyancy, and trapping force on particle
// Gravity = 0 in X, 0 in Y, +9.81m/s^2 in Z+
// +Z is down, -Z is up
// Y is aligned with the direction of motion of the laser
inline void calculateF(float3& part_pos, float3& l_pos, 
                       force_array_t& force_data, float3& output_force) {
  
  // XY forces are symmetric, distance from particle -> laser
  float relYPos = sqrtf( (part_pos.x - l_pos.x)*(part_pos.x - l_pos.x) + 
                         (part_pos.y - l_pos.y)*(part_pos.y - l_pos.y));  
  // Z forces are not symmetric
  float relZPos = part_pos.z - l_pos.z;  
	
  // We are told to ignore laser forces outside of the grid,
  // so just update the particle with gravity and buoyancy
  if(relYPos < constants::grid_Y_extent_low || 
     relYPos + constants::grid_spacing > constants::grid_Y_extent_high ||
	 relZPos < constants::grid_Z_extent_low || 
     relZPos + constants::grid_spacing > constants::grid_Z_extent_high) {
    float3 gb = make_float3(0.0f, 0.0f, (constants::langevin_F_g.z - constants::langevin_F_b.z)*9.81);
    output_force = gb;
    //cout << "Out of bounds" << endl;
    return;

  }
	
  TinyVector<float, 2> interp_res;
  interpolation::bi_interp_lookup(force_data, relYPos, relZPos, interp_res);

  // Print out the force 
  //  cerr << "F[" << relYPos * 1e6 << " microns, " << relZPos * 1e6 << " microns ] = " 
  //       << "< " << interp_res(0) << ", " << interp_res(1) << " >\n";
	
  // Need to move XY forces from local (laser) coordinate system to world coordinate system
  // The calculated Y force is in the direction of the laser
  float xDif = part_pos.x - l_pos.x;
  float yDif = part_pos.y - l_pos.y;


  // float mag = sqrt(xDif*xDif + yDif*yDif);

  // if(mag > 1e-20) {
  //   xDif /= mag;
  //   yDif /= mag;
  // } else {
  //   xDif = 0.0;
  //   yDif = 0.0;
  // }
    
  // float xc = xDif*interp_res(0);
  // float yc = yDif*interp_res(0);
  // float zc = interp_res(1);
  // float3 fInterpVec = make_float3(xc, yc, zc);
	

  float angle = atan2f(yDif, xDif);
  float xForce = interp_res(0) * cosf(angle);
  float yForce = interp_res(0) * sinf(angle);
  float zForce = interp_res(1);

  if( fabsf(sqrtf(xDif*xDif + yDif*yDif)) < 1e-18 )
    {	xForce = 0; yForce = 0; }
  if( fabsf(interp_res(1)) < 1e-20 )
    { zForce = 0; }
  float3 fInterpVec = make_float3(xForce, yForce, zForce);



  // v represets the F_ext force vector in the paper -- this is the sum
  // of gravity (in +Z direction), buoyancy forces (in -Z direction), and
  // the trapping force estimated above
  float3 gb = make_float3(0.0f, 0.0f, (constants::langevin_F_g.z - constants::langevin_F_b.z)*9.81);
  output_force =  gb + fInterpVec;
	
  return;
}

void do_particle_time_step(particle_list& particles, size_t& num_particles, 
                           float3& laserbeam, force_array_t& force_data, bool do_print) {

  // -gamma/m    * V(t)
  double c1 = -constants::langevin_gamma * constants::inv_sphere_mass;
  // 1/m * sqrt(xi / timestep)    * N[0,1]
  double c2 = constants::inv_sphere_mass * sqrt( constants::langevin_xi_squared / constants::time_step );
	
  // 1/m * (f_g + f_b + f_t)
	
  for(size_t i = 0; i < num_particles; ++i) {
    float3 vel = particles.vel[i];
    float3 pos = particles.pos[i];
    float3 acc = particles.acc[i];

    float3 first_term = vel * c1;
		
    double rand01 = norm_dist.operator () <boost::lagged_fibonacci19937>((engine));
    double rand02 = norm_dist.operator () <boost::lagged_fibonacci19937>((engine));
    double rand03 = norm_dist.operator () <boost::lagged_fibonacci19937>((engine));

    float3 second_term = make_float3(c2 * rand01, c2 * rand02, c2 * rand03);
		
    float3 third_term = make_float3(0.0, 0.0, 0.0);
    calculateF(pos, laserbeam, force_data, third_term);
    third_term *= constants::inv_sphere_mass;
		
    if (do_print) {
  
      // // XY forces are symmetric, distance from particle -> laser
      // float relYPos = sqrtf( (pos.x - 0.0)*(pos.x - 0.0) + 
      //                        (pos.y - 0.0)*(pos.y - 0.0));  
      // // Z forces are not symmetric
      // float relZPos = pos.z - 0.0;  
		
      // TinyVector<float, 2> interp_res;
      // interpolation::bi_interp_lookup(force_data, relYPos, relZPos, interp_res);

      // cout << "F[" << relYPos << ", " << relZPos << "] = " << interp_res(0) << ", " << interp_res(1) << "\n";
      // cout << "First: " << first_term.x << " " << first_term.y << " " << first_term.z << endl;
      // cout << "Second: " << second_term.x << " " << second_term.y << " " << second_term.z << endl;
      // cout << "Third: " << third_term.x << " " << third_term.y << " " << third_term.z << endl << endl;
    }
    // a(t+dt) = v(t+dt) - v(t) / dt = c1+c2+c3
    float3 accel = first_term + second_term + third_term;
		
    // Evolve particle's physics using second order Verlet velocity integrator
    // Explanation: http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
    particles.pos[i] = 
      pos + 
      vel * constants::time_step + 
      acc * constants::time_step * constants::time_step*0.5;

    particles.vel[i] = 
      vel + 
      (acc + accel) * 0.5 * constants::time_step;

    particles.acc[i] = accel;
  }

}

// Initiates particles within the grid boundaries
// Just initiates them in a line parallel to Y axis for now,
// with no initial velocity
vector<particle*>* setup_particles(int num_particles, float x_position, float y_position, float z_position)
{
  vector<particle*> *particles = new vector<particle*>();
	
  double y_start = -1;
  for(int i=0; i<num_particles; i++)
	{
      particle* p = new particle(x_position, y_position, z_position);
      particles->push_back(p);
	}
	
  return particles;
}

bool isTrapped(const float3& particle, const float3& laser)
{
  //cout << particle.x << " " << particle.y << " " << particle.z << endl;
  if ( fabs(particle.x - laser.x) < 5e-7 && fabs(particle.y - laser.y) < 5e-7 
       && (particle.z - laser.z) >= 0.0 && (particle.z - laser.z) < 10e-7 ) {
    return true;
  } else {
    return false;
  }
}

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

void teardown_particles(particle_list& l) {
  delete [] l.pos;
  delete [] l.vel;
  delete [] l.acc;
  delete [] l.active;
  delete l.num_particles;
}

void forcemap_to_array(io::forcemap* forcemap, force_array_t& forcearray) {

  forcearray.resize(81, 113);
  double ymin = 0.0;
  double ymax = 20.0;
  double zmin = -20.0;
  double zmax = 8.0;

  io::forcemap::iterator it;
  for (it = forcemap->begin(); it != forcemap->end(); ++it) {
    float cy = (it->first.x * 1e6);
    float cz = (it->first.y * 1e6);
    int y = ceil( ((cy-ymin)/(ymax-ymin)) * (81.0-1));
    int z = ceil( ((cz-zmin)/(zmax-zmin)) * (113.0-1));
    forcearray(y,z)[0] = static_cast<double>(it->second.x);
    forcearray(y,z)[1] = static_cast<double>(it->second.y);
  }
  //  std::abort();
}

int main(int argc, char** argv) 
{
	//runs simulation several times to compute the trapping probability
	checkTrappingProbability(50);

	// Keep track of seed value!
	unsigned seed_value = 2345;
	srand(seed_value);
	cout << "# Seed: " << seed_value << endl;
	try 
	{
		#ifndef	_WINDOWS_64_
			po::options_description desc("Allowed Options");
			desc.add_options()
				("help,h", "produce this message")
				("y_start,l", po::value<float>(), "pdb file")
				("y_end,r", po::value<float>(), "chain file")
				("z_start,t", po::value<float>(), "configuration file")
				("z_end,b", po::value<float>(), "curvature output file")
				("duration,d", po::value<float>(), "duration of the simulation")
				("num_particles,n", po::value<size_t>(), "# of particles")
				;

			ArgParser ap( desc );
			ap.parse_args( argc, argv );
			float total_time = ap["duration"].as<float>();
			float y_start = ap["y_start"].as<float>();
			float y_end = ap["y_end"].as<float>();
			float z_start = ap["z_start"].as<float>();
			float z_end = ap["z_end"].as<float>();
			size_t num_particles = ap["num_particles"].as<size_t>();//atoi(argv[1]);//1521;	


			//-l 0e-6 -r .99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100 >> c
			//-l 1.0e-6 -r 1.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100 >
			//-l 2.0e-6 -r 2.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100 >
			//-l 3.0e-6 -r 3.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100 >
			//-l 4.0e-6 -r 4.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100 >
			//-l 5.0e-6 -r 5.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100 >
			//-l 6.0e-6 -r 6.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100 >
			//-l 7.0e-6 -r 7.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100 >
			//-l 8.0e-6 -r 8.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100 >
			//-l 9.0e-6 -r 9.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100 >
			//-l 10.0e-6 -r 10.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100
			//-l 11.0e-6 -r 11.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100
			//-l 12.0e-6 -r 12.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100
			//-l 13.0e-6 -r 13.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100
			//-l 14.0e-6 -r 14.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100
			//-l 15.0e-6 -r 15.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100
			//-l 16.0e-6 -r 16.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100
			//-l 17.0e-6 -r 17.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100
			//-l 18.0e-6 -r 18.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100
			//-l 19.0e-6 -r 19.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100
			//-l 20.0e-6 -r 20.99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100

			// Load pre-computed force data, [Y Z ForceY ForceZ]
			#ifndef	WIN32
				string forces_path = "/chimerahomes/sujal/tweezers/force_5_microns_silica.txt";	
			#else
				string forces_path = "..\\..\\force_5_microns_silica.txt";	
			#endif

		#else
			// Load pre-computed force data, [Y Z ForceY ForceZ]
			string forces_path = "..\\..\\force_5_microns_silica.txt";	
		#endif
		io::forcemap* force_data = io::load_force_data(forces_path,constants::dist_data_scale,constants::force_data_scale);

		#ifdef	_WINDOWS_64_
		ofstream out("..\\outputCPU.txt");
		out<<setprecision(9);
		for(unsigned y=0;y<21;y++)
		{
			float total_time = 1.0f;
			float y_start = y*1.0e-6f;
			float y_end = y_start+.99e-6f;
			float z_start = -20.0e-6f;
			float z_end = 8.1e-6f;
			size_t num_particles = 100;	

			out<<"\n#"<<y<<"\t";
			out<<"-l "<<y_start<<" ";
			out<<"-r "<<y_end<<" ";
			out<<"-t "<<z_start<<" ";
			out<<"-b "<<z_end<<" ";
			out<<"-d "<<total_time<<" ";
			out<<"-n "<<num_particles<<"\n";

		#endif

			int total_iterations = (total_time / constants::time_step) + 1;
			//io::print_loaded_force_data(force_data);

			double print_time = total_time+1;  //print every .25ms

			float steps_per_micron = constants::grid_spacing;//2.5e-7f;
			float test_step = constants::grid_spacing;//2.5e-7f;
			size_t num_y_steps = ceil(static_cast<float>(std::abs(y_end - y_start)) / steps_per_micron);
			size_t num_z_steps = ceil(static_cast<float>(std::abs(z_end - z_start)) / steps_per_micron);

			cout << "# Max time: " << total_time << " seconds" << endl;
			cout << "# Timestep: " << constants::time_step * 1e6 << " microsec" << endl;
			cout << "# Num particles: " << num_particles << " per particle per run" << endl;
			cout << "# y steps: " << num_y_steps << endl;
			cout << "# z steps: " << num_z_steps << endl;
			cout << "# Iterations: " << total_iterations << " per particle per run" << endl;
			cout << "# y: [" << y_start << ", " << y_end << "]" << endl;
			cout << "# z: [" << z_start << ", " << z_end << "]" << endl;
			cout << "# Test step: " << test_step * 1e6 << " microns." << endl;


			force_array_t force_array;
			forcemap_to_array(force_data, force_array);

			// tracks pos, vel, acc, and active in float3*s
			particle_list particles = setup_particles(num_particles);

			boost::counting_iterator<size_t> yb;
			boost::counting_iterator<size_t> ye = boost::make_counting_iterator<size_t>(num_y_steps);

			boost::counting_iterator<size_t> zb;
			boost::counting_iterator<size_t> ze = boost::make_counting_iterator<size_t>(num_z_steps);

			for(zb = boost::make_counting_iterator<size_t>(0); zb != ze; ++zb) 
			{
				for(yb = boost::make_counting_iterator<size_t>(0); yb != ye; ++yb) 
				{
					float zpos = z_start + ((*zb) * steps_per_micron);
					float ypos = y_start + ((*yb) * steps_per_micron);

					float x_position = 0e-6;			// x_idx constant at 0 microns
					float y_position = ypos;//y_idx * 1e-6;  // y_idx from 0 to 20 microns
					float z_position = zpos;//-3e-6;			// z_idx constant at -15 microns

					// Initiate some random particles
					for(int i=0; i<num_particles; i++) 
					{
						particles.pos[i] = make_float3(x_position,y_position,z_position);
						particles.vel[i] = make_float3(0.0f, 0.0f, 0.0f);
						particles.acc[i] = make_float3(0.0f, 0.0f, 0.0f);
						particles.active[i] = true;
					}

					*(particles.num_particles) = num_particles;

					// Initialize laser with beam focus at origin
					float3 laserbeam = make_float3(0.0f, 0.0f, 0.0f);

					// Evolve the system
					clock_t start,finish;
					double time;
					start = clock();

					double print_counter = print_time;
					bool doprint = false;
					for(double curtime=0; curtime < total_time; curtime+=constants::time_step) 
					{

						if(print_counter >= print_time) 
						{
							doprint = true;
						} else {
							doprint = false;
						}

						do_laser_time_step(laserbeam);
						do_particle_time_step(particles, num_particles, laserbeam, force_array, doprint);

						/*
						if(print_counter >= print_time) {
						//io::print_time_step(particles, laserbeam, curtime);
						for(size_t i = 0; i < num_particles; ++i) {
						cout << particles.pos[i].x*1e6 << " " << particles.pos[i].y*1e6 << " " << particles.pos[i].z*1e6 << endl;
						print_counter = 0;
						}
						}
						*/
						print_counter+=constants::time_step;
					}

					// Check to see how many particles were trapped
					int trapped_count = 0;
					for(int p_idx=0; p_idx<num_particles; p_idx++) 
					{
						float3 p_pos = particles.pos[p_idx];
						float3 l_pos = laserbeam;
						if (particles.active[p_idx]) 
						{
							if( isTrapped(p_pos, l_pos) )
							{
								++trapped_count;
							}
						}
					}

					// y_dist z_dist %trapped
					#ifdef	_WINDOWS_64_
						out<<y_position << " " << z_position << " "<< (double) trapped_count / (double) num_particles << endl;
					#else
						cout << y_position << " " << z_position << " "<< (double) trapped_count / (double) num_particles << endl;
					#endif
					//cout << flush;

					finish = clock();
					time = (double(finish)-double(start))/CLOCKS_PER_SEC;

					//cout << num_particles << " " << time << endl;		
					cout << "# Elapsed time: " << (time*1000) << " ms" << endl;

				} // end z_idx loop
			} // end y_idx loop
			// Cleanup particles
			teardown_particles(particles);
		#ifdef	_WINDOWS_64_
		}
		out.close();
		#endif
		//// Cleanup force data
		io::delete_loaded_force_data(force_data);
	}
	catch( exception& e ) 
	{
		cerr << "Caught Exception: [" << e.what() << "]" << endl;
		abort();
	}
	#ifdef	_WINDOWS_64_
	cout<<"Press any key to continue..."<<endl;
	getch();
	#endif
	return 0;
}

//runs simulation several times to compute the trapping probability
void checkTrappingProbability(unsigned int simulationCount)
{
	unsigned seed_value = 2345;
	srand(seed_value);

	float total_time = 1.0f;
	float y_start = 5.5e-06;//1.875e-05;
	float z_start = 3.75e-06;//-8.25e-06f;
	float y_end = y_start+2.5e-8f;
	float z_end = z_start+2.5e-8f;

	int total_iterations = (total_time / constants::time_step) + 1;
	//io::print_loaded_force_data(force_data);

	double print_time = total_time+1;  //print every .25ms

	float steps_per_micron = constants::grid_spacing;//2.5e-7f;
	float test_step = constants::grid_spacing;//2.5e-7f;
	size_t num_y_steps = ceil(static_cast<float>(std::abs(y_end - y_start)) / steps_per_micron);
	size_t num_z_steps = ceil(static_cast<float>(std::abs(z_end - z_start)) / steps_per_micron);

	//cout << "# Max time: " << total_time << " seconds" << endl;
	//cout << "# Timestep: " << constants::time_step * 1e6 << " microsec" << endl;
	//cout << "# Num particles: " << num_particles << " per particle per run" << endl;
	//cout << "# y steps: " << num_y_steps << endl;
	//cout << "# z steps: " << num_z_steps << endl;
	//cout << "# Iterations: " << total_iterations << " per particle per run" << endl;
	//cout << "# y: [" << y_start << ", " << y_end << "]" << endl;
	//cout << "# z: [" << z_start << ", " << z_end << "]" << endl;
	//cout << "# Test step: " << test_step * 1e6 << " microns." << endl;
	

	string forces_path = "..\\..\\force_5_microns_silica.txt";	
	io::forcemap* force_data = io::load_force_data(forces_path,constants::dist_data_scale,constants::force_data_scale);

	ofstream out("trapProb.txt");
	out<<setprecision(9);

	for(unsigned x=5;x<9;x++)
	{
		size_t num_particles=pow(2.0,(double)x);
		num_particles=num_particles;

		for(unsigned y=0;y<simulationCount;y++)
		{
			force_array_t force_array;
			forcemap_to_array(force_data, force_array);

			// tracks pos, vel, acc, and active in float3*s
			particle_list particles = setup_particles(num_particles);

			boost::counting_iterator<size_t> yb;
			boost::counting_iterator<size_t> ye = boost::make_counting_iterator<size_t>(num_y_steps);

			boost::counting_iterator<size_t> zb;
			boost::counting_iterator<size_t> ze = boost::make_counting_iterator<size_t>(num_z_steps);

			for(zb = boost::make_counting_iterator<size_t>(0); zb != ze; ++zb) 
			{
				for(yb = boost::make_counting_iterator<size_t>(0); yb != ye; ++yb) 
				{
					float zpos = z_start + ((*zb) * steps_per_micron);
					float ypos = y_start + ((*yb) * steps_per_micron);

					float x_position = 0e-6;			// x_idx constant at 0 microns
					float y_position = ypos;//y_idx * 1e-6;  // y_idx from 0 to 20 microns
					float z_position = zpos;//-3e-6;			// z_idx constant at -15 microns

					// Initiate some random particles
					for(int i=0; i<num_particles; i++) 
					{
						particles.pos[i] = make_float3(x_position,y_position,z_position);
						particles.vel[i] = make_float3(0.0f, 0.0f, 0.0f);
						particles.acc[i] = make_float3(0.0f, 0.0f, 0.0f);
						particles.active[i] = true;
					}

					*(particles.num_particles) = num_particles;

					// Initialize laser with beam focus at origin
					float3 laserbeam = make_float3(0.0f, 0.0f, 0.0f);

					// Evolve the system
					clock_t start,finish;
					double time;
					start = clock();

					double print_counter = print_time;
					bool doprint = false;
					for(double curtime=0; curtime < total_time; curtime+=constants::time_step) 
					{

						if(print_counter >= print_time) 
						{
							doprint = true;
						} else {
							doprint = false;
						}

						do_laser_time_step(laserbeam);
						do_particle_time_step(particles, num_particles, laserbeam, force_array, doprint);
						print_counter+=constants::time_step;
					}

					// Check to see how many particles were trapped
					int trapped_count = 0;
					for(int p_idx=0; p_idx<num_particles; p_idx++) 
					{
						float3 p_pos = particles.pos[p_idx];
						float3 l_pos = laserbeam;
						if (particles.active[p_idx]) 
						{
							if( isTrapped(p_pos, l_pos) )
							{
								++trapped_count;
							}
						}
					}
					cout <<"Num Particles: "<<num_particles<<"\t Trapped particles : "<<trapped_count<<"\tRatio : "<<(double) trapped_count / (double) num_particles << endl;
					out <<num_particles<<"\t"<<trapped_count<<"\t"<<((double) trapped_count/(double)num_particles)<<endl;
					finish = clock();
					time = (double(finish)-double(start))/CLOCKS_PER_SEC;
				} // end z_idx loop
			} // end y_idx loop
			// Cleanup particles
			teardown_particles(particles);
		}
	}
	out.close();
	// Cleanup force data
	io::delete_loaded_force_data(force_data);
	cout<<"Press any key to continue..."<<endl;
	getch();
	exit(0);
}
