#ifndef __TWEEZERS_CONSTANTS_H__
#define __TWEEZERS_CONSTANTS_H__

#include <cmath>

#include <vector_types.h>
#include <vector_functions.h>

//sujal
//added for visual studio
#ifdef _MSC_VER
	#ifndef M_PI
		#define M_PI 3.141592653589793238462643 
	#endif
	#ifndef M_SQRT2
		#define M_SQRT2 1.4142135623730950488016887 
	#endif
#endif

namespace tweezers
{
	// Parameters used in Ashkin's ray optics model
	// Right-hand coordinate system such that +z points vertically downwards
	// Force values (picoNewtons) on the +YZ plane using a rectangular grid with uniform spacing of 0.25 micron
	namespace constants
	{
		// Silica sphere details
		const double sphere_diameter = 5.0e-6;   // 5 microns
		const double sphere_radius = sphere_diameter * 0.5;
		const double sphere_n = 1.46;         //refractive index
		const double sphere_density = 2000.0; // 2000 kg/m^3 
		const double sphere_mass = sphere_density * (4.0/3.0) * M_PI * pow(sphere_radius, 3.0);
		const double inv_sphere_mass = 1.0 / sphere_mass;
		
		// Fluid (water) details
		const double fluid_viscosity = 1.002e-3; // SI units
		const double fluid_n = 1.33;             //refractive index
		const double fluid_density = 1000.0; // 1000 kg/m^3 -> kg/micron^3
		const double fluid_temperature = 293.15;    //K
		// Fluid (oil) details
		const double oil_n = 1.51f; //refractive index
		
		// Gaussian laser beam, circular polarization, propagates downward
		const double laser_power = 0.1;          //W
		const double laser_wavelength = 532;     //nm 
		const double laser_waist_radius = 0.4;   //microns
		const double laser_half_cone_angle = 67; //degrees
		
		// Lens details
		const double lens_aperture = 1.4;      //f-number, focal length / aperture diameter
		const double lens_focal_length = 3.33; //mm
		
		// Natural constants
		const float3 g = {0,0,9.81};      //gravity, m/s^2  gravity is in positive Z direction
		const double c = 3.0e8; //speed of light, m/s
		const double boltzmann = 1.380650424e-23; //SI = m^2*kg / s^2*K
		//const double boltzmann = 1.38065e-16;
		
		// Grid details
		const double grid_spacing = 0.25e-6;     //0.25 microns
		const double grid_inv_spacing = 1.0/grid_spacing;
		const double grid_Y_extent_low = 0.0e-6;   //0 microns
		const double grid_Y_extent_high = 20.0e-6; //20 microns
		const double grid_Z_extent_low = -20.0e-6; //-15 microns
		const double grid_Z_extent_high = 8.0e-6;  //8 microns
	
		// Evolution Constants
		const double langevin_gamma = 6.0 * M_PI * fluid_viscosity * sphere_radius;
		const double langevin_xi_squared = 2.0 * langevin_gamma * boltzmann * fluid_temperature;

		const float3 langevin_F_g = make_float3((4.0/3.0) * M_PI * sphere_density * pow(sphere_radius, 3.0),
                                                (4.0/3.0) * M_PI * sphere_density * pow(sphere_radius, 3.0),
                                                (4.0/3.0) * M_PI * sphere_density * pow(sphere_radius, 3.0));

		const float3 langevin_F_b = make_float3((4.0/3.0) * M_PI * fluid_density * pow(sphere_radius, 3.0),
                                                (4.0/3.0) * M_PI * fluid_density * pow(sphere_radius, 3.0),
                                                (4.0/3.0) * M_PI * fluid_density * pow(sphere_radius, 3.0));
		
		// Closest multiple of 100ns that is less than the relaxation time = m / gamma
		const double time_step = static_cast<int>((sphere_mass / langevin_gamma) * 1.0e7) * (double) 1.0e-7;

        const double c1 = -langevin_gamma * inv_sphere_mass;
        const double c2 = inv_sphere_mass * sqrt( langevin_xi_squared / time_step );
		//const double time_step = 5.0e-7;
		
		// Force data's scale
		const double dist_data_scale = 1.0e-6;  // distance are given in microns
		const double force_data_scale = 1.0e-12;  // forces given in picoNewtons
	}
}

#endif

