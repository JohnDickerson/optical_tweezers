/*********************************************************************************
 * Author - Arvind Balijepalli
 *
 * Change Log
 *
 *
 *********************************************************************************/

#include "cpuSphericalParticleSim.hpp"
#include <cstdlib>
	
cpuSphericalParticleSim::cpuSphericalParticleSim(double gamma, float T, double m, double dt) {
	// set sim params
	mGamma					= gamma;
	mEtta					= sqrt ((2*gamma*Kb*T)/dt);

	mM						= m;			
	mDt						= dt;
}

cpuSphericalParticleSim::cpuSphericalParticleSim() {
	cerr << "Cannot use default constructor for class cpuSphericalParticleSim\n";
	exit(-1);
}

cpuSphericalParticleSim::~cpuSphericalParticleSim() {
}
