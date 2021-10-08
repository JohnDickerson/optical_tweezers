/*********************************************************************************
 * Author - Arvind Balijepalli
 *
 * Change Log
 *
 *
 *********************************************************************************/
#ifndef __CPU_SPHERICAL_PARTICLE_SIM_
#define __CPU_SPHERICAL_PARTICLE_SIM_

#include "CPUCompute.hpp"
#include "SphericalParticle.hpp"
#include "RNG.hpp"

class cpuSphericalParticleSim {
public:
	cpuSphericalParticleSim();
	cpuSphericalParticleSim(double gamma, float T, double m, double dt);
	~cpuSphericalParticleSim();
	
	template <class T>
	sParticleparams3D<T> CPUPass(sParticleparams3D<T> startParams, sPoint3D<T> Fext);
	
protected:
	string 	mDataPath;
	bool	mTimingOnly;
	
	int *mAbortFlag;
	
	// pre-computed random number list
	ifstream msRandStream;
	vector< sPoint3D<double> > mvRandDeviates;
	int 	mCurrRandDevIdx;
	
	// simulation params
	double	mGamma;
	double	mEtta;
	double	mM;
	double	mDt;
};

template <class T>
sParticleparams3D<T> cpuSphericalParticleSim::CPUPass(sParticleparams3D<T> startParams, sPoint3D<T> Fext) {

	CPUCompute<T> s;	

	sParticleparams3D<T> t;
	
	sPoint3D<T> rand, gamma;
	RNG rnd(0,1,NORMAL);

	gamma.x=gamma.y=gamma.z		= (T) mGamma;
		
	
	rand.x=(T) rnd.rand();
	rand.y=(T) rnd.rand();
	rand.z=(T) rnd.rand();

    //    cerr << "\n\n mM = " << mM << "\n\n";
    //    cerr << "\n\n Fext = " << Fext.x << ", " << Fext.y << ", " << Fext.z  << "\n\n";
	t=s.driverVerlet(
			startParams,					// values at the start of timestep dt
			(T) mDt,						// dt
			gamma,							// gamma
			(T) mM,							// m
			(T) mEtta,						// Etta=sqrt(2*gamma*Kb*T/dt)
			Fext,							// Fext
			rand							// random number
	);
	
	return t;
}
#endif
