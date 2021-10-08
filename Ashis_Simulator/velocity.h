#ifndef VELOCITY_H
#define VELOCITY_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <blitz/array.h>
#include <blitz/tinyvec.h>

#include "cpuSphericalParticleSim.hpp"

using blitz::Array;
using blitz::TinyVector;

class Velocity {
public:
	Velocity ();

	double *calcVelInsideTrap (double *a, double *b, double *sv, double *sa, double radius, double interval, double t,      
                               const Array< blitz::TinyVector<double, 2>, 2 >& forcemap);


	double *calcPosnInsideTrap ();

	double *calcAccnInsideTrap ();
	double *calcForceInsideTrap ();

	double calcTimeInsideTrap ();

	~Velocity ();

private:
	int integrationCounter;					// must be odd
	double timeSim;

	double ** temp;

	double riParticle;
	double riOil;
	double riWater;

	double densityParticle;
	double densityWater;

	double mass;
	double massEff;

	double viscosity;
	double drag;
	double g;
	double absT;
		
	double c;
	double lambda;
	double power;
	double NA;

	double phiMax;
	double a;
	double omega0;

	double theta;					// Angle of incidence
	double r;

	double ra;						// radial variable
	double rmax;
	double f;
	
	double beta;						// azimuthal angle variable
	double alpha;

	double gammaD;
	double mu;

	double betaD;
	double d;

	double sDash;
	double sDDash;

	double pi;

	double Rs;
	double Rp;

	double Ts;
	double Tp;

	double R;
	double T;

	double qs;
	double qg;

	double *dF;
	double *F;

	double *vel;
	double *posn;
	double *accn;
	    
    void allocateArrays ();

	void calculateAngles (double radius);
	void calculateSingleRay ();

	void integrate (int i, int j);

	void calculateForceAshkin (double radius);

    void lookupForce(double *a, double *b,  
                     const Array< blitz::TinyVector<double, 2>, 2 >& forcemap );

	void calculateForce (double *p, double *q, double radius);
};

#endif
