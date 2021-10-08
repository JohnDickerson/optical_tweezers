#ifndef __CPUCOMPUTE_HPP_
#define __CPUCOMPUTE_HPP_

#include "SphericalParticle.hpp"
#include "PhysicalConstants.hpp"
#include <cmath>
#include <vector>
#include <cstdio>
#include <iostream>

using namespace std;

template <class Type>
class CPUCompute {
protected:

public:
	FILE *f1;

	CPUCompute();
	~CPUCompute();
	

	sPoint3D<Type> acceleration(sPoint3D<Type> V, Type dt, sPoint3D<Type> gamma, Type m, Type etta, sPoint3D<Type> Fext, sPoint3D<Type> rand);
	Type acceleration(Type V, Type dt, Type gamma, Type m, Type etta, Type Fext, Type rand);

	sParticleparams3D<Type> driverVerlet(sParticleparams3D<Type> valT, Type dt, sPoint3D<Type> gamma, Type m, Type etta, sPoint3D<Type> Fext, sPoint3D<Type> rand);
	sParticleparams1D<Type> driverVerlet(sParticleparams1D<Type> valT, Type dt, Type gamma, Type m, Type etta, Type Fext, Type rand);

};

template <class Type>
CPUCompute<Type>::CPUCompute() {}

template <class Type>
CPUCompute<Type>::~CPUCompute() {}

template <class Type>
sPoint3D<Type> CPUCompute<Type>::acceleration(sPoint3D<Type> V, Type dt, sPoint3D<Type> gamma, Type m, Type etta, sPoint3D<Type> Fext, sPoint3D<Type> rand) {
	sPoint3D<Type> acc;

	acc.x=acceleration(V.x, dt, gamma.x, m, etta, Fext.x, rand.x);
	acc.y=acceleration(V.y, dt, gamma.y, m, etta, Fext.y, rand.y);
	acc.z=acceleration(V.z, dt, gamma.z, m, etta, Fext.z, rand.z);

	return acc;
}

template <class Type>
Type CPUCompute<Type>::acceleration(Type V, Type dt, Type gamma, Type m, Type etta, Type Fext, Type rand) {
  float F = (-gamma/m) * V +
			rand/m * etta +
			Fext/m;
  // cerr << "Fext = " << Fext << "\n";
  // cerr << "m = " << m << "\n";
  // cerr << "rand = " << rand << "\n";
  // cerr << "etta = " << etta << "\n";
  // cerr << "Langevin = " << (-gamma/m) * V 
  //      << ", Brownian = "  << rand/m * etta 
  //      << ", Fext / m = " << Fext / m;
  
  return	F;
	/*return	(-gamma/m) * V +
			Fext/m;*/
}

template <class Type>
sParticleparams3D<Type> CPUCompute<Type>::driverVerlet(sParticleparams3D<Type> valT, Type dt, sPoint3D<Type> gamma, Type m, Type etta, sPoint3D<Type> Fext, sPoint3D<Type> rand) {
	sParticleparams3D<Type> final;
	sParticleparams1D<Type> x, y, z, t;
	
	// A kludgy way to do the 3D simulation is to split it into three 1-D steps. Since there is no coupling
	// between the three dimensions in a free particle this shouldnt matter, but there may be timing issuse
	// Profile this piece of code to see what the timing implications are.
	x.t=valT.t;		y.t=valT.t;		z.t=valT.t;
	x.r=valT.r.x;	y.r=valT.r.y;	z.r=valT.r.z;
	x.v=valT.v.x;	y.v=valT.v.y;	z.v=valT.v.z;
	x.a=valT.a.x;	y.a=valT.a.y;	z.a=valT.a.z;

    //    cerr << "Posn: " << x.r << ", " << y.r << ", " << z.r << "\n";
    //    cerr << "Forces: \n";
    //	cerr << "x: ";
	t=driverVerlet(x, dt, gamma.x, m, etta, Fext.x, rand.x);
	final.r.x=t.r;	final.v.x=t.v;	final.a.x=t.a;
    //    cerr << "\n";

    //	cerr << "y: ";
	t=driverVerlet(y, dt, gamma.y, m, etta, Fext.y, rand.y);
	final.r.y=t.r;	final.v.y=t.v;	final.a.y=t.a;
    //    cerr << "\n";

    //	cerr << "z: ";
	t=driverVerlet(z, dt, gamma.z, m, etta, Fext.z, rand.z);
	final.r.z=t.r;	final.v.z=t.v;	final.a.z=t.a;
    //    cerr << "\n";
	
	final.t=t.t;
	
    //    cerr << "\n\n";

	return final;
}

template <class Type>
sParticleparams1D<Type> CPUCompute<Type>::driverVerlet(sParticleparams1D<Type> valT, Type dt, Type gamma, Type m, Type etta, Type Fext, Type rand) {

	/*t = valT[[1]]; r = valT[[2]]; v = valT[[3]]; a = valT[[4]];
	tt = t + dt;
	rr = r + dt v  + (1/2) dt^2 a;
	aa = acceleration[v, dt, gamma, m, T, Fext, rand]; 
	vv = ( v + ((1/2) dt a) ) + (1/2) dt aa; 
	
	Return[{tt, rr, vv,aa}]*/
	
	sParticleparams1D<Type> final;
	
	final.t = valT.t + dt;
	final.r = valT.r + (valT.v * dt)  + (0.5 * valT.a * dt * dt);

	final.a = acceleration(valT.v, dt, gamma, m, etta, Fext, rand);
	//final.r = valT.r + (valT.v * dt)  + (0.25 * (valT.a + final.a) * dt * dt); 
	final.v = valT.v + ((valT.a + final.a) * (0.5*dt)); 
	
	return final;
}

#endif
