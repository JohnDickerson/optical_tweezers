/*
	Random number generator from Numerical Recipes
*/

#ifndef __RNG_HPP_
#define __RNG_HPP_

#include <cmath>
#include <iostream>

#include <fstream>

#ifdef WIN32
#include <windows.h>
#else					// Linux
#include <sys/time.h>
#endif

// Random Number generator constants
//#define IA 16807 
//#define IM 2147483647 
//#define AM (1.0/IM) 
//#define IQ 127773 
//#define IR 2836 
//#define NTAB 32 
//#define NDIV (1+(IM-1)/NTAB) 
//#define EPS 1.2e-7 
//#define RNMX (1.0-EPS) 

// Random Number distributions
#define UNIFORM 	0xafcb6801
#define NORMAL		0xafcb6802


class RNG {
private:
	double m_mean, m_sd;
	int m_distribution;
	bool m_init_flag;
	
	double gset;		// RNG variables that need to persist
	double retValue;
	int iset; 			// normal rand flag
	
	int mRandSeed;
#ifdef WIN32
	LARGE_INTEGER timerFrequency;
	LARGE_INTEGER timeCount;
#else
	struct timeval timeCount;
	struct timezone tz;
#endif	
	// typedef a pointer
	typedef double (RNG::*t_randfuncptr) (int &);

	// generate a uniform deviate random number between 0 and 1
	// adapted from code on pg 284, Numerical recipes in C++
	double uniform_rand(int &idum);					
	
	// transformation to Normal deviates
	// adapted from code on pg 293, Numerical recipes in C++
	double normal_rand(int &idum);
	
	// function pointer to the random number generator.
	// The adddress of this pointer is set at initialization.
	// Subsequent calls to RNG::rand() will be redirected to the function this var points to.
	t_randfuncptr m_randptr;
	
public:
	RNG();
	RNG(double _mean, double _sd, int _distribution);
	~RNG();
	int seed();
	//int seed(int s);
	int init(double _mean, double _sd, int _distribution);
	double rand();
};
#endif
