#include "RNG.hpp"

RNG::RNG() {
	// Default behavior is to generate 
	// a unifrom random number between 0 and 1
	m_init_flag=false;
	iset=0;
	init(0, 1, UNIFORM); 
}

RNG::~RNG() {};

RNG::RNG(double _mean, double _sd, int _distribution) {
	m_init_flag=false;
	iset=0; 
	init(_mean, _sd, _distribution);
}

int RNG::init(double _mean, double _sd, int _distribution) {
	m_mean=_mean;
	m_sd=_sd;
	m_distribution=_distribution;
	
	m_randptr = (m_distribution==NORMAL) ? &RNG::normal_rand : &RNG::uniform_rand;
	
#ifdef WIN32
	QueryPerformanceFrequency(&timerFrequency);		// timer frequency for seed
#endif

	m_init_flag=true;

	// seed the RNG
	seed();
	
	return 0;
}

int RNG::seed() {
	double time;
	if(!m_init_flag) {
		std::cerr << "Error: Random number generator not initialized.\n";
		return -1;
	} 
#ifdef WIN32	
	QueryPerformanceCounter(&timeCount);
	// Calculate elapsed time 
	time = ((double)timeCount.QuadPart)
		/((double)timerFrequency.QuadPart);
#else
	gettimeofday(&timeCount, &tz);
	time = (double) (timeCount.tv_sec+timeCount.tv_usec);
#endif
	mRandSeed = (int) time*-1;
	
	return 0;
}

/*int RNG::seed(int s) {
	// make sure the seed is always a negative number.
	mRandSeed = (int) -1 * abs(s);
}*/


/*
Long period (>2e18) random number generator of L'Ecuyer with Bayes-Durham shuffle 
and added safdguards.  Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the 
endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter 
idum between successive deviates in a sequence. RNMX should approximate the largest  
floating value that is less than 1.
*/
double RNG::uniform_rand(int &idum) {
	const int IM1=2147483563, IM2=2147483399;
	const int IA1=40014, IA2=40692, IQ1=53668, IQ2=52774;
	const int IR1=12211, IR2=3791, NTAB=32, IMM1=IM1-1;
	const int NDIV = 1 + IMM1/NTAB;
	const double EPS=3.0e-16, RNMX=1.0-EPS, AM=1.0/(double)IM1;
	static int idum2=123456789, iy=0;
	static int iv[NTAB];
	int j, k;
	double temp;
	
	if (idum <= 0) {						// Initialize
		idum=(idum==0 ? 1 : -idum);			// Be sure to prevent idum = 0
		idum2=idum;
		for (j=NTAB+7; j>=0; j--) {			// Load the shuffle table (after 8 warm-ups)
			k=idum/IQ1;
			idum=IA1*(idum-k*IQ1)-k*IR1;
			if (idum < 0) idum +=IM1;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=idum/IQ1;								// Start here when not initializing
	idum=IA1*(idum-k*IQ1)-k*IR1;			// Compute idum=(IA1*idum) % IM1 without
	if (idum < 0) idum +=IM1;				//		overflows by Schrage's method.
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;			// Compute idum2=(IA2*idum2) % IM2 likewise
	if (idum2 < 0) idum2 +=IM2;
	j=iy/NDIV;								// Will be in the range 0..NTAB-1
	iy=iv[j]-idum2;							// Here idum is shuffled. idum and idum2 are 
	iv[j] = idum;							// 		combined to generate the output.
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX ) 
		return RNMX;						// Because users don't expect endpoint values
	else 
		return temp;
}

/*	Returns a normally distributed deviate with zero mean and unit variance, 
	using ran1(idum) as the source of uniform deviates.
*/
double RNG::normal_rand(int &idum) { 
	double fac,rsq,v1,v2; 
	
	if (idum < 0) 
		iset=0; 									// Reinitialize. 
				
		if (iset == 0) {		 					// We don t have an extra deviate handy, so 
			do { 
				v1=2.0*uniform_rand(idum)-1.0;		// pick two uniform numbers in the square extending from -1 to +1 in each direction, 
				v2=2.0*uniform_rand(idum)-1.0;
				//rsq=(v1*v1-m_mean*m_mean)+(v2*v2-m_mean*m_mean);
				rsq=v1*v1+v2*v2; 					//see if they are in the unit circle,
			} while (rsq >= 1.0 || rsq == 0.0); 	//		and if they are not, try again.
			fac=sqrt(-2.0*log(rsq)/rsq); 
			
			/* Now make the Box-Muller transformation 
			 * to get two normal deviates. Return one 
			 * and save the other for next time.
			*/
			gset=v1*fac; 
			iset=1; 								//Set flag. 
			retValue = v2*fac; 
		} else { 									// We have an extra deviate handy, 
			iset=0; 								//		so unset the flag, 
			retValue = gset; 						// and return it. 
		} 
		return retValue * m_sd + m_mean;			// First Multiply by the sd 
													// 		shift by the mean before returning 
}

double RNG::rand() {
//	std::cout << mRandSeed << std::endl;
	return (this->*m_randptr) (mRandSeed);
}
