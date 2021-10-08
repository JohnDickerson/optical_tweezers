#include <memory>
#include "velocity.h"
#include <blitz/blitz.h>
#include <cerrno>
#include <cstring>
#include <cstdio>

Velocity::Velocity() 
{
	integrationCounter = 41;			// Accuracy can be enhanced by increasing this no.
	
	pi = 3.1415;
	c = 3e8;

	riParticle = 1.46;						// Amorphous silica
	riOil = 1.51;
	riWater = 1.33;

	power = 0.1;
	a = 1.0;
	lambda = 0.532e-6;

	NA = 1.4;							// Oil immersion objective
	f = 3.33e-3;
	absT = 293.15;						// Room temperature

	viscosity = 1.002e-3;
	densityParticle = 2000;
	densityWater = 1000;
	g = 9.81;

	dF = NULL;
	F = NULL;
	temp = NULL;

	vel = NULL;
	posn = NULL;
	accn = NULL;

	allocateArrays ();
}

Velocity::~Velocity() 
{
	if (dF != NULL)
		delete [] dF;

	if (F != NULL)
		delete [] F;

	if (vel != NULL)
		delete [] vel;

	if (posn != NULL)
		delete [] posn;

	if (accn != NULL)
		delete [] accn;

	if (temp != NULL)
		for (int i = 0; i < this->integrationCounter; ++i)
			delete [] temp[i];
	delete [] temp;
}

void Velocity::allocateArrays () 
{
	dF = new double [3];
	F = new double [3];

	vel = new double [3];
	posn = new double [3];
	accn = new double [3];
	
	temp = new double *[this->integrationCounter];
	for (int i = 0; i < this->integrationCounter; ++i)
		temp[i] = new double [3 * this->integrationCounter];

	for (int k = 0; k < 3; ++k)
	{
		this->dF[k] = 0;
		this->F[k] = 0;
	}
}

double* Velocity::calcVelInsideTrap (double *a, double *b, 
                                     double *sv, double *sa, 
                                     double radius, double interval, double t,
                                     const Array< blitz::TinyVector<double, 2>, 2 >& force_array)

{


  for (int i = 0; i < this->integrationCounter; ++i)
    for (int j = 0; j < 3 * this->integrationCounter; ++j)
      temp[i][j] = 0.0;

  for (size_t j = 0; j < 3; ++j) {
	dF[j] = 0.0;
	F[j] = 0.0;
    //	vel[j] = 0.0;
    //	posn[j] = 0.0;
    //	accn[j] = 0.0;
  }


	sParticleparams3D<double> simOutput;	
	sPoint3D<double> Fext;	

	this->mass = (4.0f/3.0f) * this->pi * this->densityParticle * radius * radius * radius;
	this->massEff = (4.0f/3.0f) * this->pi * (this->densityParticle - this->densityWater) * radius * radius * radius;
	this->drag = 6 * this->pi * this->viscosity * radius;

	simOutput.t = t;

	simOutput.r.x = b[0];
	simOutput.r.y = b[1];
	simOutput.r.z = b[2];
		
	simOutput.v.x = sv[0];
	simOutput.v.y = sv[1];
	simOutput.v.z = sv[2];

	simOutput.a.x = sa[0];
	simOutput.a.y = sa[1];
	simOutput.a.z = sa[2];

    
    calculateForce (a, b, radius);
    //lookupForce(a, b, force_array);

	Fext.x = F[0];
	Fext.y = F[1];
	Fext.z = F[2] + this->massEff * this->g;			// Buoyancy effect included

	auto_ptr <cpuSphericalParticleSim> cpuComputeObject (new cpuSphericalParticleSim(
			drag,
			absT,
			mass,
			interval
		)
	);

	simOutput = cpuComputeObject->CPUPass<double> (simOutput, Fext);

	timeSim = simOutput.t;

	posn[0] = simOutput.r.x;
	posn[1] = simOutput.r.y;
	posn[2] = simOutput.r.z;

	vel[0] = simOutput.v.x;
	vel[1] = simOutput.v.y;
	vel[2] = simOutput.v.z;

	accn[0] = simOutput.a.x;
	accn[1] = simOutput.a.y;
	accn[2] = simOutput.a.z;
    
    // cerr << "posn = <" << posn[0] << ", " << posn[1] << ", " << posn[2] << "> \n";
    // cerr << "vel = <" << vel[0] << ", " << vel[1] << ", " << vel[2] << "> \n";
    // cerr << "accn = <" << accn[0] << ", " << accn[1] << ", " << accn[2] << "> \n";

    // cerr << "this->F = <" << this->F[0] << ", " << this->F[1] << ", " << this->F[2] << "> \n";
    // cerr << "Fext = <" << Fext.x << ", " << Fext.y << ", " << Fext.z << "> \n";

	return vel;
}

double *Velocity::calcPosnInsideTrap ()
{
	return posn;
}

double *Velocity::calcAccnInsideTrap ()
{
	return accn;
}

double *Velocity::calcForceInsideTrap ()
{
	return F;
}

double Velocity::calcTimeInsideTrap ()
{
	return timeSim;
}

bool bi_interp_lookup ( const Array< blitz::TinyVector<double, 2>, 2 >& force_array,
                        double r, double z, TinyVector<double, 2>& interp_res ) {

  double rmin = 0.0f;
  double rmax = 20.0f;
  double zmin = -20.0f;
  double zmax = 8.0f;
    
  size_t nysamp = force_array.extent(1);
  size_t nzsamp = force_array.extent(0);

  int steps_per_micron = 4;
  
  double ind_r = ((r - rmin) / (rmax - rmin)) * (nysamp-1);
  double ind_z = ((z - zmin) / (zmax - zmin)) * (nzsamp-1);
    
  double low_r = static_cast<double>(static_cast<int>(r * steps_per_micron)) / (steps_per_micron);
  double low_z = static_cast<double>(static_cast<int>(z * steps_per_micron)) / (steps_per_micron);
  
  if (r > rmax || r < rmin || z > zmax || z < zmin) {
    interp_res[0] = 0.0;
    interp_res[1] = 0.0;
    return false;
  }

  int r_low = floor(ind_r);
  int r_hi = ceil(ind_r);
  int z_low = floor(ind_z);
  int z_hi = ceil(ind_z);

  double pr_low = ((static_cast<double>(r_low)/(nysamp-1)) * (rmax-rmin)) + rmin;
  double pr_hi = ((static_cast<double>(r_hi)/(nysamp-1)) * (rmax-rmin)) + rmin;
  double pz_low = ((static_cast<double>(z_low)/(nzsamp-1)) * (zmax-zmin)) + zmin;
  double pz_hi = ((static_cast<double>(z_hi)/(nzsamp-1)) * (zmax-zmin)) + zmin;

  cerr << "\n\n" 
       << "pr_low = " << pr_low << " "
       << "pr_hi = " << pr_hi << " " 
       << "pz_low = " << pz_low << " " 
       << "pz_hi = " << pz_hi << "\n";

  double ir_low = ((pr_low - rmin) / (rmax - rmin)) * (nysamp-1);
  double ir_hi = ((pr_hi - rmin) / (rmax - rmin)) * (nysamp-1);
  double iz_low = ((pz_low - zmin) / (zmax - zmin)) * (nzsamp-1);
  double iz_hi = ((pz_hi - zmin) / (zmax - zmin)) * (nzsamp-1);

  cerr << "ir_low = " << ir_low << " "
       << "ir_hi = " << ir_hi << " " 
       << "iz_low = " << iz_low << " " 
       << "iz_hi = " << iz_hi << "\n\n";


  TinyVector<double, 2> interp_rhi(0.0);
  TinyVector<double, 2> interp_rlow(0.0);
  TinyVector<double, 2> interp_z(0.0);
  TinyVector<double, 2> temp(0.0);

  double i0 = static_cast<double>(pr_hi - r) / (pr_hi - pr_low);
  double i1 = static_cast<double>(r - pr_low) / (pr_hi - pr_low);

  if ( (r_hi == r_low) ) {
    i0 = 1.0;
    i1 = 0.0;
  }

  double j0 = static_cast<double>(pz_hi - z) / (pz_hi - pz_low);
  double j1 = static_cast<double>(z - pz_low) / (pz_hi - pz_low);

  if ( (z_hi == z_low) ) {
    j0 = 1.0;
    j1 = 0.0;
  }

  cerr << "ind_r = " << ind_r << ", ind_z = " << ind_z 
       << ", r_low = " << r_low << ", r_hi = " << r_hi 
       << ", z_low = " << z_low << ", z_hi = " << z_hi << "\n";
  cerr << "i0 = " << i0 << ", i1 = " << i1 << ", j0 = " << j0 << ", j1 = " << j1 << "\n";


  interp_rlow = force_array(z_low, r_low);
  interp_rlow *= i0;
  temp = force_array(z_low, r_hi);
  temp *= i1;
  interp_rlow += temp;
  temp[0] = temp[1] = 0.0;

  cerr << "interp_rlow = " << interp_rlow << " = "  
       << i0 << " *  f00 = " << force_array(z_low, r_low) << " "
       << i1  << " * f01 = " << force_array(z_low, r_hi) << "\n";
  
  interp_rhi = force_array(z_hi, r_low);
  interp_rhi *= i0;
  temp = force_array(z_hi, r_hi);
  temp *= i1;
  interp_rhi += temp;
  temp[0] = temp[1] = 0.0;

  cerr << "interp_hi = " << interp_rhi << " = "  
       << i0 << " *  f10 = " << force_array(z_hi, r_low) << " "
       << i1  << " * f11 = " << force_array(z_hi, r_hi) << "\n";

  interp_z = interp_rlow;
  interp_z *= j0;
  temp = interp_rhi;
  temp *= j1;
  interp_z += temp;
  temp[0] = temp[1] = 0.0;

  interp_res[0] = interp_z[0] * 1e-12;
  interp_res[1] = interp_z[1] * 1e-12;


  cerr << "[ " << r << ", " << z << "] ==> "
       << interp_res << "\n";


  return true;
}


void Velocity::lookupForce( double *a, double *b,
                            const Array< blitz::TinyVector<double, 2>, 2 >& force_array ) {

	this->sDash = sqrt ((b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]));

	this->sDDash = b[2] - a[2];						// Assuming upward beam propagation and right-hand coordinate system
    
	if (this->sDDash == 0)
		this->sDDash = 0.0005 * 1e-6;
	
	if (this->NA > this->riOil)
		this->NA = this->riOil;

	this->phiMax = asin (this->NA / this->riOil);
	this->rmax = tan (this->phiMax) * this->f;						// Radial integration limit 
	this->omega0 = this->a * this->rmax;


    TinyVector<double, 2> interp_res;
    bi_interp_lookup( force_array, this->sDash * 1e6, this->sDDash * 1e6, interp_res );

    F[1] = interp_res[0];
    F[2] = interp_res[1];
    cerr << "LOOKUP FORCES = " << F[1] << ", " << F[2] << "\n";

    F[0] = 0.0;
    F[1] = 0.0;
    F[2] = 0.0;

    calculateForceAshkin (2.5e-6);

    cerr << "CALC FORCES = " << F[1] << ", " << F[2] << "\n";


	double angle = atan2 ((b[1] - a[1]), (b[0] - a[0]));
	
	F[0] = F[1] * cos (angle);
	F[1] = F[1] * sin (angle);

	if (abs(sDash) < 1e-18)
	{
		F[0] = 0;
		F[1] = 0;
	}

	for (int i = 0; i < 3; ++i)
	{
		if (fabs(F[i]) < 1e-17)
			F[i] = 0;
	}

}

void Velocity::calculateForce (double *a, double *b, double radius)
{
	this->sDash = sqrt ((b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]));

	this->sDDash = b[2] - a[2];						// Assuming upward beam propagation and right-hand coordinate system
    
	if (this->sDDash == 0)
		this->sDDash = 0.0005 * 1e-6;
	
	if (this->NA > this->riOil)
		this->NA = this->riOil;

	this->phiMax = asin (this->NA / this->riOil);
	this->rmax = tan (this->phiMax) * this->f;						// Radial integration limit 
	this->omega0 = this->a * this->rmax;

	calculateForceAshkin (radius);
	
	double angle = atan2 ((b[1] - a[1]), (b[0] - a[0]));

    cout << " " << F[1]*1e12 << " " << F[2]*1e12 << "\n";
	
	F[0] = F[1] * cos (angle);
	F[1] = F[1] * sin (angle);

	if (abs(sDash) < 1e-18)
	{
		F[0] = 0;
		F[1] = 0;
	}

	for (int i = 0; i < 3; ++i)
	{
		if (fabs(F[i]) < 1e-17)
			F[i] = 0;
	}

}

void Velocity::calculateForceAshkin (double radius)
{

	double sum[3];

	double h1 = (this->pi - 0) / (this->integrationCounter - 1);
	double h2 = (this->rmax - 0) / (this->integrationCounter - 1);
	
	for (int i = 0; i < this->integrationCounter; ++i)
	{
		this->beta = 0 + i * pi / (this->integrationCounter - 1);

		for (int k = 0; k < 3; ++k)
			sum[k] = 0;

		for (int j = 0; j < this->integrationCounter; ++j)
		{
			this->ra = 0.005 * 1e-6 + j * this->rmax / (this->integrationCounter - 1);
            
			calculateAngles (radius);
		
			calculateSingleRay ();
			integrate (i, j);					// Using Romberg's integration rule: Source: http://mathworld.wolfram.com/Newton-CotesFormulas.html

			if (j == 0 || j == (this->integrationCounter - 1))
			{
				for (int k = 1; k < 3; ++k)
					sum[k] += h2 / 3 * this->temp[i][3 * j + k];
			}

			else
			{
				if (j % 2 == 1)
				{
					for (int k = 1; k < 3; ++k)
						sum[k] += h2 / 3 * 4 * this->temp[i][3 * j + k];
				}

				if (j % 2 == 0)
				{
					for (int k = 1; k < 3; ++k)
						sum[k] += h2 / 3 * 2 * this->temp[i][3 * j + k];
				}
			}
		}								// End of inner summation (integration)
		
		if (i == 0 || i == (this->integrationCounter - 1))
		{
			for (int k = 1; k < 3; ++k)
				this->F[k] += h1 / 3 * sum[k];
		}

		else
		{
			if (i % 2 == 1)
			{
				for (int k = 1; k < 3; ++k)
					this->F[k] += h1 / 3 * 4 * sum[k];
			}

			if (i % 2 == 0)
			{
				for (int k = 1; k < 3; ++k)
					this->F[k] += h1 / 3 * 2 * sum[k];
			}
		}
	}								// End of double summation
                                    // (integration)
    
}

void Velocity::calculateAngles (double radius)
{
	double x = this->ra * sin (this->beta);
	double y = this->ra * cos (this->beta);
	double z = - this->f;

	double norm = sqrt (x * x + y * y + z * z);

	this->alpha = acos ((x * x + y * y) / (norm * this->ra));

	this->betaD = atan2 (this->sDash * sin (this->beta), (this->sDash * cos (this->beta) + this->sDDash / tan (this->alpha)));

	this->d = this->sDDash * cos (this->betaD) / tan (this->alpha) + this->sDash * cos (this->beta - this->betaD);

	this->gammaD = acos (cos (this->alpha) * cos (this->betaD));

	double zz = tan (this->alpha) / tan (this->gammaD);

	if (zz < -1.0)
		zz = - 1.0;
	if (zz > 1.0)
		zz = 1.0;

	this->mu = acos (zz);

	double zzz = this->d * sin (this->gammaD) / radius;

	if (zzz < -1.0)
		zzz = - 1.0;
	if (zzz > 1.0)
		zzz = 1.0;

	this->theta = asin (zzz);					
	
	this->r = asin (this->riWater * sin (theta) / this->riParticle);				// Snell's law
}

void Velocity::calculateSingleRay ()
{
    this->Rs = (sin (this->r - this->theta) / sin (this->r + this->theta)) * (sin (this->r - this->theta) / sin (this->r + this->theta));

    this->Rp = (tan (this->r - this->theta) / tan (this->r + this->theta)) * (tan (this->r - this->theta) / tan (this->r + this->theta));

    this->Ts = 1 - this->Rs;
    this->Tp = 1 - this->Rp;

    this->R = (this->Rs + this->Rp) / 2;					// Circularly polarized light
    this->T = (this->Ts + this->Tp) / 2;

    this->qs = 1 + R * cos (2 * theta) - (T * T * (cos (2 * theta - 2 * r) + R * cos (2 * theta))) / (1 + R * R + 2 * R * cos (2 * r));

    this->qg = (R * sin (2 * theta) - (T * T * (sin (2 * theta - 2 * r) + R * sin (2 * theta))) / (1 + R * R + 2 * R * cos (2 * r)));

	this->dF[2] = qs * sin (alpha) + qg * cos (mu) * cos (alpha);

    this->dF[1] = -qs * cos (alpha) * cos (beta) + qg * cos (mu) * sin (alpha) * cos (beta) + qg * sin (mu) * sin (beta);
    this->dF[0] = -qs * cos (alpha) * sin (beta) + qg * cos (mu) * sin (alpha) * sin (beta) - qg * sin (mu) * cos (beta);
}

void Velocity::integrate (int i, int j)			
{
	double common = 4 * this->riWater * this->power / (this->c * this->pi * this->omega0 * this->omega0); 
	
	for (int k = 0; k < 3; ++k) {
		this->temp[i][3 * j + k] = common * this->ra * exp(-2 * this->ra * this->ra / (this->omega0 * this->omega0)) * this->dF[k];
	}
}
