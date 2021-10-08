#include "velocity.h"
#include "io.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>



using namespace blitz;

Velocity OT;
int trapped = 10;

double NA = 1.4;
double riOil = 1.51;

double angle, aspectRatio;

double computeTimeInterval (double radius)			
{
	double viscosity = 1.002e-3;
	double pi = 3.1415;
	double density = 2000;

	angle = asin (NA / riOil);
	aspectRatio = tan (angle);

	double mass = 4 * pi * density * radius * radius * radius / 3;
	double drag = 6 * pi * viscosity * radius;

	int x = (int) ((mass / drag) * 1e7);
	int y[] = {1, 4, 5, 10, 20, 25, 40, 50, 80, 100, 200, 250, 400, 500, 1000};
	int quo = 0;

	for (int i = 0; i < 14; ++i)
	{
		if ((x >= y[i]) && (x < y[i+1]))
			quo = y[i];
	}

	double interval = (x - (x % quo)) * 1e-7;				// Lowest
                                                            // possible
                                                            // value
                                                            // is 0.1
                                                            // microsecond 
    //std::cerr << interval  << "\n";
	return interval;
}

bool checkTrapping (double *o, double *n, int num)
{
	if ((num >= 100) && (abs(o[0] - n[0]) < 5e-7) && (abs(o[1] - n[1]) < 5e-7) && ((n[2] - o[2]) >= 0) && ((n[2] - o[2]) < 5e-7))	
	{
		trapped = 0;
		return true;
	}

	else
		return false;
}

bool endCalc (double *f, double *s, double radius, int num)
{
	double z = abs (f[2] - s[2]);
	double r = z * aspectRatio + 0.4e-6;				// 0.4 micron beam focal waist radius
		
	double dist = sqrt ((f[0] - s[0]) * (f[0] - s[0]) + (f[1] - s[1]) * (f[1] - s[1]));

	if (num >= 100)
	{
		if (z > 3e-6)									// Suppress termination criterion when particle near beam focus
		{
			if (dist > (r + radius))					
			{
				trapped = 1;
				return true;
			}
		}
	}

	if ((num >= 100) && ((s[2] - f[2]) >= 15e-6))		// Particle going to settle at the bottom
	{
		trapped = 2;
		return true;
	}
	
	if (trapped == 10)
		return false;
}

int main (int argc, char *argv[])
{
	FILE *f;
		
	double *laserVel = new double [3];
	double *sphereVel = new double [3];
	
	double *focalPosn = new double [3];
	double *spherePosn = new double [3];

	double *oldCoord = new double [3];
	double *OTforces = new double [3];
	double *sphereAccn = new double [3];

	double timeSim = 0.0;

	for (int i = 0; i < 3; i++)
	{
		focalPosn[i] = atof(argv[1+i]) * 1e-6;
		spherePosn[i] = atof(argv[4+i]) / 10000 * 1e-6;

		laserVel[i] = atof(argv[7+i]) / 10000 * 1e-3;
		sphereVel[i] = atof(argv[10+i]) / 1000 * 1e-3;
		
		sphereAccn[i] = 0.0;
	}

   
    //std::cerr << "laser@ : " << focalPosn[0] 
    //          << ", " << focalPosn[1] 
    //          << ", " << focalPosn[2] << "\n"; 


    //std::cerr << "laser velocity : " << laserVel[0] 
    //          << ", " << laserVel[1] 
    //          << ", " << laserVel[2] << "\n"; 

    //std::cerr << "Sphere@ : " << spherePosn[0] 
    //          << ", " << spherePosn[1] 
    //          << ", " << spherePosn[2] << "\n"; 

    //std::cerr << "Sphere velocity : " << sphereVel[0] 
    //          << ", " << sphereVel[1] 
    //          << ", " << sphereVel[2] << "\n"; 

	double radius = atof(argv[13]) / 1000 * 1e-6;
	double interval = computeTimeInterval (radius);
			    
    // std::cerr << "Sphere radius : " << radius << "\n";

	f = fopen(argv[14], "a");
	int iterNo = 0;

    double start_y = spherePosn[1] * 1e6;
    double start_z = spherePosn[2] * 1e6;

    double end_y = (atof(argv[15]) / 10000 * 1e-6) * 1e6;
    double end_z = (atof(argv[16]) / 10000 * 1e-6) * 1e6;

    double step_size = (atof(argv[17]) / 10000 * 1e-6) * 1e6;

    std::string force_in(argv[18]);
    Array< TinyVector<double,2>, 2> force_table;
    load_force_data(force_in, 1e-6, 1e-12, force_table); 

    std::string prob_out(argv[19]);
    std::ofstream prob_file;
    prob_file.open(prob_out.c_str());

    std::cerr << "start_y = " << start_y << ", end_y = " << end_y
              << ", start_z = " << start_z << ", end_z = " << end_z 
              << ", step_size = " << step_size << "\n";
    
    int steps_per_micron = 1.0 / step_size;
    cerr << "steps_per_micron = "  << steps_per_micron << "\n";

    for (double y = start_y; y < end_y; y += step_size) {
      for (double z = start_z; z < end_z; z += step_size) {

        double use_y = static_cast<double>(static_cast<int>(y * steps_per_micron)) / (steps_per_micron);
        double use_z = static_cast<double>(static_cast<int>(z * steps_per_micron)) / (steps_per_micron);
        cerr << "y = " << y << ", z = " << z << "\n"; 


        float trapped_count = 0.0;

        size_t max_pc_count = 1;

        for (size_t pc = 0; pc < max_pc_count; ++pc) {
          // New sphere pos
          spherePosn[0] = 0.0;
          spherePosn[1] = (use_y * 1e-6);
          spherePosn[2] = (use_z * 1e-6);

          // Zero out Vel, Acc, iter
          for (int j = 0; j < 3; ++j) {
            sphereVel[j] = 0.0;
            sphereAccn[j] = 0.0;
            iterNo = 0;
          }
        
          trapped = 10;

          while ((! checkTrapping (focalPosn, spherePosn, iterNo)) && (! endCalc (focalPosn, spherePosn, radius, iterNo)) && (iterNo <= 0)) { // 1000000)) {

            //            if(iterNo % 10000 == 0) {
            //              printf("\n\n ITERATION NUMBER: %d \n\n", iterNo);
            //            }

            for (int j = 0; j < 3; ++j) {
              focalPosn[j] += laserVel[j] * interval;
            }

            //        ofile << spherePosn[1]*1e6 << " " << spherePosn[2]*1e6;
            
            cout << use_y << " " << use_z << " ";

            sphereVel = OT.calcVelInsideTrap (focalPosn, spherePosn, sphereVel, sphereAccn, radius, interval, timeSim, force_table);
            spherePosn = OT.calcPosnInsideTrap ();

            sphereAccn = OT.calcAccnInsideTrap ();
            OTforces = OT.calcForceInsideTrap ();

            timeSim = OT.calcTimeInsideTrap ();

            //
            //  }
            //  }


   		
            //            if (iterNo % 100 == 0)
            //              fprintf(f, "%1.12lf %1.12lf %1.12lf\n", spherePosn[0], spherePosn[1], spherePosn[2]);

            ++iterNo;

          } //end WHILE 
          
          trapped_count += (trapped == 0);
        }
        cerr << y*1e-6 << " " << z*1e-6 << " " << static_cast<double>(trapped_count) / max_pc_count << "\n";
        prob_file << y*1e-6 << " " << z*1e-6 << " " << static_cast<double>(trapped_count) / max_pc_count << "\n";
        prob_file << std::flush;
      }
    }
    prob_file.close();

	fprintf (f, "trap %d iteration %d\n", trapped, iterNo);
	fclose (f);								

	exit (0);

	if (focalPosn != NULL)
		delete focalPosn;

	if (spherePosn != NULL)
		delete spherePosn;

	if (oldCoord != NULL)
		delete oldCoord;

	if (laserVel != NULL)
		delete laserVel;

	if (sphereVel != NULL)
		delete sphereVel;

    return 0;
}

