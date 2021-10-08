#include "laser.h"

namespace tweezers
{
  laser::laser() :
    physobj(0.0, 0.0, 0.0) {}

	laser::laser(double x, double y, double z, double vx, double vy, double vz, double ax, double ay, double az)
		: physobj(x,y,z,vx,vy,vz,ax,ay,az)
	{
		
	}

	laser::~laser()
	{
	
	}
}
