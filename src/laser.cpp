#include "laser.h"

namespace tweezers
{
  laser::laser() :
    physobj(0.0, 0.0, 0.0) {}

	laser::laser(float x, float y, float z, float vx, float vy, float vz, float ax, float ay, float az)
		: physobj(x,y,z,vx,vy,vz,ax,ay,az)
	{
		
	}

	laser::~laser()
	{
	
	}
}
