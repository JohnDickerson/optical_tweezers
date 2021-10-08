#include "particle.h"

namespace tweezers
{
	particle::particle(double x, double y, double z, double vx, double vy, double vz, double ax, double ay, double az)
      : physobj(x,y,z,vx,vy,vz,ax,ay,az), active(true)
	{}
	
	particle::~particle()
	{

	}
	
}
