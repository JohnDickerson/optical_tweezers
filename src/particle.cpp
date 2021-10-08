#include "particle.h"

namespace tweezers
{
	particle::particle(float x, float y, float z, float vx, float vy, float vz, float ax, float ay, float az)
      : physobj(x,y,z,vx,vy,vz,ax,ay,az), active(true)
	{}
	
	particle::~particle()
	{

	}
	
}
