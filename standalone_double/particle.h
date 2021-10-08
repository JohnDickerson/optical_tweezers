#ifndef __TWEEZERS_PARTICLE_H__
#define __TWEEZERS_PARTICLE_H__

#include "vector3d.h"
#include "physobj.h"
#include "constants.h"

namespace tweezers
{
	class particle : public physobj
	{
		public:
			particle(double x, double y, double z, double vx=0.0f, double vy=0.0f, double vz=0.0f, double ax=0.0f, double ay=0.0f, double az=0.0f);
			~particle();
		
			bool isActive() { return active; }
			void setActive(bool active) { this->active = active; }
		private:
			bool active;
	};
}
#endif


