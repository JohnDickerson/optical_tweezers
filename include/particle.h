#ifndef __TWEEZERS_PARTICLE_H__
#define __TWEEZERS_PARTICLE_H__

#include "vector3d.h"
#include "physobj.h"

namespace tweezers
{
	class particle : public physobj
	{
		public:
			particle(float x, float y, float z, float vx=0.0f, float vy=0.0f, float vz=0.0f, float ax=0.0f, float ay=0.0f, float az=0.0f);
			~particle();
		
			bool isActive() { return active; }
			void setActive(bool active) { this->active = active; }
		private:
			bool active;
	};
}
#endif


