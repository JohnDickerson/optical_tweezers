#ifndef __TWEEZERS_PHYSOBJ_H__
#define __TWEEZERS_PHYSOBJ_H__

#include <iostream>

#include "constants.h"
#include "vector3d.h"

using namespace std;

namespace tweezers
{
	class physobj
	{
		public:
			physobj(double x, double y, double z);	
			physobj(double x, double y, double z, double vx, double vy, double vz, double ax, double ay, double az);
			~physobj();
	
			double3S pos;
			double3S oldpos;
			double3S vel;
			double3S acc;
	
			friend ostream &operator << (ostream &outs, physobj &p);
	
		private:
	};
}

#endif

