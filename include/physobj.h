#ifndef __TWEEZERS_PHYSOBJ_H__
#define __TWEEZERS_PHYSOBJ_H__

#include <iostream>

#include <vector_types.h>

#include "vector3d.h"

using namespace std;

namespace tweezers
{
	class physobj
	{
		public:
			physobj(float x, float y, float z);	
			physobj(float x, float y, float z, float vx, float vy, float vz, float ax, float ay, float az);
			~physobj();
	
			float3 pos;
			float3 oldpos;
			float3 vel;
			float3 acc;
	
			friend ostream &operator << (ostream &outs, physobj &p);
	
		private:
	};
}

#endif

