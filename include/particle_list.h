#ifndef __PARTICLE_LIST_H__
#define __PARTICLE_LIST_H__

//#include <cutil_inline.h>
#include <boost/cstdint.hpp>
#include <vector_types.h>

namespace tweezers
{
	class particle_list {
    public:
      float3* pos;
      float3* vel;
      float3* acc;
      bool* active;
      int* num_particles;

    private:

	};
}

#endif

