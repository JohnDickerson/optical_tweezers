#ifndef __PARTICLE_LIST_H__
#define __PARTICLE_LIST_H__

#include <boost/cstdint.hpp>
#include "constants.h"

namespace tweezers
{
	class particle_list {
    public:
      double3S* pos;
      double3S* vel;
      double3S* acc;
      bool* active;
      int* num_particles;

    private:

	};
}

#endif

