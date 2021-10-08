#ifndef __LASER_LIST_H__
#define __LASER_LIST_H__

//#include <cutil_inline.h>
#include <boost/cstdint.hpp>
#include <vector_types.h>

namespace tweezers
{
	class laser_list {
    public:
      double3S* pos;
      double3S* vel;
      int* num_lasers;

    private:

	};
}

#endif

