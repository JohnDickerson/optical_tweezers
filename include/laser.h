#ifndef __TWEEZERS_LASER_H__
#define __TWEEZERS_LASER_H__

#include "vector3d.h"
#include "physobj.h"

namespace tweezers
{
  class laser : public physobj
  {
  public:
    laser();
    laser(float x, float y, float z, float vx=0.0, float vy=0.0, float vz=0.0, float ax=0.0, float ay=0.0, float az=0.0);
    ~laser();


  private:

  };
}

#endif

