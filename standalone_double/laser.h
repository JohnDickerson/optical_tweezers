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
    laser(double x, double y, double z, double vx=0.0, double vy=0.0, double vz=0.0, double ax=0.0, double ay=0.0, double az=0.0);
    ~laser();


  private:

  };
}

#endif

