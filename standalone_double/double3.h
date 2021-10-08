#ifndef __double3S_H__
#define __double3S_H__

namespace  tweezers {
  class double3S {
    
  public:
    double x;
    double y;
    double z;
    
    friend double3S make_double3S(double _x, double _y, double _z);
    
    double3S& operator=(const double3S &rhs);

    double3S& operator+=(const double3S &rhs);
    const double3S operator+(const double3S &d) const; 

    double3S& operator*=(const double &c);
    const double3S operator*(const double &c);
  };

  inline double3S make_double3S (double _x, double _y, double _z) 
  { 
	  double3S d; 
	  d.x=_x;
	  d.y=_y;
	  d.z=_z;
	  return d; 
  }
}



#endif
