#include "double3.h"

namespace tweezers {



double3S& double3S::operator=(const double3S &rhs) {
  
  if(this != &rhs) {
    x = rhs.x; y = rhs.y; z = rhs.z;
  }
  
  return *this;
}

double3S& double3S::operator+=(const double3S &rhs) {
  x += rhs.x; y += rhs.y; z += rhs.z;
  return *this;
}

const double3S double3S::operator+(const double3S &d) const {
  double3S res = *this;
  res += d;
  return res;
}

double3S& double3S::operator*=(const double &c) {
  x *= c; y *= c; z *= c;
  return *this;
}

const double3S double3S::operator*(const double &c) {
  double3S res = make_double3S(this->x, this->y, this->z);
  res *= c;
  return res;
}

}
