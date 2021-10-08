#ifndef _vector3d_H_
#define _vector3d_H_

#include <iostream>

using namespace std;

class vector3d
{
	public:
		vector3d();
		vector3d(double x, double y, double z);
		~vector3d();

		double dot(const vector3d &v) const;
		double angle(vector3d &v);
		double magnitude();
		vector3d cross(const vector3d &v) const;
		
		vector3d & normalize();
		vector3d & rotX(const double angle);
		vector3d & rotY(const double angle);
		vector3d & rotZ(const double angle);
		vector3d & rotate(const double angle, const vector3d &unit_v);
		
		vector3d operator + (const vector3d &v) const;
		vector3d operator - (const vector3d &v) const;
		vector3d operator * (const double &scalar) const;
		vector3d operator / (const double &scalar) const;
		vector3d operator - () const;
		
		vector3d & operator = (const vector3d &v);
		vector3d & operator += (const vector3d &v);
		vector3d & operator -= (const vector3d &v);
		vector3d & operator *= (const double &scalar);
		vector3d & operator /= (const double &scalar);
		
		double x;
		double y;
		double z;
		
		friend ostream &operator << (ostream &outs, vector3d &v);
	
	private:
		double mag;
		bool validMagnitude;
};

#endif


