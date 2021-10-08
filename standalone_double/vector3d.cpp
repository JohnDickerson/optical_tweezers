#include <cmath>
#include "vector3d.h"

vector3d::vector3d()
 : x(0), y(0), z(0), mag(-1), validMagnitude(false) {}

vector3d::vector3d(double x, double y, double z)
{
	this->x = x;
	this->y = y;
	this->z = z;
	this->mag = -1;
	this->validMagnitude = false;
}

vector3d::~vector3d()
{
	
}


// Returns the dot product of two vectors
double vector3d::dot(const vector3d &v) const 
{
	return (this->x * v.x) + (this->y * v.y) + (this->z * v.z);
}

// Returns the cross product of two vectors
vector3d vector3d::cross(const vector3d &v) const
{
	vector3d result;
	result.x = this->y * v.z - this->z * v.y;
	result.y = this->z * v.x - this->x * v.z;
	result.z = this->x * v.y - this->y * v.x;
	return result;
}

// Calculates the angle (IN RADIANS) between two vectors
double vector3d::angle(vector3d &v)
{
	// a dot b = ||a|| ||b|| cos(theta)
	return acos( this->dot(v) / this->magnitude() * v.magnitude() );
}

// Updates and returns the magnitude of this vector
double vector3d::magnitude()
{
	// If we have invalidated our previously calculated magnitude, recalculate
	if(!validMagnitude)
	{
		mag = sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
		validMagnitude = true;
	}
	return mag;
}

// Normalizes the vector -- same direction, magnitude is now 1
vector3d & vector3d::normalize()
{
	if(!validMagnitude || mag != 1)
	{
		x /= magnitude();
		y /= magnitude();
		z /= magnitude();

		// We are guaranteed that the magnitude of a normalized vector is 1
		mag = 1;
		validMagnitude = true;
	}
	return *this;
}

// Rotates this vector about the x-axis by angle radians
vector3d & vector3d::rotX(const double angle)
{
	double cosAngle = cos(angle);
	double sinAngle = sin(angle);
	this->y = this->y*cosAngle - this->z*sinAngle;
	this->z = this->y*sinAngle + this->z*cosAngle;
	return *this;
}

// Rotates this vector about the y-axis by angle radians
vector3d & vector3d::rotY(const double angle)
{
	double cosAngle = cos(angle);
	double sinAngle = sin(angle);
	this->x = this->x*cosAngle - this->z*sinAngle;
	this->z = this->x*sinAngle + this->z*cosAngle;
	return *this;
}

// Rotates this vector about the z-axis by angle radians
vector3d & vector3d::rotZ(const double angle)
{
	double cosAngle = cos(angle);
	double sinAngle = sin(angle);
	this->x = this->x*cosAngle - this->y*sinAngle;
	this->y = this->x*sinAngle + this->y*cosAngle;
	return *this;
}

// Rotates this vector about an arbitrary unit vector by angle radians
vector3d & vector3d::rotate(const double angle, const vector3d &unit_v)
{
	double cosAngle = cos(angle);
	double sinAngle = sin(angle);
	
	this->x = this->x * (unit_v.x * unit_v.x + (1 - unit_v.x * unit_v.x) * cosAngle) + 
		this->y * (unit_v.x * unit_v.y * (1 - cosAngle) - unit_v.z * sinAngle) + 
		this->z * (unit_v.x * unit_v.z * (1 - cosAngle) + unit_v.y * sinAngle);
	
	this->y = this->x * (unit_v.x * unit_v.y * (1 - cosAngle) + unit_v.z * sinAngle) + 
		this->y * (unit_v.y * unit_v.y + (1 - unit_v.y * unit_v.y) * cosAngle) +
		this->z * (unit_v.y * unit_v.z * (1 - cosAngle) - unit_v.x * sinAngle);

	this->z = this->x * (unit_v.x * unit_v.z * (1 - cosAngle) - unit_v.y * sinAngle) +
		this->y * (unit_v.y * unit_v.z * (1 - cosAngle) + unit_v.x * sinAngle) + 
		this->z * (unit_v.z * unit_v.z + (1 - unit_v.z * unit_v.z) * cosAngle);
		
	return *this;
}

// Returns sum of two vectors
vector3d vector3d::operator + (const vector3d &v) const
{
	vector3d result;
	result.x = this->x + v.x;
	result.y = this->y + v.y;
	result.z = this->z + v.z;
	return result;
}

// Returns difference of two vectors
vector3d vector3d::operator - (const vector3d &v) const
{
	vector3d result;
	result.x = this->x - v.x;
	result.y = this->y - v.y;
	result.z = this->z - v.z;
	return result;
}

// Returns vector with each element scaled by scalar
vector3d vector3d::operator * (const double &scalar) const
{
	vector3d result;
	result.x = this->x * scalar;
	result.y = this->y * scalar;
	result.z = this->z * scalar;
	return result;
}

// Returns vector with each element divided by scalar
vector3d vector3d::operator / (const double &scalar) const
{
	vector3d result;
	result.x = this->x / scalar;
	result.y = this->y / scalar;
	result.z = this->z / scalar;
	return result;
}

// Returns vector with each element negated
vector3d vector3d::operator - () const
{
	vector3d result;
	result.x = -this->x;
	result.y = -this->y;
	result.z = -this->z;
	
	// A change in direction does not affect magnitude
	result.mag = this->mag;
	result.validMagnitude = this->validMagnitude;
	return result;
}
	
// Sets the underlying math of this vector to equal the passed-in vector
vector3d & vector3d::operator = (const vector3d &v)
{
	this->x = v.x;
	this->y = v.y;
	this->z = v.z;
	this->mag = v.mag;
	this->validMagnitude = v.validMagnitude;
	return *this;
}

// vector3d-adds a vector to this
vector3d & vector3d::operator += (const vector3d &v)
{
	this->validMagnitude = false;
	this->x += v.x;
	this->y += v.y;
	this->z += v.z;
	return *this;
}

// vector3d-subtracts a vector from this
vector3d & vector3d::operator -= (const vector3d &v)
{
	this->validMagnitude = false;
	this->x -= v.x;
	this->y -= v.y;
	this->z -= v.z;
	return *this;
}

// Multiplies each coordinate by scalar
vector3d & vector3d::operator *= (const double &scalar)
{
	this->validMagnitude = false;
	this->x *= scalar;
	this->y *= scalar;
	this->z *= scalar;
	return *this;
}

// Divides each coordinate by scalar
vector3d & vector3d::operator /= (const double &scalar)
{
	this->validMagnitude = false;
	this->x /= scalar;
	this->y /= scalar;
	this->z /= scalar;
	return *this;
}

// Prints the vector in the form <x, y, z>
ostream &operator << (ostream &outs, vector3d &v)
{
	outs << "<" << v.x << ", " << v.y << ", " << v.z << ">";
	return outs;
}
