#include "physobj.h"

namespace tweezers
{
	physobj::physobj(float x, float y, float z)
	{
		physobj(x,y,z,0,0,0,0,0,0);
	}
	
	physobj::physobj(float x, float y, float z, float vx, float vy, float vz, float ax, float ay, float az)
	{
		pos.x = x; pos.y = y; pos.z = z;
		oldpos.x = x; oldpos.y = y; oldpos.z = z;
		vel.x = vx; vel.y = vy; vel.z = vz;
		acc.x = ax; acc.y = ay; acc.z = az;
	}
	
	physobj::~physobj()
	{
		//delete pos;
		//delete vel;
		//delete acc;
	}
	
	ostream &operator << (ostream &outs, physobj &p)
	{
		outs << "<" << p.pos.x << ", " << p.pos.y << ", " << p.pos.z << ">, ";
		outs << "<" << p.vel.x << ", " << p.vel.y << ", " << p.vel.z << ">, ";
		outs << "<" << p.acc.x << ", " << p.acc.y << ", " << p.acc.z << ">";
		return outs;
	}



}
