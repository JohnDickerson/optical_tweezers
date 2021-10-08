#include "physobj.h"

namespace tweezers
{
	physobj::physobj(double x, double y, double z)
	{
		physobj(x,y,z,0,0,0,0,0,0);
	}
	
	physobj::physobj(double x, double y, double z, double vx, double vy, double vz, double ax, double ay, double az)
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
