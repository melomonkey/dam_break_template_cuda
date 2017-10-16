#ifndef struct_particle_h_
#define struct_particle_h_

#include <vector>
using namespace std; 

struct ATTRIBUTION
{
	double val[3]; 
}; 

struct PARTICLE
{
	int id; // type

	ATTRIBUTION mass; 
	ATTRIBUTION density; 
    ATTRIBUTION viscosity; 
	ATTRIBUTION pressure; 
	ATTRIBUTION temperature;
	ATTRIBUTION energy; 
	ATTRIBUTION velX; // velocityX
	ATTRIBUTION velY; // velocityY
	ATTRIBUTION velZ; // velocityZ
	ATTRIBUTION accX; // accelerationX
	ATTRIBUTION accY; // accelerationY
	ATTRIBUTION accZ; // accelerationZ
	ATTRIBUTION coorX; // coordinationX
	ATTRIBUTION coorY; // coordinationY
	ATTRIBUTION coorZ; // coordinationZ
	ATTRIBUTION smthR; // smoothedRadius
	ATTRIBUTION Radius; 
	ATTRIBUTION time; 

	double normal_vec[3]; 
	
};

#endif

