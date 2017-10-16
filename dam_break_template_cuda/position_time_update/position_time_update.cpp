#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\basic_func\basic_func.h"
#include <vector>
using namespace std;

void position_time_update_schm1(vector<PARTICLE> & particles, double dt)
{
	// Position update
	for(int i = 0; i < particles.size(); ++ i)
	{
		if(particles[i].id == 0 || particles[i].id == 1 || particles[i].id == 2 || particles[i].id == 3) continue; 

		// Coordination X update
		particles[i].coorX.val[2] = particles[i].coorX.val[1]; 
		particles[i].coorX.val[1] = particles[i].coorX.val[0]; 
		particles[i].coorX.val[0] += dt * particles[i].velX.val[0]; 
		

		// Coordination Y update
		particles[i].coorY.val[2] = particles[i].coorY.val[1]; 
		particles[i].coorY.val[1] = particles[i].coorY.val[0]; 
		particles[i].coorY.val[0] += dt * particles[i].velY.val[0]; 
		

		// Coordination Z update
		particles[i].coorZ.val[2] = particles[i].coorZ.val[1]; 
		particles[i].coorZ.val[1] = particles[i].coorZ.val[0]; 
		particles[i].coorZ.val[0] += dt * particles[i].velZ.val[0];

		// * Box constrain
		if(BOX_CONTAIN == 1)
		{
			if(particles[i].coorX.val[0] < BOX_X_MIN) particles[i].coorX.val[0] = BOX_X_MIN;
			if(particles[i].coorX.val[0] > BOX_X_MAX) particles[i].coorX.val[0] = BOX_X_MAX;

			if(particles[i].coorY.val[0] < BOX_Y_MIN) particles[i].coorY.val[0] = BOX_Y_MIN;
			if(particles[i].coorY.val[0] > BOX_Y_MAX) particles[i].coorY.val[0] = BOX_Y_MAX;

			if(particles[i].coorZ.val[0] < BOX_Z_MIN) particles[i].coorZ.val[0] = BOX_Z_MIN; 
			if(particles[i].coorZ.val[0] > BOX_Z_MAX) particles[i].coorZ.val[0] = BOX_Z_MAX;

			if(CASE_DIM == 0) particles[i].coorZ.val[0] = 0.0;
		}	
	}

	// Time update
	for(int i = 0; i < particles.size(); ++ i)
	{
		particles[i].time.val[2] = particles[i].time.val[1]; 
		particles[i].time.val[1] = particles[i].time.val[0];
		particles[i].time.val[0] += dt; 
	}
}
