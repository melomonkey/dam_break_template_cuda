#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\basic_func\basic_func.h"
#include <vector>
using namespace std;

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__global__ void _postion_time_update_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * neighbors_list, PARTICLE * particles, double dt, unsigned int offset)
{
	unsigned int pos = blockIdx.x + offset;

	if (pos > neigh_list_length - 1) return;

	unsigned int calcu_particle = neighbors_list[0 + neigh_list_width * pos] - 1;

	if (particles[calcu_particle].id == 0 || particles[calcu_particle].id == 1 || particles[calcu_particle].id == 2 || particles[calcu_particle].id == 3) return;

	// Coordination X update
	particles[calcu_particle].coorX.val[2] = particles[calcu_particle].coorX.val[1];
	particles[calcu_particle].coorX.val[1] = particles[calcu_particle].coorX.val[0];
	particles[calcu_particle].coorX.val[0] += dt * particles[calcu_particle].velX.val[0];


	// Coordination Y update
	particles[calcu_particle].coorY.val[2] = particles[calcu_particle].coorY.val[1];
	particles[calcu_particle].coorY.val[1] = particles[calcu_particle].coorY.val[0];
	particles[calcu_particle].coorY.val[0] += dt * particles[calcu_particle].velY.val[0];


	// Coordination Z update
	particles[calcu_particle].coorZ.val[2] = particles[calcu_particle].coorZ.val[1];
	particles[calcu_particle].coorZ.val[1] = particles[calcu_particle].coorZ.val[0];
	particles[calcu_particle].coorZ.val[0] += dt * particles[calcu_particle].velZ.val[0];

	// * Box constrain
	if (BOX_CONTAIN == 1)
	{
		if (particles[calcu_particle].coorX.val[0] < BOX_X_MIN) particles[calcu_particle].coorX.val[0] = BOX_X_MIN;
		if (particles[calcu_particle].coorX.val[0] > BOX_X_MAX) particles[calcu_particle].coorX.val[0] = BOX_X_MAX;

		if (particles[calcu_particle].coorY.val[0] < BOX_Y_MIN) particles[calcu_particle].coorY.val[0] = BOX_Y_MIN;
		if (particles[calcu_particle].coorY.val[0] > BOX_Y_MAX) particles[calcu_particle].coorY.val[0] = BOX_Y_MAX;

		if (particles[calcu_particle].coorZ.val[0] < BOX_Z_MIN) particles[calcu_particle].coorZ.val[0] = BOX_Z_MIN;
		if (particles[calcu_particle].coorZ.val[0] > BOX_Z_MAX) particles[calcu_particle].coorZ.val[0] = BOX_Z_MAX;

		if (CASE_DIM == 0) particles[calcu_particle].coorZ.val[0] = 0.0;
	}

	particles[calcu_particle].time.val[2] = particles[calcu_particle].time.val[1];
	particles[calcu_particle].time.val[1] = particles[calcu_particle].time.val[0];
	particles[calcu_particle].time.val[0] += dt;
}

void position_time_update_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * _neigh_list_cuda, PARTICLE * _particles_cuda, double dt, int blcks = 50)
{
	int offset;

	int _blcks = blcks; // how much GPUs launch at once
	int _thrds = 1;

	int cycle = (neigh_list_length / blcks) + 2;
	dim3 grid(_blcks, 1, 1);
	dim3 thrd(_thrds, 1, 1);

	offset = 0;
	for (int i = 0; i < cycle; ++i)
	{
		_postion_time_update_cuda << <grid, thrd >> >(neigh_list_length, neigh_list_width, _neigh_list_cuda, _particles_cuda, dt, offset);
		offset += blcks;
	}
}