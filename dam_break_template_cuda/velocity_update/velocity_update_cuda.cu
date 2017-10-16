#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\basic_func\basic_func.h"
#include <vector>
using namespace std;

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__global__ void _velocity_update_schm1_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * neighbors_list, PARTICLE * particles, double dt, unsigned int offset)
{
	unsigned int pos = blockIdx.x + offset;

	if (pos > neigh_list_length - 1) return; 

	unsigned int calcu_particle = neighbors_list[0 + neigh_list_width * pos] - 1;

	double h = particles[calcu_particle].smthR.val[0];

	// Velocity X update
	particles[calcu_particle].velX.val[2] = particles[calcu_particle].velX.val[1];
	particles[calcu_particle].velX.val[1] = particles[calcu_particle].velX.val[0];
	particles[calcu_particle].velX.val[0] += dt * particles[calcu_particle].accX.val[0];

	// Velocity Y update
	particles[calcu_particle].velY.val[2] = particles[calcu_particle].velY.val[1];
	particles[calcu_particle].velY.val[1] = particles[calcu_particle].velY.val[0];
	particles[calcu_particle].velY.val[0] += dt * particles[calcu_particle].accY.val[0];

	// Velocity Z update
	particles[calcu_particle].velZ.val[2] = particles[calcu_particle].velZ.val[1];
	particles[calcu_particle].velZ.val[1] = particles[calcu_particle].velZ.val[0];
	particles[calcu_particle].velZ.val[0] += dt * particles[calcu_particle].accZ.val[0];

	if (BOX_CONTAIN == 1)
	{
		if (particles[calcu_particle].coorX.val[0] <= BOX_X_MIN) particles[calcu_particle].velX.val[0] = 0.0;
		if (particles[calcu_particle].coorX.val[0] >= BOX_X_MAX) particles[calcu_particle].velX.val[0] = 0.0;

		if (particles[calcu_particle].coorY.val[0] <= BOX_Y_MIN) particles[calcu_particle].velY.val[0] = 0.0;
		if (particles[calcu_particle].coorY.val[0] >= BOX_Y_MAX) particles[calcu_particle].velY.val[0] = 0.0;

		if (particles[calcu_particle].coorZ.val[0] <= BOX_Z_MIN) particles[calcu_particle].velZ.val[0] = 0.0;
		if (particles[calcu_particle].coorZ.val[0] >= BOX_Z_MAX) particles[calcu_particle].velZ.val[0] = 0.0;

		if (CASE_DIM == 0) particles[calcu_particle].velZ.val[0] = 0.0;
	}

	{
		particles[calcu_particle].velX.val[0] = constrain_d_gpu(particles[calcu_particle].velX.val[0], VEL_X_LIM_L, VEL_X_LIM_U);
		particles[calcu_particle].velY.val[0] = constrain_d_gpu(particles[calcu_particle].velY.val[0], VEL_Y_LIM_L, VEL_Y_LIM_U);
		particles[calcu_particle].velZ.val[0] = constrain_d_gpu(particles[calcu_particle].velZ.val[0], VEL_Z_LIM_L, VEL_Z_LIM_U);
	}
}

void velocity_update_schm1_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * _neigh_list_cuda, PARTICLE * _particles_cuda, double dt, int blcks = 50)
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
		_velocity_update_schm1_cuda << <grid, thrd >> >(neigh_list_length, neigh_list_width, _neigh_list_cuda, _particles_cuda, dt, offset);
		offset += blcks;
	}
}