#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\basic_func\basic_func.h"
#include <vector>
using namespace std;

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__global__ void _pressure_update_schm1(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * neighbors_list, PARTICLE * particles, unsigned int offset)
{
	unsigned int pos = blockIdx.x + offset; 

	if (pos > neigh_list_length - 1) return; 

	unsigned int calcu_particle = neighbors_list[0 + neigh_list_width * pos] - 1;

	double density_ratio = particles[calcu_particle].density.val[0] / REF_DENSITY;
	particles[calcu_particle].pressure.val[0] = COMPRS_LIM_CONST_B * (pow((density_ratio), 7.0) - 1.0);

	particles[calcu_particle].pressure.val[0] = constrain_d_gpu(particles[calcu_particle].pressure.val[0], PRESSURE_LIM_L, PRESSURE_LIM_U);
}

void pressure_update_schm1_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * _neigh_list_cuda, PARTICLE * _particles_cuda, int blcks = 50)
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
		_pressure_update_schm1 << <grid, thrd >> >(neigh_list_length, neigh_list_width, _neigh_list_cuda, _particles_cuda, offset);
		offset += blcks;
	}
}