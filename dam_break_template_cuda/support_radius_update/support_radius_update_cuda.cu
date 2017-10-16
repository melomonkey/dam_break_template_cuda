#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\basic_func\basic_func.h"
#include <vector>
using namespace std;

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__global__ void _support_radius_update_schm1_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * neighbors_list, PARTICLE * particles, double dt, unsigned int offset)
{
	double dimension;
	if (CASE_DIM == 0) dimension = 2.0;
	if (CASE_DIM == 1) dimension = 3.0;

	unsigned int pos = blockIdx.x + offset;

	if (pos > neigh_list_length - 1) return; 

	unsigned int calcu_particle = neighbors_list[0 + neigh_list_width * pos] - 1;

	double density_dev;
	density_dev = particles[calcu_particle].density.val[1];

	double tmp_h = -(1.0 / dimension) * (particles[calcu_particle].smthR.val[0] / particles[calcu_particle].density.val[0]) * density_dev;

	particles[calcu_particle].smthR.val[2] = particles[calcu_particle].smthR.val[1];
	particles[calcu_particle].smthR.val[1] = particles[calcu_particle].smthR.val[0];
	particles[calcu_particle].smthR.val[0] += tmp_h * dt;
}

void support_radius_update_schm1_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * _neigh_list_cuda, PARTICLE * _particles_cuda, double dt, int blcks = 500)
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
		_support_radius_update_schm1_cuda << <grid, thrd >> >(neigh_list_length, neigh_list_width, _neigh_list_cuda, _particles_cuda, dt, offset);
		offset += blcks;
	}
}