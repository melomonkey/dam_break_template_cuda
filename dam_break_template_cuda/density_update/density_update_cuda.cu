#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\basic_func\basic_func.h"
#include <vector>
using namespace std;

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__global__ void _density_update_schm1_cuda_part1(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * neighbors_list, PARTICLE * particles, unsigned int offset)
{
	unsigned int pos = blockIdx.x + offset;

	if (pos > neigh_list_length - 1) return;

	unsigned int calcu_particle = neighbors_list[0 + neigh_list_width * pos] - 1;

    double local_x, local_y, local_z;
	local_x = particles[calcu_particle].coorX.val[0];
	local_y = particles[calcu_particle].coorY.val[0];
	local_z = particles[calcu_particle].coorZ.val[0];
	double h = particles[calcu_particle].smthR.val[0];

	double tmp_dev = 0.0;
	for (int j = 1; j < neigh_list_width; ++j)
	{
		if (neighbors_list[j + pos * neigh_list_width] == 0) break;

		unsigned int label_ij = neighbors_list[j + pos * neigh_list_width] - 1;

		if (particles[label_ij].id == 0 || particles[label_ij].id == 1 || particles[label_ij].id == 2 || particles[label_ij].id == 3) continue;

		double tmp_x, tmp_y, tmp_z;
		tmp_x = particles[label_ij].coorX.val[0];
		tmp_y = particles[label_ij].coorY.val[0];
		tmp_z = particles[label_ij].coorZ.val[0];

		double dis_x, dis_y, dis_z;
		dis_x = local_x - tmp_x;
		dis_y = local_y - tmp_y;
		dis_z = local_z - tmp_z;

		double coeff0 = particles[calcu_particle].density.val[0] * particles[label_ij].mass.val[0] / particles[label_ij].density.val[0];

		tmp_dev += coeff0 * (particles[calcu_particle].velX.val[0] - particles[label_ij].velX.val[0]) * kernel_function_1dev_gpu(dis_x, h);
		tmp_dev += coeff0 * (particles[calcu_particle].velY.val[0] - particles[label_ij].velY.val[0]) * kernel_function_1dev_gpu(dis_y, h);
		if (CASE_DIM == 1) tmp_dev += coeff0 * (particles[calcu_particle].velZ.val[0] - particles[label_ij].velZ.val[0]) * kernel_function_1dev_gpu(dis_z, h);
	}

	double tmp = tmp_dev;

	particles[calcu_particle].density.val[1] = tmp;
}

__global__ void _density_update_schm1_cuda_part2(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * neighbors_list, PARTICLE * particles, double dt, unsigned int offset)
{
	unsigned int pos = blockIdx.x + offset;

	if (pos > neigh_list_length - 1) return;

	unsigned int calcu_particle = neighbors_list[0 + neigh_list_width * pos] - 1;

	double change = particles[calcu_particle].density.val[1] * dt;
	change = constrain_d_gpu(change, -DENS_CHANGE_LIM, DENS_CHANGE_LIM);

	particles[calcu_particle].density.val[0] += change;
	particles[calcu_particle].density.val[0] = constrain_d_gpu(particles[calcu_particle].density.val[0], DENS_LIM_L, DENS_LIM_U);
}

void density_update_schm1_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * _neigh_list_cuda, PARTICLE * _particles_cuda, double dt, int blcks = 50)
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
		_density_update_schm1_cuda_part1 << <grid, thrd >> >(neigh_list_length, neigh_list_width, _neigh_list_cuda, _particles_cuda, offset);
		offset += blcks;
	}

	offset = 0;
	for (int i = 0; i < cycle; ++i)
	{
		_density_update_schm1_cuda_part2 << <grid, thrd >> >(neigh_list_length, neigh_list_width, _neigh_list_cuda, _particles_cuda, dt, offset);
		offset += blcks;
	}
}