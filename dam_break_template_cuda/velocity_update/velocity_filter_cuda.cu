#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\basic_func\basic_func.h"
#include <vector>
using namespace std;

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__global__ void _velocity_filter_schm1_cuda_part1(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * neighbors_list, PARTICLE * particles, unsigned int offset)
{
	unsigned int pos = blockIdx.x + offset;

	if (pos > neigh_list_length - 1) return;

	unsigned int calcu_particle = neighbors_list[0 + neigh_list_width * pos] - 1;

	double h = particles[calcu_particle].smthR.val[0];

	// The revisements of velocities in XYZ direction
	double tmp_velX_rev = 0.0,
		tmp_velY_rev = 0.0,
		tmp_velZ_rev = 0.0;

	double epsilon = 0.2;

	for (int j = 1; j < neigh_list_width; ++j)
	{
		if (neighbors_list[j + pos * neigh_list_width] == 0) break;

		int label_ij = neighbors_list[j + neigh_list_width * pos] - 1;

		double dis_x, dis_y, dis_z;
		dis_x = particles[label_ij].coorX.val[0] - particles[calcu_particle].coorX.val[0];
		dis_y = particles[label_ij].coorY.val[0] - particles[calcu_particle].coorY.val[0];
		dis_z = particles[label_ij].coorZ.val[0] - particles[calcu_particle].coorZ.val[0];

		double coeff0 = epsilon * particles[label_ij].mass.val[0];
		coeff0 /= (.5 * (particles[label_ij].density.val[0] + particles[calcu_particle].density.val[0]));

		double tmp_x, tmp_y, tmp_z;
		tmp_x = coeff0 * (particles[label_ij].velX.val[0] - particles[calcu_particle].velX.val[0]) * kernel_function_gpu(dis_x, h);
		tmp_y = coeff0 * (particles[label_ij].velY.val[0] - particles[calcu_particle].velY.val[0]) * kernel_function_gpu(dis_y, h);
		tmp_z = coeff0 * (particles[label_ij].velZ.val[0] - particles[calcu_particle].velZ.val[0]) * kernel_function_gpu(dis_z, h);

		tmp_velX_rev += tmp_x;
		tmp_velY_rev += tmp_y;
		tmp_velZ_rev += tmp_z;
	}

	particles[calcu_particle].velX.val[2] = particles[calcu_particle].velX.val[0] + tmp_velX_rev;
	particles[calcu_particle].velY.val[2] = particles[calcu_particle].velY.val[0] + tmp_velY_rev;
	particles[calcu_particle].velZ.val[2] = particles[calcu_particle].velZ.val[0] + tmp_velZ_rev;
}

__global__ void _velocity_filter_schm1_cuda_part2(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * neighbors_list, PARTICLE * particles, unsigned int offset)
{
	unsigned int pos = blockIdx.x + offset;

	if (pos > neigh_list_length - 1) return;

	unsigned int calcu_particle = neighbors_list[0 + neigh_list_width * pos] - 1;

	particles[calcu_particle].velX.val[0] = particles[calcu_particle].velX.val[2];
	particles[calcu_particle].velY.val[0] = particles[calcu_particle].velY.val[2];
	particles[calcu_particle].velZ.val[0] = particles[calcu_particle].velZ.val[2];
}

void velocity_filter_schm1_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * _neigh_list_cuda, PARTICLE * _particles_cuda, int blcks = 500)
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
		_velocity_filter_schm1_cuda_part1 << <grid, thrd >> >(neigh_list_length, neigh_list_width, _neigh_list_cuda, _particles_cuda, offset);
		offset += blcks;
	}

	offset = 0;
	for (int i = 0; i < cycle; ++i)
	{
		_velocity_filter_schm1_cuda_part2 << <grid, thrd >> >(neigh_list_length, neigh_list_width, _neigh_list_cuda, _particles_cuda, offset);
		offset += blcks;
	}
}