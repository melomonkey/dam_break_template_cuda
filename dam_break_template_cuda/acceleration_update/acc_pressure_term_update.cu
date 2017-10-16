#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\basic_func\basic_func.h"
#include <vector>
using namespace std;

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__global__ void _acc_pressure_term_update_schm1_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * neighbors_list, PARTICLE * particles, unsigned int offset)
{
	unsigned int pos = blockIdx.x + offset;

	if (pos > neigh_list_length - 1) return; 

	unsigned int calcu_particle = neighbors_list[0 + neigh_list_width * pos] - 1;

	double local_x, local_y, local_z;
	local_x = particles[calcu_particle].coorX.val[0];
	local_y = particles[calcu_particle].coorY.val[0];
	local_z = particles[calcu_particle].coorZ.val[0];
	double h = particles[calcu_particle].smthR.val[0];

	double acc_x = 0.0,
		acc_y = 0.0,
		acc_z = 0.0;

	double _acc[3];
	for (int m = 0; m < 3; ++m) _acc[m] = 0.0;

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

		// Interactive particles
		if (particles[label_ij].id >= 4)
		{
			double coeff0;
			coeff0 = -(particles[calcu_particle].pressure.val[0] + particles[label_ij].pressure.val[0]) / (particles[calcu_particle].density.val[0] * particles[label_ij].density.val[0]);
			coeff0 *= particles[label_ij].mass.val[0];
			double de_kernel_val[3] = { kernel_function_1dev_gpu(dis_x, h), kernel_function_1dev_gpu(dis_y, h), kernel_function_1dev_gpu(dis_z, h) };

			_acc[0] += coeff0 * de_kernel_val[0];
			_acc[1] += coeff0 * de_kernel_val[1];
			_acc[2] += coeff0 * de_kernel_val[2];
		}
	}

	{
		acc_x += _acc[0];
		acc_y += _acc[1];
		acc_z += _acc[2];
	}

	// Body force
	{
		double volume_acc[3];

		if (CASE_DIM == 0) { volume_acc[0] = 0.0; volume_acc[1] = -GRAVITY; volume_acc[2] = 0.0; };
		if (CASE_DIM == 1) { volume_acc[0] = 0.0; volume_acc[1] = 0.0; volume_acc[2] = -GRAVITY; };

		acc_x += volume_acc[0];
		acc_y += volume_acc[1];
		acc_z += volume_acc[2];
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	{
		acc_x = constrain_d_gpu(acc_x, ACC_X_LIM_L, ACC_X_LIM_U);
		acc_y = constrain_d_gpu(acc_y, ACC_Y_LIM_L, ACC_Y_LIM_U);
		acc_z = constrain_d_gpu(acc_z, ACC_Z_LIM_L, ACC_Z_LIM_U);
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	particles[calcu_particle].accX.val[2] = particles[calcu_particle].accX.val[1];
	particles[calcu_particle].accX.val[1] = particles[calcu_particle].accX.val[0];
	particles[calcu_particle].accX.val[0] = acc_x;

	particles[calcu_particle].accY.val[2] = particles[calcu_particle].accY.val[1];
	particles[calcu_particle].accY.val[1] = particles[calcu_particle].accY.val[0];
	particles[calcu_particle].accY.val[0] = acc_y;

	particles[calcu_particle].accZ.val[2] = particles[calcu_particle].accZ.val[1];
	particles[calcu_particle].accZ.val[1] = particles[calcu_particle].accZ.val[0];
	particles[calcu_particle].accZ.val[0] = acc_z;
}

void acc_pressure_term_update_schm1_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * _neigh_list_cuda, PARTICLE * _particles_cuda, double dt, int blcks = 50)
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
		_acc_pressure_term_update_schm1_cuda << <grid, thrd >> >(neigh_list_length, neigh_list_width, _neigh_list_cuda, _particles_cuda, offset);
		offset += blcks;
	}
}