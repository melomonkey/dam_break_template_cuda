#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\basic_func\basic_func.h"
#include <vector>
using namespace std;

void density_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, double dt)
{
	for (int i = 0; i < neighbors_list.size(); ++i)
	{
		int calcu_particle = neighbors_list[i][0];

		double local_x, local_y, local_z;
		local_x = particles[calcu_particle].coorX.val[0];
		local_y = particles[calcu_particle].coorY.val[0];
		local_z = particles[calcu_particle].coorZ.val[0];
		double h = particles[calcu_particle].smthR.val[0];

		double tmp_dev = 0.0;
		for (int j = 1; j < neighbors_list[i].size(); ++j)
		{
			int label_ij = neighbors_list[i][j];

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

			tmp_dev += coeff0 * -(particles[label_ij].velX.val[0] - particles[calcu_particle].velX.val[0]) * kernel_function_1dev(dis_x, h);
			tmp_dev += coeff0 * -(particles[label_ij].velY.val[0] - particles[calcu_particle].velY.val[0]) * kernel_function_1dev(dis_y, h);
			if (CASE_DIM == 1) tmp_dev += coeff0 * -(particles[calcu_particle].velZ.val[0] - particles[label_ij].velZ.val[0]) * kernel_function_1dev(dis_z, h);
		}

		double tmp = tmp_dev;

		particles[calcu_particle].density.val[1] = tmp;
	}

	for (int i = 0; i < neighbors_list.size(); ++i)
	{
		int calcu_particle = neighbors_list[i][0];

		double change = particles[calcu_particle].density.val[1] * dt;
		change = constrain_d(change, -DENS_CHANGE_LIM, DENS_CHANGE_LIM);

		particles[calcu_particle].density.val[0] += change;
		particles[calcu_particle].density.val[0] = constrain_d(particles[calcu_particle].density.val[0], DENS_LIM_L, DENS_LIM_U);
	}
}

void density_update_schm1_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, double dt, int threads_num = 5)
{
	int each_thread_cnt = neighbors_list.size() / threads_num + 1;

#pragma omp parallel for
	for (int ii = 0; ii < threads_num; ++ii)
	{
		int beg_idx, end_idx;
		beg_idx = ii * each_thread_cnt;
		end_idx = (ii + 1) * each_thread_cnt;

		if (beg_idx < 0) beg_idx = 0;
		if (end_idx > neighbors_list.size() - 1) end_idx = neighbors_list.size() - 1;

		for (int i = beg_idx; i <= end_idx; ++i)
		{
			int calcu_particle = neighbors_list[i][0];

			double local_x, local_y, local_z;
			local_x = particles[calcu_particle].coorX.val[0];
			local_y = particles[calcu_particle].coorY.val[0];
			local_z = particles[calcu_particle].coorZ.val[0];
			double h = particles[calcu_particle].smthR.val[0];

			double tmp_dev = 0.0;
			for (int j = 1; j < neighbors_list[i].size(); ++j)
			{
				int label_ij = neighbors_list[i][j];

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

				tmp_dev += coeff0 * -(particles[label_ij].velX.val[0] - particles[calcu_particle].velX.val[0]) * kernel_function_1dev(dis_x, h);
				tmp_dev += coeff0 * -(particles[label_ij].velY.val[0] - particles[calcu_particle].velY.val[0]) * kernel_function_1dev(dis_y, h);
				if (CASE_DIM == 1) tmp_dev += coeff0 * -(particles[label_ij].velZ.val[0] - particles[calcu_particle].velZ.val[0]) * kernel_function_1dev(dis_z, h);
			}

			double tmp = tmp_dev;

			particles[calcu_particle].density.val[1] = tmp;
		}
	}

#pragma omp parallel for
	for (int ii = 0; ii < threads_num; ++ii)
	{
		int beg_idx, end_idx;
		beg_idx = ii * each_thread_cnt;
		end_idx = (ii + 1) * each_thread_cnt;

		if (beg_idx < 0) beg_idx = 0;
		if (end_idx > neighbors_list.size() - 1) end_idx = neighbors_list.size() - 1;

		for (int i = beg_idx; i <= end_idx; ++i)
		{
			int calcu_particle = neighbors_list[i][0];

			double change = particles[calcu_particle].density.val[1] * dt;
			change = constrain_d(change, -DENS_CHANGE_LIM, DENS_CHANGE_LIM);

			particles[calcu_particle].density.val[0] += change;
			particles[calcu_particle].density.val[0] = constrain_d(particles[calcu_particle].density.val[0], DENS_LIM_L, DENS_LIM_U);
		}
	}
}

void density_update_schm2(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles)
{
	for (int i = 0; i < neighbors_list.size(); ++i)
	{
		int calcu_particle = neighbors_list[i][0];

		double local_x, local_y, local_z;
		local_x = particles[calcu_particle].coorX.val[0];
		local_y = particles[calcu_particle].coorY.val[0];
		local_z = particles[calcu_particle].coorZ.val[0];
		double h = particles[calcu_particle].smthR.val[0];

		double tmp = 0.0;
		for (int j = 1; j < neighbors_list[i].size(); ++j)
		{
			int label_ij = neighbors_list[i][j];

			if (particles[label_ij].id == 0 || particles[label_ij].id == 1 || particles[label_ij].id == 2 || particles[label_ij].id == 3) continue;

			double tmp_x, tmp_y, tmp_z;
			tmp_x = particles[label_ij].coorX.val[0];
			tmp_y = particles[label_ij].coorY.val[0];
			tmp_z = particles[label_ij].coorZ.val[0];

			double dis_x, dis_y, dis_z;
			dis_x = local_x - tmp_x;
			dis_y = local_y - tmp_y;
			dis_z = local_z - tmp_z;

			double coeff0 = particles[label_ij].mass.val[0];

			tmp += coeff0 * kernel_function(dis_x, h);
			tmp += coeff0 * kernel_function(dis_y, h);
			if (CASE_DIM == 1) tmp += coeff0 * kernel_function(dis_z, h);
		}
		particles[calcu_particle].density.val[0] = tmp;

		particles[calcu_particle].density.val[0] = constrain_d(particles[calcu_particle].density.val[0], DENS_LIM_L, DENS_LIM_U);
	}

}

void density_update_schm2_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, int threads_num = 5)
{
	int each_thread_cnt = neighbors_list.size() / threads_num + 1;

#pragma omp parallel for
	for (int ii = 0; ii < threads_num; ++ii)
	{
		int beg_idx, end_idx;
		beg_idx = ii * each_thread_cnt;
		end_idx = (ii + 1) * each_thread_cnt;

		if (beg_idx < 0) beg_idx = 0;
		if (end_idx > neighbors_list.size() - 1) end_idx = neighbors_list.size() - 1;

		for (int i = beg_idx; i <= end_idx; ++i)
		{
			int calcu_particle = neighbors_list[i][0];

			double local_x, local_y, local_z;
			local_x = particles[calcu_particle].coorX.val[0];
			local_y = particles[calcu_particle].coorY.val[0];
			local_z = particles[calcu_particle].coorZ.val[0];
			double h = particles[calcu_particle].smthR.val[0];

			double tmp = 0.0;
			for (int j = 1; j < neighbors_list[i].size(); ++j)
			{
				int label_ij = neighbors_list[i][j];

				if (particles[label_ij].id == 0 || particles[label_ij].id == 1 || particles[label_ij].id == 2 || particles[label_ij].id == 3) continue;

				double tmp_x, tmp_y, tmp_z;
				tmp_x = particles[label_ij].coorX.val[0];
				tmp_y = particles[label_ij].coorY.val[0];
				tmp_z = particles[label_ij].coorZ.val[0];

				double dis_x, dis_y, dis_z;
				dis_x = local_x - tmp_x;
				dis_y = local_y - tmp_y;
				dis_z = local_z - tmp_z;

				double coeff0 = particles[label_ij].mass.val[0];

				tmp += coeff0 * kernel_function(dis_x, h);
				tmp += coeff0 * kernel_function(dis_y, h);
				if (CASE_DIM == 1) tmp += coeff0 * kernel_function(dis_z, h);
			}
			particles[calcu_particle].density.val[0] = tmp;

			particles[calcu_particle].density.val[0] = constrain_d(particles[calcu_particle].density.val[0], DENS_LIM_L, DENS_LIM_U);
		}
	}
}