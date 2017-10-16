#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\basic_func\basic_func.h"
#include <vector>
using namespace std;

void velocity_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, double dt)
{
	for (int i = 0; i < neighbors_list.size(); ++i)
	{
		int calcu_particle = neighbors_list[i][0];

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
			particles[calcu_particle].velX.val[0] = constrain_d(particles[calcu_particle].velX.val[0], VEL_X_LIM_L, VEL_X_LIM_U);
			particles[calcu_particle].velY.val[0] = constrain_d(particles[calcu_particle].velY.val[0], VEL_Y_LIM_L, VEL_Y_LIM_U);
			particles[calcu_particle].velZ.val[0] = constrain_d(particles[calcu_particle].velZ.val[0], VEL_Z_LIM_L, VEL_Z_LIM_U);
		}

	}
}

void velocity_update_schm1_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, double dt, int threads_num = 5)
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
				particles[calcu_particle].velX.val[0] = constrain_d(particles[calcu_particle].velX.val[0], VEL_X_LIM_L, VEL_X_LIM_U);
				particles[calcu_particle].velY.val[0] = constrain_d(particles[calcu_particle].velY.val[0], VEL_Y_LIM_L, VEL_Y_LIM_U);
				particles[calcu_particle].velZ.val[0] = constrain_d(particles[calcu_particle].velZ.val[0], VEL_Z_LIM_L, VEL_Z_LIM_U);
			}
		}
	}
}

void velocity_update_schm2(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, double dt)
{
	for (int i = 0; i < neighbors_list.size(); ++i)
	{
		int calcu_particle = neighbors_list[i][0];

		double h = particles[calcu_particle].smthR.val[0];

		// Velocity X update
		particles[calcu_particle].velX.val[2] = particles[calcu_particle].velX.val[1];
		particles[calcu_particle].velX.val[1] = particles[calcu_particle].velX.val[0];
		particles[calcu_particle].velX.val[0] += dt * (.6 * particles[calcu_particle].accX.val[0] + .4 * particles[calcu_particle].accX.val[1]);

		// Velocity Y update
		particles[calcu_particle].velY.val[2] = particles[calcu_particle].velY.val[1];
		particles[calcu_particle].velY.val[1] = particles[calcu_particle].velY.val[0];
		particles[calcu_particle].velY.val[0] += dt * (.6 * particles[calcu_particle].accY.val[0] + .4 * particles[calcu_particle].accY.val[1]);

		// Velocity Z update
		particles[calcu_particle].velZ.val[2] = particles[calcu_particle].velZ.val[1];
		particles[calcu_particle].velZ.val[1] = particles[calcu_particle].velZ.val[0];
		particles[calcu_particle].velZ.val[0] += dt * (.6 * particles[calcu_particle].accZ.val[0] + .4 * particles[calcu_particle].accZ.val[1]);

		if (BOX_CONTAIN == 1)
		{
			/*if(particles[calcu_particle].coorX.val[0] <= BOX_X_MIN) particles[calcu_particle].velX.val[0] = 0.0;
			if(particles[calcu_particle].coorX.val[0] >= BOX_X_MAX) particles[calcu_particle].velX.val[0] = 0.0;

			if(particles[calcu_particle].coorY.val[0] <= BOX_Y_MIN) particles[calcu_particle].velY.val[0] = 0.0;
			if(particles[calcu_particle].coorY.val[0] >= BOX_Y_MAX) particles[calcu_particle].velY.val[0] = 0.0;

			if(particles[calcu_particle].coorZ.val[0] <= BOX_Z_MIN) particles[calcu_particle].velZ.val[0] = 0.0;
			if(particles[calcu_particle].coorZ.val[0] >= BOX_Z_MAX) particles[calcu_particle].velZ.val[0] = 0.0;*/

			if (CASE_DIM == 0) particles[calcu_particle].velZ.val[0] = 0.0;
		}

		{
			particles[calcu_particle].velX.val[0] = constrain_d(particles[calcu_particle].velX.val[0], VEL_X_LIM_L, VEL_X_LIM_U);
			particles[calcu_particle].velY.val[0] = constrain_d(particles[calcu_particle].velY.val[0], VEL_Y_LIM_L, VEL_Y_LIM_U);
			particles[calcu_particle].velZ.val[0] = constrain_d(particles[calcu_particle].velZ.val[0], VEL_Z_LIM_L, VEL_Z_LIM_U);
		}

	}
}