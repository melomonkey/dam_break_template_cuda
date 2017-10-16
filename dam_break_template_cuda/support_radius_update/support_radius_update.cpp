#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\basic_func\basic_func.h"
#include <vector>
using namespace std;

void support_radius_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles)
{
	double dimension;
	if (CASE_DIM == 0) dimension = 2.0; 
	if (CASE_DIM == 1) dimension = 3.0; 

	for (int i = 0; i < neighbors_list.size(); ++i)
	{
		int calcu_particle = neighbors_list[i][0];

		double ref_density = particles[calcu_particle].density.val[2];
		double ref_smthR = particles[calcu_particle].smthR.val[2];

		particles[calcu_particle].smthR.val[0] = ref_smthR * pow((ref_density / particles[calcu_particle].density.val[0]), 1.0 / dimension);
	}
}

void support_radius_update_schm2(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, double dt)
{
	double dimension;
	if (CASE_DIM == 0) dimension = 2.0; 
	if (CASE_DIM == 1) dimension = 3.0; 

	for (int i = 0; i < neighbors_list.size(); ++i)
	{
		int calcu_particle = neighbors_list[i][0];

		double density_dev;
		density_dev = particles[calcu_particle].density.val[1];

		double tmp_h = -(1.0 / dimension) * (particles[calcu_particle].smthR.val[0] / particles[calcu_particle].density.val[0]) * density_dev;

		particles[calcu_particle].smthR.val[2] = particles[calcu_particle].smthR.val[1];
		particles[calcu_particle].smthR.val[1] = particles[calcu_particle].smthR.val[0];
		particles[calcu_particle].smthR.val[0] += tmp_h * dt;

	}
}

void support_radius_update_schm2_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, double dt, int threads_num = 5)
{
	double dimension;
	if (CASE_DIM == 0) dimension = 2.0;
	if (CASE_DIM == 1) dimension = 3.0;

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

			double density_dev;
			density_dev = particles[calcu_particle].density.val[1];

			double tmp_h = -(1.0 / dimension) * (particles[calcu_particle].smthR.val[0] / particles[calcu_particle].density.val[0]) * density_dev;

			particles[calcu_particle].smthR.val[2] = particles[calcu_particle].smthR.val[1];
			particles[calcu_particle].smthR.val[1] = particles[calcu_particle].smthR.val[0];
			particles[calcu_particle].smthR.val[0] += tmp_h * dt;

		}
	}
}