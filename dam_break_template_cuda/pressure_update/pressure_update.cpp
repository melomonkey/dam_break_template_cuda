#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\basic_func\basic_func.h"
#include <vector>
using namespace std;

void pressure_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles)
{
	for(int i = 0; i < neighbors_list.size(); ++ i)
	{	
		int calcu_particle = neighbors_list[i][0];

		// Pressure update
		particles[calcu_particle].pressure.val[2] = particles[neighbors_list[i][0]].pressure.val[1]; 
	    particles[calcu_particle].pressure.val[1] = particles[neighbors_list[i][0]].pressure.val[0];
 
		double density_ratio = particles[calcu_particle].density.val[0] / REF_DENSITY; 
		particles[calcu_particle].pressure.val[0] = COMPRS_LIM_CONST_B * (pow((density_ratio), 7.0) - 1.0);
 
		particles[calcu_particle].pressure.val[0] = constrain_d(particles[calcu_particle].pressure.val[0], PRESSURE_LIM_L, PRESSURE_LIM_U);  
	}	
}

void pressure_update_schm1_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, int threads_num = 5)
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

			// Pressure update
			particles[calcu_particle].pressure.val[2] = particles[neighbors_list[i][0]].pressure.val[1];
			particles[calcu_particle].pressure.val[1] = particles[neighbors_list[i][0]].pressure.val[0];

			double density_ratio = particles[calcu_particle].density.val[0] / REF_DENSITY;
			particles[calcu_particle].pressure.val[0] = COMPRS_LIM_CONST_B * (pow((density_ratio), 7.0) - 1.0);

			particles[calcu_particle].pressure.val[0] = constrain_d(particles[calcu_particle].pressure.val[0], PRESSURE_LIM_L, PRESSURE_LIM_U);
		}
	}
}

void pressure_update_schm2(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles)
{

	for(int i = 0; i < neighbors_list.size(); ++ i)
	{	
		int calcu_particle = neighbors_list[i][0];

		// Pressure update
		particles[calcu_particle].pressure.val[2] = particles[neighbors_list[i][0]].pressure.val[1]; 
	    particles[calcu_particle].pressure.val[1] = particles[neighbors_list[i][0]].pressure.val[0];

		particles[calcu_particle].pressure.val[0] = COMPRS_LIM_CONST_K * (particles[calcu_particle].density.val[0] - REF_DENSITY);
 
		particles[calcu_particle].pressure.val[0] = constrain_d(particles[calcu_particle].pressure.val[0], PRESSURE_LIM_L, PRESSURE_LIM_U);  
	}	
}