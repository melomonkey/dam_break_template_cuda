#include "struct_particle.h"
#include "configuration.h"
#include "basic_func.h"
#include "global_variables.h"
#include <vector>
using namespace std;

void energy_pressure_term_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles, double dt)
{
	for(int i = 0; i < neighbors_list.size(); ++ i)
	{
		int calcu_particle = neighbors_list[i][0];

		double local_x, local_y, local_z;
		local_x = (*particles)[calcu_particle].coorX.val[0]; 
		local_y = (*particles)[calcu_particle].coorY.val[0];
		local_z = (*particles)[calcu_particle].coorZ.val[0];
		double h = (*particles)[calcu_particle].smthR.val[0];

		double energy_deriv = 0.0; 

		for(int j = 1; j < neighbors_list[i].size(); ++ j)
		{
			int label_ij = neighbors_list[i][j];

			double tmp_x, tmp_y, tmp_z; 
			tmp_x = (*particles)[label_ij].coorX.val[0]; 
			tmp_y = (*particles)[label_ij].coorY.val[0];
			tmp_z = (*particles)[label_ij].coorZ.val[0];

			double dis_x, dis_y, dis_z;  
			dis_x = local_x - tmp_x; 
			dis_y = local_y - tmp_y; 
			dis_z = local_z - tmp_z; 

			// Interactive particles
			if((*particles)[label_ij].id >= 4)
			{
				double coeff0; 
				coeff0 = -((*particles)[label_ij].mass.val[0] * ((*particles)[calcu_particle].pressure.val[0] + REF_PRESSURE + (*particles)[label_ij].pressure.val[0] + REF_PRESSURE) / ((*particles)[calcu_particle].density.val[0] * (*particles)[label_ij].density.val[0])); 

				double de_kernel_val[3] = { kernel_function_1dev(dis_x, h), kernel_function_1dev(dis_y, h), kernel_function_1dev(dis_z, h) };

				energy_deriv += .5 * coeff0 * ((*particles)[calcu_particle].velX.val[0] - (*particles)[label_ij].velX.val[0]) * de_kernel_val[0]; 
				energy_deriv += .5 * coeff0 * ((*particles)[calcu_particle].velY.val[0] - (*particles)[label_ij].velY.val[0]) * de_kernel_val[1]; 
				if(CASE_DIM == 1) energy_deriv += .5 * coeff0 * ((*particles)[calcu_particle].velZ.val[0] - (*particles)[label_ij].velZ.val[0]) * de_kernel_val[2]; 
			}
		}
			

		/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

		(*particles)[calcu_particle].energy.val[2] = (*particles)[calcu_particle].energy.val[1]; 
		(*particles)[calcu_particle].energy.val[1] = (*particles)[calcu_particle].energy.val[0];
		(*particles)[calcu_particle].energy.val[0] = energy_deriv * dt;
	}

}

void energy_pressure_term_update_schm2(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles, double dt)
{
	for(int i = 0; i < neighbors_list.size(); ++ i)
	{
		int calcu_particle = neighbors_list[i][0];

		double local_x, local_y, local_z;
		local_x = (*particles)[calcu_particle].coorX.val[0]; 
		local_y = (*particles)[calcu_particle].coorY.val[0];
		local_z = (*particles)[calcu_particle].coorZ.val[0];
		double h = (*particles)[calcu_particle].smthR.val[0];

		double energy_deriv = 0.0; 

		for(int j = 1; j < neighbors_list[i].size(); ++ j)
		{
			int label_ij = neighbors_list[i][j];

			double tmp_x, tmp_y, tmp_z; 
			tmp_x = (*particles)[label_ij].coorX.val[0]; 
			tmp_y = (*particles)[label_ij].coorY.val[0];
			tmp_z = (*particles)[label_ij].coorZ.val[0];

			double dis_x, dis_y, dis_z;  
			dis_x = local_x - tmp_x; 
			dis_y = local_y - tmp_y; 
			dis_z = local_z - tmp_z; 

			// Interactive particles
			if((*particles)[label_ij].id >= 4)
			{
				double coeff0; 
				coeff0 = (((*particles)[calcu_particle].pressure.val[0] + REF_PRESSURE) / pow((*particles)[calcu_particle].density.val[0], 2.0)) + (((*particles)[label_ij].pressure.val[0] + REF_PRESSURE) / pow((*particles)[label_ij].density.val[0], 2.0)); 
				coeff0 = -(*particles)[label_ij].mass.val[0] * coeff0;

				double de_kernel_val[3] = { kernel_function_1dev(dis_x, h), kernel_function_1dev(dis_y, h), kernel_function_1dev(dis_z, h) };

				energy_deriv += .5 * coeff0 * ((*particles)[calcu_particle].velX.val[0] - (*particles)[label_ij].velX.val[0]) * de_kernel_val[0]; 
				energy_deriv += .5 * coeff0 * ((*particles)[calcu_particle].velY.val[0] - (*particles)[label_ij].velY.val[0]) * de_kernel_val[1]; 
				if(CASE_DIM == 1) energy_deriv += .5 * coeff0 * ((*particles)[calcu_particle].velZ.val[0] - (*particles)[label_ij].velZ.val[0]) * de_kernel_val[2]; 
			}
		}
			

		/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

		(*particles)[calcu_particle].energy.val[2] = (*particles)[calcu_particle].energy.val[1]; 
		(*particles)[calcu_particle].energy.val[1] = (*particles)[calcu_particle].energy.val[0];
		(*particles)[calcu_particle].energy.val[0] = energy_deriv * dt;
	}

}