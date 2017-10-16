#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\basic_func\basic_func.h"
#include <vector>
using namespace std;

void acc_viscosity_term_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles)
{
	for(int i = 0; i < neighbors_list.size(); ++ i)
	{
		int calcu_particle = neighbors_list[i][0];

		double local_x, local_y, local_z;
		local_x = (*particles)[calcu_particle].coorX.val[0]; 
		local_y = (*particles)[calcu_particle].coorY.val[0];
		local_z = (*particles)[calcu_particle].coorZ.val[0];
		double h = (*particles)[calcu_particle].smthR.val[0];

		double acc_x = 0.0, 
			   acc_y = 0.0, 
			   acc_z = 0.0; 

		for(int j = 1; j < neighbors_list[i].size(); ++ j)
		{
			int label_ij = neighbors_list[i][j];

			double tmp_x, tmp_y, tmp_z; 
			tmp_x = (*particles)[label_ij].coorX.val[0]; 
			tmp_y = (*particles)[label_ij].coorY.val[0];
			tmp_z = (*particles)[label_ij].coorZ.val[0];

			double dis_x, dis_y, dis_z; 
			double dis_sqr, dis; 
			dis_x = local_x - tmp_x; 
			dis_y = local_y - tmp_y; 
			dis_z = local_z - tmp_z; 
			dis_sqr = pow(dis_x, 2.0) + pow(dis_y, 2.0) + pow(dis_z, 2.0); 
			dis = pow(dis_sqr, .5); 

			// Interactive particles
			if((*particles)[label_ij].id >= 4)
			{
				double coeff1; 			
				coeff1 = (*particles)[label_ij].mass.val[0] * .5 * ((*particles)[calcu_particle].viscosity.val[0] + (*particles)[label_ij].viscosity.val[0]); 
				coeff1 /= ((*particles)[calcu_particle].density.val[0] * (*particles)[label_ij].density.val[0]);

				double de_kernael_val_dis;
				de_kernael_val_dis = kernel_function_1dev(dis, h);

				if(dis == 0.0) dis = 1e70; 

				double de_kernel_val;
		
				acc_x += coeff1 * (de_kernael_val_dis * (1.0 / dis)) * ((*particles)[calcu_particle].velX.val[0] - (*particles)[label_ij].velX.val[0]); 

				acc_y += coeff1 * (de_kernael_val_dis * (1.0 / dis)) * ((*particles)[calcu_particle].velY.val[0] - (*particles)[label_ij].velY.val[0]);

				acc_z += coeff1 * (de_kernael_val_dis * (1.0 / dis)) * ((*particles)[calcu_particle].velZ.val[0] - (*particles)[label_ij].velZ.val[0]);

			}

		}
	
		/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

		{ 
			acc_x = constrain_d(acc_x, ACC_X_LIM_L, ACC_X_LIM_U); 
			acc_y = constrain_d(acc_y, ACC_Y_LIM_L, ACC_Y_LIM_U); 
			acc_z = constrain_d(acc_z, ACC_Z_LIM_L, ACC_Z_LIM_U); 			
		}

		/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

		(*particles)[calcu_particle].accX.val[2] = (*particles)[calcu_particle].accX.val[1]; 
		(*particles)[calcu_particle].accX.val[1] = (*particles)[calcu_particle].accX.val[0];
		(*particles)[calcu_particle].accX.val[0] += acc_x;

		(*particles)[calcu_particle].accY.val[2] = (*particles)[calcu_particle].accY.val[1]; 
		(*particles)[calcu_particle].accY.val[1] = (*particles)[calcu_particle].accY.val[0];
		(*particles)[calcu_particle].accY.val[0] += acc_y; 

		(*particles)[calcu_particle].accZ.val[2] = (*particles)[calcu_particle].accZ.val[1]; 
		(*particles)[calcu_particle].accZ.val[1] = (*particles)[calcu_particle].accZ.val[0];
		(*particles)[calcu_particle].accZ.val[0] += acc_z; 
	}
}