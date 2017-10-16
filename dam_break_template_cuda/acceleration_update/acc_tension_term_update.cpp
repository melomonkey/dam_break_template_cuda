#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\basic_func\basic_func.h"
#include <vector>
using namespace std;

void acc_tension_term_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles)
{
	for(int i = 0; i < neighbors_list.size(); ++ i)
	{
		int calcu_particle = neighbors_list[i][0];

		int neigh_num = neighbors_list[i].size(); 

		double local_x, local_y, local_z;
		local_x = (*particles)[calcu_particle].coorX.val[0]; 
		local_y = (*particles)[calcu_particle].coorY.val[0];
		local_z = (*particles)[calcu_particle].coorZ.val[0];
		double h = (*particles)[calcu_particle].smthR.val[0];

		double surf_dir[3] = { 0.0, 0.0, 0.0 }; 
		double surf_dir_len = 0.0; 

		for(int j = 0; j < neighbors_list[i].size(); ++ j)
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

			double coeff0 = (*particles)[label_ij].mass.val[0] / (*particles)[label_ij].density.val[0]; 

			surf_dir[0] += coeff0 * kernel_function_1dev(dis_x, h); 
			surf_dir[1] += coeff0 * kernel_function_1dev(dis_y, h); 
			surf_dir[2] += coeff0 * kernel_function_1dev(dis_z, h); 
		}

		surf_dir_len = pow(surf_dir[0], 2.0) + pow(surf_dir[1], 2.0) + pow(surf_dir[2], 2.0); 
		surf_dir_len = powl(surf_dir_len, .5); 
		for(int k = 0; k < 3; ++ k) surf_dir[k] /= surf_dir_len;

		// Surface particles with tension
		if(surf_dir_len >= INI_TEN_DIR_LEN * 0.7)
		{
			(*particles)[calcu_particle].accX.val[0] += surf_dir[0] * INI_TENSION_ACC * (surf_dir_len / INI_TEN_DIR_LEN);
		    (*particles)[calcu_particle].accY.val[0] += surf_dir[1] * INI_TENSION_ACC * (surf_dir_len / INI_TEN_DIR_LEN);
		    (*particles)[calcu_particle].accZ.val[0] += surf_dir[2] * INI_TENSION_ACC * (surf_dir_len / INI_TEN_DIR_LEN);
		}
		
	}
}