#include "struct_particle.h"
#include "configuration.h"
#include "basic_func.h"
#include "global_variables.h"
#include <vector>
using namespace std;

void temperature_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles, double dt)
{ 
	for(int i = 0; i < neighbors_list.size(); ++ i)
	{
		int calcu_particle = neighbors_list[i][0]; 

		double local_x, local_y, local_z;
		local_x = (*particles)[calcu_particle].coorX.val[0]; 
		local_y = (*particles)[calcu_particle].coorY.val[0];
		local_z = (*particles)[calcu_particle].coorZ.val[0];
		double h = (*particles)[calcu_particle].smthR.val[0];

		double k = 0.6; 
		double c_i = 4.2e3; 
		double Q = 0.1; 
		
		double temp_deriv[3] = { 0.0, 0.0, 0.0 }; 

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

			double coeff0 = k * ((*particles)[label_ij].mass.val[0] / (*particles)[label_ij].density.val[0]); 
			coeff0 /= (*particles)[calcu_particle].density.val[0] * c_i; 

			double diff_temp = (*particles)[label_ij].temperature.val[0] - (*particles)[calcu_particle].temperature.val[0];
			
			temp_deriv[0] += coeff0 * diff_temp * kernel_function(dis_x, h) / h / h; 
			temp_deriv[1] += coeff0 * diff_temp * kernel_function(dis_y, h) / h / h;
			temp_deriv[2] += coeff0 * diff_temp * kernel_function(dis_z, h) / h / h;
		}

		double sum = temp_deriv[0] + temp_deriv[1]; 
		if(CASE_DIM == 1) sum += temp_deriv[2]; 

		(*particles)[calcu_particle].temperature.val[2] = (*particles)[calcu_particle].temperature.val[1];
		(*particles)[calcu_particle].temperature.val[1] = (*particles)[calcu_particle].temperature.val[0];
		(*particles)[calcu_particle].temperature.val[0] += sum * dt; 
	}
}

void temperature_update_schm2(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles, double dt)
{ 
	vector<double> temperature_dev_x((*particles).size(), 0.0); 
	vector<double> temperature_dev_y((*particles).size(), 0.0); 
	vector<double> temperature_dev_z((*particles).size(), 0.0); 

	for(int i = 0; i < neighbors_list.size(); ++ i)
	{
		int calcu_particle = neighbors_list[i][0]; 

		double local_x, local_y, local_z;
		local_x = (*particles)[calcu_particle].coorX.val[0]; 
		local_y = (*particles)[calcu_particle].coorY.val[0];
		local_z = (*particles)[calcu_particle].coorZ.val[0];
		double h = (*particles)[calcu_particle].smthR.val[0];

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

			double coeff0 = ((*particles)[label_ij].mass.val[0] / (*particles)[label_ij].density.val[0]); 

			double diff_temp = (*particles)[label_ij].temperature.val[0] - (*particles)[calcu_particle].temperature.val[0];

			temperature_dev_x[calcu_particle] += coeff0 * diff_temp * kernel_function_1dev(dis_x, h); 
			temperature_dev_y[calcu_particle] += coeff0 * diff_temp * kernel_function_1dev(dis_y, h); 
			if(CASE_DIM == 1) temperature_dev_z[calcu_particle] += coeff0 * diff_temp * kernel_function_1dev(dis_z, h); 
		}
	}

	for(int i = 0; i < neighbors_list.size(); ++ i)
	{
		int calcu_particle = neighbors_list[i][0]; 

		double local_x, local_y, local_z;
		local_x = (*particles)[calcu_particle].coorX.val[0]; 
		local_y = (*particles)[calcu_particle].coorY.val[0];
		local_z = (*particles)[calcu_particle].coorZ.val[0];
		double h = (*particles)[calcu_particle].smthR.val[0];

		double k = 0.6; 
		double c_i = 4.2e3; 
		double Q = 0.1; 
		
		double temp_deriv[3] = { 0.0, 0.0, 0.0 }; 

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

			double coeff0 = k * ((*particles)[label_ij].mass.val[0] / (*particles)[label_ij].density.val[0]); 
			coeff0 /= (*particles)[calcu_particle].density.val[0] * c_i; 

			double diff_temp_dev_x = temperature_dev_x[label_ij] - temperature_dev_x[calcu_particle], 
				   diff_temp_dev_y = temperature_dev_y[label_ij] - temperature_dev_y[calcu_particle], 
				   diff_temp_dev_z = temperature_dev_z[label_ij] - temperature_dev_z[calcu_particle];
			
			temp_deriv[0] += coeff0 * diff_temp_dev_x * kernel_function_1dev(dis_x, h); 
			temp_deriv[1] += coeff0 * diff_temp_dev_y * kernel_function_1dev(dis_y, h);
			temp_deriv[2] += coeff0 * diff_temp_dev_z * kernel_function_1dev(dis_z, h);
		}

		double sum = temp_deriv[0] + temp_deriv[1]; 
	    if(CASE_DIM == 1) sum += temp_deriv[2]; 

	    (*particles)[calcu_particle].temperature.val[2] = (*particles)[calcu_particle].temperature.val[1];
	    (*particles)[calcu_particle].temperature.val[1] = (*particles)[calcu_particle].temperature.val[0];
	    (*particles)[calcu_particle].temperature.val[0] += sum * dt;
	}

}

void temperature_update_schm11(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles, double dt)
{ 
	for(int i = 0; i < neighbors_list.size(); ++ i)
	{
		int calcu_particle = neighbors_list[i][0]; 

		double local_x, local_y, local_z;
		local_x = (*particles)[calcu_particle].coorX.val[0]; 
		local_y = (*particles)[calcu_particle].coorY.val[0];
		local_z = (*particles)[calcu_particle].coorZ.val[0];
		double h = (*particles)[calcu_particle].smthR.val[0];

		double k = 1.0; 
		double c_i = 1.0; 
		double Q = 0.1; 

		double temp_deriv_x = 0.0;
		double temp_deriv_y = 0.0; 
		double tmp_sum_x = 0.0;
		double tmp_sum_y = 0.0; 
		
		double temp_deriv[3] = { 0.0, 0.0, 0.0 }; 
		double tmp_sum[3] = { 0.0, 0.0, 0.0 }; 

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

			double coeff0 = ((*particles)[label_ij].mass.val[0] / (*particles)[label_ij].density.val[0]); 
			coeff0 /= (*particles)[calcu_particle].density.val[0]; 
			
			temp_deriv_x += 1.0 * ((*particles)[label_ij].temperature.val[0] - (*particles)[calcu_particle].temperature.val[0]) * kernel_function(dis_x, h); 
			tmp_sum_x += kernel_function(dis_x, h); 
			temp_deriv_y += 1.0 * ((*particles)[label_ij].temperature.val[0] - (*particles)[calcu_particle].temperature.val[0]) * kernel_function(dis_y, h);
			tmp_sum_y += kernel_function(dis_y, h); 
			//temp_deriv += coeff0 * ((*particles)[calcu_particle].temperature_dev.val[0] - (*particles)[label_ij].temperature_dev.val[0]) * kernel_function_1dev(dis_z, h);
		}

		double sum = (temp_deriv_x / tmp_sum_x) + (temp_deriv_y / tmp_sum_y); 

		(*particles)[calcu_particle].temperature.val[2] = (*particles)[calcu_particle].temperature.val[1];
		(*particles)[calcu_particle].temperature.val[1] = (*particles)[calcu_particle].temperature.val[0];
		(*particles)[calcu_particle].temperature.val[0] += sum * dt; 
	}
}