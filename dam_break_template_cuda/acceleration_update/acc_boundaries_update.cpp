#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\basic_func\basic_func.h"
#include <vector>
using namespace std;

void acc_boundaries_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles)
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

		/*
		*
		* _acc is a matrix stores the accelerations caused in 2 differernt ways.
		*
		* 1. The first line of the matrix represents the accelerations caused by the momentum bounced back from wall.
		* 2. The second line of the matrix represents the accelerations caused by the replusive force generated by the wall. They will be never bigger than the accelerations the particle pocesses.
		*
		*/
	
		double _acc[2][3]; 
		for(int m = 0; m < 2; ++ m) { for(int n = 0; n < 3; ++ n) _acc[m][n] = 0.0; } 

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

			// Wall particles
			if((*particles)[label_ij].id != 0) continue; 
			if((*particles)[label_ij].id == 0)
			{
				double tmp_vel[3]; 
				tmp_vel[0] = (*particles)[calcu_particle].velX.val[0] - (*particles)[label_ij].velX.val[0]; 
				tmp_vel[1] = (*particles)[calcu_particle].velY.val[0] - (*particles)[label_ij].velY.val[0]; 
				tmp_vel[2] = (*particles)[calcu_particle].velZ.val[0] - (*particles)[label_ij].velZ.val[0]; 

				double tmp_acc[3]; 
				tmp_acc[0] = (*particles)[calcu_particle].accX.val[0] - (*particles)[label_ij].accX.val[0]; 
				tmp_acc[1] = (*particles)[calcu_particle].accY.val[0] - (*particles)[label_ij].accY.val[0];
				tmp_acc[2] = (*particles)[calcu_particle].accZ.val[0] - (*particles)[label_ij].accZ.val[0]; 

				double tmp_dir[3];
				tmp_dir[0] = dis_x; 
				tmp_dir[1] = dis_y; 
				tmp_dir[2] = dis_z; 

				double tmp_boun_dir[3]; 
				tmp_boun_dir[0] = (*particles)[label_ij].normal_vec[0]; 
				tmp_boun_dir[1] = (*particles)[label_ij].normal_vec[1]; 
				tmp_boun_dir[2] = (*particles)[label_ij].normal_vec[2]; 

				/*double tmp_boun_dir_length = tmp_boun_dir[0] * tmp_boun_dir[0] + tmp_boun_dir[1] * tmp_boun_dir[1] + tmp_boun_dir[2] * tmp_boun_dir[2];
				tmp_boun_dir_length = pow(tmp_boun_dir_length, .5); 

				for(int k = 0; k < 3; ++ k) tmp_boun_dir[k] /= tmp_boun_dir_length; */

				double sign = 0.0; 
				for(int k = 0; k < 3; ++ k) sign += tmp_dir[k] * tmp_boun_dir[k]; 
				if(sign < 0.0) { tmp_boun_dir[0] *= -1; tmp_boun_dir[1] *= -1; tmp_boun_dir[2] *= -1; }

				double modu_v, modu_a; 
				modu_v = modulus_3_d(tmp_vel); 
				modu_a = modulus_3_d(tmp_acc); 

				double vel_prj, acc_prj; 
				vel_prj = dot_prod_3_d(tmp_vel, tmp_boun_dir);
				acc_prj = dot_prod_3_d(tmp_acc, tmp_boun_dir);

				vel_prj = (vel_prj <= 0.0) ? vel_prj : 0.0; 
				acc_prj = (acc_prj <= 0.0) ? acc_prj : 0.0; 

				if(vel_prj == 0.0) goto ACC;
				
				/* * * * * * * * * * * * * * * * * * * * * */

				double tmp_dis = dot_prod_3_d(tmp_dir, tmp_boun_dir); 

				double colli_time = tmp_dis / abs(vel_prj);

				double c0 = 1.0; 

				double _vel_boun_acc; 
				_vel_boun_acc = abs(vel_prj) / (c0 * colli_time); 

				_acc[0][0] += _vel_boun_acc * (tmp_boun_dir[0]); 
				_acc[0][1] += _vel_boun_acc * (tmp_boun_dir[1]); 
				_acc[0][2] += _vel_boun_acc * (tmp_boun_dir[2]); 

				/* * * * * * * * * * * * * * * * * * * * * */

ACC:
				if(acc_prj == 0.0) continue;

				double c1 = .93;

				double penetration_dis = 1.0 * (*particles)[calcu_particle].smthR.val[0] - dis;

				double  penetr_dis_ratio, acc_ratio; 

				penetr_dis_ratio = penetration_dis / (1.0 * (*particles)[calcu_particle].smthR.val[0]); 
				acc_ratio = penetr_dis_ratio / c1; 

				double _acc_boun_acc; 
				_acc_boun_acc = abs(acc_prj) * acc_ratio; 

				_acc[1][0] += _acc_boun_acc * (tmp_boun_dir[0]); 
				_acc[1][1] += _acc_boun_acc * (tmp_boun_dir[1]); 
				_acc[1][2] += _acc_boun_acc * (tmp_boun_dir[2]); 
			}
		}

		int tmp_acc_sign[3]; 

		tmp_acc_sign[0] = (_acc[1][0] >= 0.0) ? 1.0 : -1.0; 
		tmp_acc_sign[1] = (_acc[1][1] >= 0.0) ? 1.0 : -1.0; 
		tmp_acc_sign[2] = (_acc[1][2] >= 0.0) ? 1.0 : -1.0; 


		_acc[1][0] = constrain_d(abs(_acc[1][0]), 0.0, abs((*particles)[calcu_particle].accX.val[0])); _acc[1][0] *= tmp_acc_sign[0]; 
		_acc[1][1] = constrain_d(abs(_acc[1][1]), 0.0, abs((*particles)[calcu_particle].accY.val[0])); _acc[1][1] *= tmp_acc_sign[1];
		_acc[1][2] = constrain_d(abs(_acc[1][2]), 0.0, abs((*particles)[calcu_particle].accZ.val[0])); _acc[1][2] *= tmp_acc_sign[2];


		{
			for(int k = 0; k < 2; ++ k) acc_x += _acc[k][0]; 
			for(int k = 0; k < 2; ++ k) acc_y += _acc[k][1]; 
			for(int k = 0; k < 2; ++ k) acc_z += _acc[k][2]; 
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

		if(CASE_DIM == 0) (*particles)[calcu_particle].accZ.val[0] = 0.0;
	}
}