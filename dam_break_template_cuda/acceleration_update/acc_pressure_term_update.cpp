#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\basic_func\basic_func.h"
#include <vector>
using namespace std;

void acc_pressure_term_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles)
{
	for (int i = 0; i < neighbors_list.size(); ++i)
	{
		int calcu_particle = neighbors_list[i][0];

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

		for (int j = 1; j < neighbors_list[i].size(); ++j)
		{
			int label_ij = neighbors_list[i][j];

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
				double de_kernel_val[3] = { kernel_function_1dev(dis_x, h), kernel_function_1dev(dis_y, h), kernel_function_1dev(dis_z, h) };

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
			acc_x = constrain_d(acc_x, ACC_X_LIM_L, ACC_X_LIM_U);
			acc_y = constrain_d(acc_y, ACC_Y_LIM_L, ACC_Y_LIM_U);
			acc_z = constrain_d(acc_z, ACC_Z_LIM_L, ACC_Z_LIM_U);
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

}

void acc_pressure_term_update_schm1_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, int threads_num = 5)
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

			double acc_x = 0.0,
				acc_y = 0.0,
				acc_z = 0.0;

			double _acc[3];
			for (int m = 0; m < 3; ++m) _acc[m] = 0.0;

			for (int j = 1; j < neighbors_list[i].size(); ++j)
			{
				int label_ij = neighbors_list[i][j];

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
					coeff0 = -(particles[label_ij].mass.val[0] * (particles[calcu_particle].pressure.val[0] + REF_PRESSURE + particles[label_ij].pressure.val[0] + REF_PRESSURE) / (particles[calcu_particle].density.val[0] * particles[label_ij].density.val[0]));

					double de_kernel_val[3] = { kernel_function_1dev(dis_x, h), kernel_function_1dev(dis_y, h), kernel_function_1dev(dis_z, h) };

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
				acc_x = constrain_d(acc_x, ACC_X_LIM_L, ACC_X_LIM_U);
				acc_y = constrain_d(acc_y, ACC_Y_LIM_L, ACC_Y_LIM_U);
				acc_z = constrain_d(acc_z, ACC_Z_LIM_L, ACC_Z_LIM_U);
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
	}
}

void acc_pressure_term_update_schm2(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, double dt)
{
	for (int i = 0; i < neighbors_list.size(); ++i)
	{
		int calcu_particle = neighbors_list[i][0];

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

		for (int j = 1; j < neighbors_list[i].size(); ++j)
		{
			int label_ij = neighbors_list[i][j];

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
				coeff0 = ((particles[calcu_particle].pressure.val[0] + REF_PRESSURE) / pow(particles[calcu_particle].density.val[0], 2.0)) + ((particles[label_ij].pressure.val[0] + REF_PRESSURE) / pow(particles[label_ij].density.val[0], 2.0));
				coeff0 = -particles[label_ij].mass.val[0] * coeff0;

				double de_kernel_val[3] = { kernel_function_1dev(dis_x, h), kernel_function_1dev(dis_y, h), kernel_function_1dev(dis_z, h) };

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
			acc_x = constrain_d(acc_x, ACC_X_LIM_L, ACC_X_LIM_U);
			acc_y = constrain_d(acc_y, ACC_Y_LIM_L, ACC_Y_LIM_U);
			acc_z = constrain_d(acc_z, ACC_Z_LIM_L, ACC_Z_LIM_U);
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

}

void acc_pressure_term_update_schm2_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, double dt, int threads_num = 5)
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

			double acc_x = 0.0,
				acc_y = 0.0,
				acc_z = 0.0;

			double _acc[3];
			for (int m = 0; m < 3; ++m) _acc[m] = 0.0;

			for (int j = 1; j < neighbors_list[i].size(); ++j)
			{
				int label_ij = neighbors_list[i][j];

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
					coeff0 = ((particles[calcu_particle].pressure.val[0] + REF_PRESSURE) / pow(particles[calcu_particle].density.val[0], 2.0)) + ((particles[label_ij].pressure.val[0] + REF_PRESSURE) / pow(particles[label_ij].density.val[0], 2.0));
					coeff0 = -particles[label_ij].mass.val[0] * coeff0;

					double de_kernel_val[3] = { kernel_function_1dev(dis_x, h), kernel_function_1dev(dis_y, h), kernel_function_1dev(dis_z, h) };

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
				acc_x = constrain_d(acc_x, ACC_X_LIM_L, ACC_X_LIM_U);
				acc_y = constrain_d(acc_y, ACC_Y_LIM_L, ACC_Y_LIM_U);
				acc_z = constrain_d(acc_z, ACC_Z_LIM_L, ACC_Z_LIM_U);
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
	}
}


//void acc_pressure_term_update_schm22232(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles, double dt)
//{
//	for(int i = 0; i < neighbors_list.size(); ++ i)
//	{
//		int calcu_particle = neighbors_list[i][0];
//
//		double local_x, local_y, local_z;
//		local_x = (*particles)[calcu_particle].coorX.val[0]; 
//		local_y = (*particles)[calcu_particle].coorY.val[0];
//		local_z = (*particles)[calcu_particle].coorZ.val[0];
//		double h = (*particles)[calcu_particle].smthR.val[0];
//
//		double acc_x = 0.0, 
//			   acc_y = 0.0, 
//			   acc_z = 0.0; 
//
//		/*
//		*
//		* _acc is a matrix stores the accelerations caused in 3 differernt way.
//		*
//		* 1. The first line of the matrix represents the accelerations caused by the interaction between the general particles.
//		* 2. The second line of the matrix represents the accelerations caused by the momentum bounced back from wall.
//		* 3. The third line of the matrix represents the accelerations caused by the replusive force generated by the wall. They will be never bigger than the accelerations the particle pocesses.
//		*
//		*/
//
//		double _acc[3][3]; 
//		for(int m = 0; m < 3; ++ m) { for(int n = 0; n < 3; ++ n) _acc[m][n] = 0.0; }
//
//		for(int j = 1; j < neighbors_list[i].size(); ++ j)
//		{
//			int label_ij = neighbors_list[i][j];
//
//			double tmp_x, tmp_y, tmp_z; 
//			tmp_x = (*particles)[label_ij].coorX.val[0]; 
//			tmp_y = (*particles)[label_ij].coorY.val[0];
//			tmp_z = (*particles)[label_ij].coorZ.val[0];
//
//			double dis_x, dis_y, dis_z; 
//			double dis_sqr, dis; 
//			dis_x = local_x - tmp_x; 
//			dis_y = local_y - tmp_y; 
//			dis_z = local_z - tmp_z; 
//			dis_sqr = pow(dis_x, 2.0) + pow(dis_y, 2.0) + pow(dis_z, 2.0); 
//			dis = pow(dis_sqr, .5); 
//
//			double r_0 = (*particles)[label_ij].smthR.val[0]; 
//
//			// Interactive particles
//			if((*particles)[label_ij].id >= 4)
//			{
//				double coeff0; 
//				coeff0 = -((*particles)[label_ij].mass.val[0] * ((*particles)[calcu_particle].pressure.val[0] + REF_PRESSURE + (*particles)[label_ij].pressure.val[0] + REF_PRESSURE) / ((*particles)[calcu_particle].density.val[0] * (*particles)[label_ij].density.val[0])); 
//
//				double coeff1; 			
//				coeff1 = (*particles)[label_ij].mass.val[0] * ((*particles)[calcu_particle].viscosity.val[0] + (*particles)[label_ij].viscosity.val[0]); 
//				coeff1 /= ((*particles)[calcu_particle].density.val[0] * (*particles)[label_ij].density.val[0]);
//
//				double de_kernael_val_dis;
//				de_kernael_val_dis = kernel_function_1dev(dis, h);
//
//				if(dis == 0.0) dis = 1e70; 
//
//				double de_kernel_val;
//
//				de_kernel_val = kernel_function_1dev(dis_x, h); 		
//				_acc[0][0] += coeff0 * de_kernel_val + coeff1 * (de_kernael_val_dis * (1.0 / dis)) * ((*particles)[calcu_particle].velX.val[0] - (*particles)[label_ij].velX.val[0]); 
//
//				de_kernel_val = kernel_function_1dev(dis_y, h);
//				_acc[0][1] += coeff0 * de_kernel_val + coeff1 * (de_kernael_val_dis * (1.0 / dis)) * ((*particles)[calcu_particle].velY.val[0] - (*particles)[label_ij].velY.val[0]);
//
//				de_kernel_val = kernel_function_1dev(dis_z, h);
//				_acc[0][2] += coeff0 * de_kernel_val + coeff1 * (de_kernael_val_dis * (1.0 / dis)) * ((*particles)[calcu_particle].velZ.val[0] - (*particles)[label_ij].velZ.val[0]);
//
//			}
//
//
//			// Wall particles
//			if((*particles)[label_ij].id < 4)
//			{
//				double tmp_vel[3]; 
//				tmp_vel[0] = (*particles)[calcu_particle].velX.val[0] - (*particles)[label_ij].velX.val[0]; 
//				tmp_vel[1] = (*particles)[calcu_particle].velY.val[0] - (*particles)[label_ij].velY.val[0]; 
//				tmp_vel[2] = (*particles)[calcu_particle].velZ.val[0] - (*particles)[label_ij].velZ.val[0]; 
//
//				double tmp_acc[3]; 
//				tmp_acc[0] = (*particles)[calcu_particle].accX.val[0] - (*particles)[label_ij].accX.val[0]; 
//				tmp_acc[1] = (*particles)[calcu_particle].accY.val[0] - (*particles)[label_ij].accY.val[0];
//				tmp_acc[2] = (*particles)[calcu_particle].accZ.val[0] - (*particles)[label_ij].accZ.val[0]; 
//
//				double tmp_dir[3];
//				tmp_dir[0] = dis_x; 
//				tmp_dir[1] = dis_y; 
//				tmp_dir[2] = dis_z; 
//
//
//
//				double modu_v, modu_a, modu_n; 
//				modu_v = modulus_3_d(tmp_vel); 
//				modu_a = modulus_3_d(tmp_acc); 
//				modu_n = modulus_3_d(tmp_dir); 
//
//				double vel_prj, acc_prj; 
//				vel_prj = dot_prod_3_d(tmp_vel, tmp_dir) / modu_n;
//				acc_prj = dot_prod_3_d(tmp_acc, tmp_dir) / modu_n;
//
//				vel_prj = (vel_prj <= 0.0) ? vel_prj : 0.0; 
//				acc_prj = (acc_prj <= 0.0) ? acc_prj : 0.0; 
//
//
//				/* * * * * * * * * * * * * * * * * * * * * */
//
//				double colli_time = modu_n / abs(vel_prj);
//
//				double c0 = .8; 
//
//				double _vel_boun_acc; 
//				_vel_boun_acc = abs(vel_prj) / (c0 * colli_time); 
//
//				_acc[1][0] += _vel_boun_acc * (tmp_dir[0] / modu_n); 
//				_acc[1][1] += _vel_boun_acc * (tmp_dir[1] / modu_n); 
//				_acc[1][2] += _vel_boun_acc * (tmp_dir[2] / modu_n); 
//
//				/* * * * * * * * * * * * * * * * * * * * * */
//				
//				double c1 = 0.1;
//
//				double penetration_dis = 2.0 * (*particles)[calcu_particle].smthR.val[0] - dis;
//
//				double  penetr_dis_ratio, acc_ratio; 
//
//				penetr_dis_ratio = penetration_dis / (2.0 * (*particles)[calcu_particle].smthR.val[0]); 
//				acc_ratio = penetr_dis_ratio / c1; 
//
//				double _acc_boun_acc; 
//				_acc_boun_acc = abs(acc_prj) * acc_ratio; 
//
//				_acc[2][0] += _acc_boun_acc * (tmp_dir[0] / modu_n); 
//				_acc[2][1] += _acc_boun_acc * (tmp_dir[1] / modu_n); 
//				_acc[2][2] += _acc_boun_acc * (tmp_dir[2] / modu_n); 
//			}
//		}
//
//		int tmp_acc_sign[3]; 
//
//		tmp_acc_sign[0] = (_acc[2][0] >= 0.0) ? 1.0 : -1.0; 
//		tmp_acc_sign[1] = (_acc[2][1] >= 0.0) ? 1.0 : -1.0; 
//		tmp_acc_sign[2] = (_acc[2][2] >= 0.0) ? 1.0 : -1.0; 
//
//
//		_acc[2][0] = constrain_d(abs(_acc[2][0]), 0.0, abs((*particles)[calcu_particle].accX.val[0])); _acc[2][0] *= tmp_acc_sign[0]; 
//		_acc[2][1] = constrain_d(abs(_acc[2][1]), 0.0, abs((*particles)[calcu_particle].accY.val[0])); _acc[2][1] *= tmp_acc_sign[1];
//		_acc[2][2] = constrain_d(abs(_acc[2][2]), 0.0, abs((*particles)[calcu_particle].accZ.val[0])); _acc[2][2] *= tmp_acc_sign[2];
//
//		// Body force
//		{
//			double volume_acc[3]; 
//
//			if(CASE_DIM == 0) { volume_acc[0] = 0.0; volume_acc[1] = -GRAVITY; volume_acc[2] = 0.0; };
//			if(CASE_DIM == 1) { volume_acc[0] = 0.0; volume_acc[1] = 0.0; volume_acc[2] = -GRAVITY; }; 
//
//			acc_x += volume_acc[0]; 
//			acc_y += volume_acc[1]; 
//			acc_z += volume_acc[2]; 
//		}
//
//
//		{
//			for(int k = 0; k < 3; ++ k) acc_x += _acc[k][0]; 
//			for(int k = 0; k < 3; ++ k) acc_y += _acc[k][1]; 
//			for(int k = 0; k < 3; ++ k) acc_z += _acc[k][2]; 
//		}
//
//		/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//
//		{ 
//			acc_x = constrain_d(acc_x, ACC_X_LIM_L, ACC_X_LIM_U); 
//			acc_y = constrain_d(acc_y, ACC_Y_LIM_L, ACC_Y_LIM_U); 
//			acc_z = constrain_d(acc_z, ACC_Z_LIM_L, ACC_Z_LIM_U); 			
//		}
//
//		/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//
//
//
//		(*particles)[calcu_particle].accX.val[2] = (*particles)[calcu_particle].accX.val[1]; 
//		(*particles)[calcu_particle].accX.val[1] = (*particles)[calcu_particle].accX.val[0];
//		(*particles)[calcu_particle].accX.val[0] = acc_x;
//
//		(*particles)[calcu_particle].accY.val[2] = (*particles)[calcu_particle].accY.val[1]; 
//		(*particles)[calcu_particle].accY.val[1] = (*particles)[calcu_particle].accY.val[0];
//		(*particles)[calcu_particle].accY.val[0] = acc_y; 
//
//		(*particles)[calcu_particle].accZ.val[2] = (*particles)[calcu_particle].accZ.val[1]; 
//		(*particles)[calcu_particle].accZ.val[1] = (*particles)[calcu_particle].accZ.val[0];
//		(*particles)[calcu_particle].accZ.val[0] = acc_z; 
//	}
//
//}
//
//void acceleration_update_schm2(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles, double dt)
//{
//	for(int i = 0; i < neighbors_list.size(); ++ i)
//	{
//		int calcu_particle = neighbors_list[i][0];
//
//		double local_x, local_y, local_z;
//		local_x = (*particles)[calcu_particle].coorX.val[0]; 
//		local_y = (*particles)[calcu_particle].coorY.val[0];
//		local_z = (*particles)[calcu_particle].coorZ.val[0];
//		double h = (*particles)[calcu_particle].smthR.val[0];
//
//		double acc_x = 0.0, 
//			   acc_y = 0.0, 
//			   acc_z = 0.0; 
//
//		/*
//		*
//		* _acc is a matrix stores the accelerations caused in 3 differernt way.
//		*
//		* 1. The first line of the matrix represents the accelerations caused by the interaction between the general particles.
//		* 2. The second line of the matrix represents the accelerations caused by the momentum bounced back from wall.
//		* 3. The third line of the matrix represents the accelerations caused by the replusive force generated by the wall. They will be never bigger than the accelerations the particle pocesses.
//		*
//		*/
//
//		double _acc[3][3]; 
//		for(int m = 0; m < 3; ++ m) { for(int n = 0; n < 3; ++ n) _acc[m][n] = 0.0; }
//
//		for(int j = 1; j < neighbors_list[i].size(); ++ j)
//		{
//			int label_ij = neighbors_list[i][j];
//
//			double tmp_x, tmp_y, tmp_z; 
//			tmp_x = (*particles)[label_ij].coorX.val[0]; 
//			tmp_y = (*particles)[label_ij].coorY.val[0];
//			tmp_z = (*particles)[label_ij].coorZ.val[0];
//
//			double dis_x, dis_y, dis_z; 
//			double dis_sqr, dis; 
//			dis_x = local_x - tmp_x; 
//			dis_y = local_y - tmp_y; 
//			dis_z = local_z - tmp_z; 
//			dis_sqr = pow(dis_x, 2.0) + pow(dis_y, 2.0) + pow(dis_z, 2.0); 
//			dis = pow(dis_sqr, .5); 
//
//			double r_0 = (*particles)[label_ij].smthR.val[0]; 
//
//			// Interactive particles
//			if((*particles)[label_ij].id >= 4)
//			{
//				double coeff0; 
//                coeff0 = (((*particles)[calcu_particle].pressure.val[0] + REF_PRESSURE) / pow((*particles)[calcu_particle].density.val[0], 2.0)) + (((*particles)[label_ij].pressure.val[0] + REF_PRESSURE) / pow((*particles)[label_ij].density.val[0], 2.0)); 
//				coeff0 = -(*particles)[label_ij].mass.val[0] * coeff0; 
//				double coeff1; 			
//				coeff1 = (*particles)[label_ij].mass.val[0] * ((*particles)[calcu_particle].viscosity.val[0] + (*particles)[label_ij].viscosity.val[0]); 
//				coeff1 /= ((*particles)[calcu_particle].density.val[0] * (*particles)[label_ij].density.val[0]);
//
//				double de_kernael_val_dis;
//				de_kernael_val_dis = kernel_function_1dev(dis, h);
//
//				if(dis == 0.0) dis = 1e70; 
//
//				double de_kernel_val;
//
//				de_kernel_val = kernel_function_1dev(dis_x, h) / h; 		
//				_acc[0][0] += coeff0 * de_kernel_val + coeff1 * (de_kernael_val_dis * (1.0 / dis)) * ((*particles)[calcu_particle].velX.val[0] - (*particles)[label_ij].velX.val[0]); 
//
//				de_kernel_val = kernel_function_1dev(dis_y, h) / h;
//				_acc[0][1] += coeff0 * de_kernel_val + coeff1 * (de_kernael_val_dis * (1.0 / dis)) * ((*particles)[calcu_particle].velY.val[0] - (*particles)[label_ij].velY.val[0]);
//
//				de_kernel_val = kernel_function_1dev(dis_z, h) / h;
//				_acc[0][2] += coeff0 * de_kernel_val + coeff1 * (de_kernael_val_dis * (1.0 / dis)) * ((*particles)[calcu_particle].velZ.val[0] - (*particles)[label_ij].velZ.val[0]);
//
//			}
//
//
//			// Wall particles
//			if((*particles)[label_ij].id < 4)
//			{
//				double tmp_vel[3]; 
//				tmp_vel[0] = (*particles)[calcu_particle].velX.val[0] - (*particles)[label_ij].velX.val[0]; 
//				tmp_vel[1] = (*particles)[calcu_particle].velY.val[0] - (*particles)[label_ij].velY.val[0]; 
//				tmp_vel[2] = (*particles)[calcu_particle].velZ.val[0] - (*particles)[label_ij].velZ.val[0]; 
//
//				double tmp_acc[3]; 
//				tmp_acc[0] = (*particles)[calcu_particle].accX.val[0] - (*particles)[label_ij].accX.val[0]; 
//				tmp_acc[1] = (*particles)[calcu_particle].accY.val[0] - (*particles)[label_ij].accY.val[0];
//				tmp_acc[2] = (*particles)[calcu_particle].accZ.val[0] - (*particles)[label_ij].accZ.val[0]; 
//
//				double tmp_dir[3];
//				tmp_dir[0] = dis_x; 
//				tmp_dir[1] = dis_y; 
//				tmp_dir[2] = dis_z; 
//
//
//
//				double modu_v, modu_a, modu_n; 
//				modu_v = modulus_3_d(tmp_vel); 
//				modu_a = modulus_3_d(tmp_acc); 
//				modu_n = modulus_3_d(tmp_dir); 
//
//				double vel_prj, acc_prj; 
//				vel_prj = dot_prod_3_d(tmp_vel, tmp_dir) / modu_n;
//				acc_prj = dot_prod_3_d(tmp_acc, tmp_dir) / modu_n;
//
//				vel_prj = (vel_prj <= 0.0) ? vel_prj : 0.0; 
//				acc_prj = (acc_prj <= 0.0) ? acc_prj : 0.0; 
//
//
//				/* * * * * * * * * * * * * * * * * * * * * */
//
//				double colli_time = modu_n / abs(vel_prj);
//
//				double c0 = .8; 
//
//				double _vel_boun_acc; 
//				_vel_boun_acc = abs(vel_prj) / (c0 * colli_time); 
//
//				_acc[1][0] += _vel_boun_acc * (tmp_dir[0] / modu_n); 
//				_acc[1][1] += _vel_boun_acc * (tmp_dir[1] / modu_n); 
//				_acc[1][2] += _vel_boun_acc * (tmp_dir[2] / modu_n); 
//
//				/* * * * * * * * * * * * * * * * * * * * * */
//				
//				double c1 = 0.1;
//
//				double penetration_dis = 2.0 * (*particles)[calcu_particle].smthR.val[0] - dis;
//
//				double  penetr_dis_ratio, acc_ratio; 
//
//				penetr_dis_ratio = penetration_dis / (2.0 * (*particles)[calcu_particle].smthR.val[0]); 
//				acc_ratio = penetr_dis_ratio / c1; 
//
//				double _acc_boun_acc; 
//				_acc_boun_acc = abs(acc_prj) * acc_ratio; 
//
//				_acc[2][0] += _acc_boun_acc * (tmp_dir[0] / modu_n); 
//				_acc[2][1] += _acc_boun_acc * (tmp_dir[1] / modu_n); 
//				_acc[2][2] += _acc_boun_acc * (tmp_dir[2] / modu_n); 
//			}
//		}
//
//		int tmp_acc_sign[3]; 
//
//		tmp_acc_sign[0] = (_acc[2][0] >= 0.0) ? 1.0 : -1.0; 
//		tmp_acc_sign[1] = (_acc[2][1] >= 0.0) ? 1.0 : -1.0; 
//		tmp_acc_sign[2] = (_acc[2][2] >= 0.0) ? 1.0 : -1.0; 
//
//
//		_acc[2][0] = constrain_d(abs(_acc[2][0]), 0.0, abs((*particles)[calcu_particle].accX.val[0])); _acc[2][0] *= tmp_acc_sign[0]; 
//		_acc[2][1] = constrain_d(abs(_acc[2][1]), 0.0, abs((*particles)[calcu_particle].accY.val[0])); _acc[2][1] *= tmp_acc_sign[1];
//		_acc[2][2] = constrain_d(abs(_acc[2][2]), 0.0, abs((*particles)[calcu_particle].accZ.val[0])); _acc[2][2] *= tmp_acc_sign[2];
//
//		// Body force
//		{
//			double volume_acc[3]; 
//
//			if(CASE_DIM == 0) { volume_acc[0] = 0.0; volume_acc[1] = -GRAVITY; volume_acc[2] = 0.0; };
//			if(CASE_DIM == 1) { volume_acc[0] = 0.0; volume_acc[1] = 0.0; volume_acc[2] = -GRAVITY; }; 
//
//			acc_x += volume_acc[0]; 
//			acc_y += volume_acc[1]; 
//			acc_z += volume_acc[2]; 
//		}
//
//
//		{
//			for(int k = 0; k < 3; ++ k) acc_x += _acc[k][0]; 
//			for(int k = 0; k < 3; ++ k) acc_y += _acc[k][1]; 
//			for(int k = 0; k < 3; ++ k) acc_z += _acc[k][2]; 
//		}
//
//		/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//
//		{ 
//			acc_x = constrain_d(acc_x, ACC_X_LIM_L, ACC_X_LIM_U); 
//			acc_y = constrain_d(acc_y, ACC_Y_LIM_L, ACC_Y_LIM_U); 
//			acc_z = constrain_d(acc_z, ACC_Z_LIM_L, ACC_Z_LIM_U); 			
//		}
//
//		/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//
//
//
//		(*particles)[calcu_particle].accX.val[2] = (*particles)[calcu_particle].accX.val[1]; 
//		(*particles)[calcu_particle].accX.val[1] = (*particles)[calcu_particle].accX.val[0];
//		(*particles)[calcu_particle].accX.val[0] = acc_x;
//
//		(*particles)[calcu_particle].accY.val[2] = (*particles)[calcu_particle].accY.val[1]; 
//		(*particles)[calcu_particle].accY.val[1] = (*particles)[calcu_particle].accY.val[0];
//		(*particles)[calcu_particle].accY.val[0] = acc_y; 
//
//		(*particles)[calcu_particle].accZ.val[2] = (*particles)[calcu_particle].accZ.val[1]; 
//		(*particles)[calcu_particle].accZ.val[1] = (*particles)[calcu_particle].accZ.val[0];
//		(*particles)[calcu_particle].accZ.val[0] = acc_z; 
//	}
//
//}
//
//void acceleration_update_schm3(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles, double dt)
//{
//	for(int i = 0; i < neighbors_list.size(); ++ i)
//	{
//		int calcu_particle = neighbors_list[i][0];
//
//		double local_x, local_y, local_z;
//		local_x = (*particles)[calcu_particle].coorX.val[0]; 
//		local_y = (*particles)[calcu_particle].coorY.val[0];
//		local_z = (*particles)[calcu_particle].coorZ.val[0];
//		double h = (*particles)[calcu_particle].smthR.val[0];
//
//		double acc_x = 0.0, 
//			   acc_y = 0.0, 
//			   acc_z = 0.0; 
//
//		for(int j = 1; j < neighbors_list[i].size(); ++ j)
//		{
//			int label_ij = neighbors_list[i][j];
//
//			double tmp_x, tmp_y, tmp_z; 
//			tmp_x = (*particles)[label_ij].coorX.val[0]; 
//			tmp_y = (*particles)[label_ij].coorY.val[0];
//			tmp_z = (*particles)[label_ij].coorZ.val[0];
//
//			double dis_x, dis_y, dis_z; 
//			double dis_sqr, dis; 
//			dis_x = local_x - tmp_x; 
//			dis_y = local_y - tmp_y; 
//			dis_z = local_z - tmp_z; 
//			dis_sqr = pow(dis_x, 2.0) + pow(dis_y, 2.0) + pow(dis_z, 2.0); 
//			dis = pow(dis_sqr, .5); 
//
//			double r_0 = (*particles)[label_ij].smthR.val[0]; 
//
//			double _acc[3][3]; 
//			for(int m = 0; m < 3; ++ m)
//			{
//				for(int n = 0; n < 3; ++ n)
//					_acc[m][n] = 0.0; 
//			}
//
//			// 压力和粘性力
//			{
//				double coeff0; 
//				coeff0 = -((*particles)[label_ij].mass.val[0] * ((*particles)[calcu_particle].pressure.val[0] + REF_PRESSURE + (*particles)[label_ij].pressure.val[0] + REF_PRESSURE) / ((*particles)[calcu_particle].density.val[0] * (*particles)[label_ij].density.val[0])); 
//
//				double coeff1; 			
//				coeff1 = (*particles)[label_ij].mass.val[0] * ((*particles)[calcu_particle].viscosity.val[0] + (*particles)[label_ij].viscosity.val[0]); 
//				coeff1 /= ((*particles)[calcu_particle].density.val[0] * (*particles)[label_ij].density.val[0]);
//
//				double de_kernael_val_dis;
//				de_kernael_val_dis = kernel_function_1dev(dis, h);
//
//				if(dis == 0.0) dis = 1e70; 
//
//				double de_kernel_val;
//				de_kernel_val = kernel_function_1dev(dis_x, h); 		
//				_acc[0][0] = coeff0 * de_kernel_val + coeff1 * (de_kernael_val_dis * (1.0 / dis)) * ((*particles)[calcu_particle].velX.val[0] - (*particles)[label_ij].velX.val[0]); 
//
//				de_kernel_val = kernel_function_1dev(dis_y, h);
//				_acc[0][1] = coeff0 * de_kernel_val + coeff1 * (de_kernael_val_dis * (1.0 / dis)) * ((*particles)[calcu_particle].velY.val[0] - (*particles)[label_ij].velY.val[0]);
//
//				de_kernel_val = kernel_function_1dev(dis_z, h);
//				_acc[0][2] = coeff0 * de_kernel_val + coeff1 * (de_kernael_val_dis * (1.0 / dis)) * ((*particles)[calcu_particle].velZ.val[0] - (*particles)[label_ij].velZ.val[0]);
//			}
//
//			// 壁面粒子
//			{
//				if(r_0 / dis < 1.0)
//				{
//					double n1 = 12.0; double n2 = 6.0;
//					_acc[1][0] = -(WALL_BOUNCE_CONST * (pow((r_0 / dis), n1) - pow((r_0 / dis), n2)) * dis_x) / dis_sqr; 
//					_acc[1][1] = -(WALL_BOUNCE_CONST * (pow((r_0 / dis), n1) - pow((r_0 / dis), n2)) * dis_y) / dis_sqr; 
//					_acc[1][2] = -(WALL_BOUNCE_CONST * (pow((r_0 / dis), n1) - pow((r_0 / dis), n2)) * dis_z) / dis_sqr; 
//				}
//			}
//
//			
//			// 相应型号对应边界虚粒子更新规则
//			{
//			
//				if((*particles)[label_ij].id == 0) { acc_x += _acc[1][0]; acc_y += _acc[1][1]; acc_z += _acc[1][2]; }
//				if((*particles)[label_ij].id == 1) { acc_x += _acc[1][0]; }
//				if((*particles)[label_ij].id == 2) { acc_y += _acc[1][1]; }
//				if((*particles)[label_ij].id == 3) { acc_z += _acc[1][2]; }
//				if((*particles)[label_ij].id == 4) { acc_x += _acc[0][0]; acc_y += _acc[0][1]; acc_z += _acc[0][2]; }
//			}
//		}
//
//		// 体积力
//		{
//			double volume_acc[3] = { 0.0, -GRAVITY, 0.0 }; 
//
//			acc_x += volume_acc[0]; 
//			acc_y += volume_acc[1]; 
//			acc_z += volume_acc[2]; 
//		}
//
//
//
//		/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//
//		{ 
//			acc_x = constrain_d(acc_x, ACC_X_LIM_L, ACC_X_LIM_U); 
//			acc_y = constrain_d(acc_y, ACC_Y_LIM_L, ACC_Y_LIM_U); 
//			acc_z = constrain_d(acc_z, ACC_Z_LIM_L, ACC_Z_LIM_U); 			
//		}
//
//		/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//
//
//
//		(*particles)[calcu_particle].accX.val[2] = (*particles)[calcu_particle].accX.val[1]; 
//		(*particles)[calcu_particle].accX.val[1] = (*particles)[calcu_particle].accX.val[0];
//		(*particles)[calcu_particle].accX.val[0] = acc_x;
//
//		(*particles)[calcu_particle].accY.val[2] = (*particles)[calcu_particle].accY.val[1]; 
//		(*particles)[calcu_particle].accY.val[1] = (*particles)[calcu_particle].accY.val[0];
//		(*particles)[calcu_particle].accY.val[0] = acc_y; 
//
//		(*particles)[calcu_particle].accZ.val[2] = (*particles)[calcu_particle].accZ.val[1]; 
//		(*particles)[calcu_particle].accZ.val[1] = (*particles)[calcu_particle].accZ.val[0];
//		(*particles)[calcu_particle].accZ.val[0] = acc_z; 
//	}
//
//}