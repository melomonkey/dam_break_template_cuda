#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\basic_func\basic_func.h"
#include <vector>
using namespace std;

void velocity_filter1(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles)
{
	for(int i = 0; i < neighbors_list.size(); ++ i)
	{
		int calcu_particle = neighbors_list[i][0]; 

		double h = particles[calcu_particle].smthR.val[0]; 

		// The revisements of velocities in XYZ direction
		double tmp_velX_rev = 0.0,
		       tmp_velY_rev = 0.0,
			   tmp_velZ_rev = 0.0; 

		double epsilon = .2; 

		for(int j = 1; j < neighbors_list[i].size(); ++ j)
		{
			int label_ij = neighbors_list[i][j];

			double dis_x, dis_y, dis_z; 
			dis_x = particles[label_ij].coorX.val[0] - particles[calcu_particle].coorX.val[0]; 
			dis_y = particles[label_ij].coorY.val[0] - particles[calcu_particle].coorY.val[0]; 
			dis_z = particles[label_ij].coorZ.val[0] - particles[calcu_particle].coorZ.val[0];  

			double coeff0 = epsilon * particles[label_ij].mass.val[0];
			coeff0 /= (.5 * (particles[label_ij].density.val[0] + particles[calcu_particle].density.val[0]));

			double tmp_x, tmp_y, tmp_z; 
			tmp_x = coeff0 * (particles[label_ij].velX.val[0] - particles[calcu_particle].velX.val[0]) * kernel_function(dis_x, h); 
			tmp_y = coeff0 * (particles[label_ij].velY.val[0] - particles[calcu_particle].velY.val[0]) * kernel_function(dis_y, h); 
			tmp_z = coeff0 * (particles[label_ij].velZ.val[0] - particles[calcu_particle].velZ.val[0]) * kernel_function(dis_z, h); 

			tmp_velX_rev += tmp_x;
			tmp_velY_rev += tmp_y; 
			tmp_velZ_rev += tmp_z; 
		}

		particles[calcu_particle].velX.val[2] = particles[calcu_particle].velX.val[0] + tmp_velX_rev; 
		particles[calcu_particle].velY.val[2] = particles[calcu_particle].velY.val[0] + tmp_velY_rev; 
		particles[calcu_particle].velZ.val[2] = particles[calcu_particle].velZ.val[0] + tmp_velZ_rev; 
	}

	for(int i = 0; i < neighbors_list.size(); ++ i)
	{
		int calcu_particle = neighbors_list[i][0]; 

		particles[calcu_particle].velX.val[0] = particles[calcu_particle].velX.val[2]; 
		particles[calcu_particle].velY.val[0] = particles[calcu_particle].velY.val[2]; 
		particles[calcu_particle].velZ.val[0] = particles[calcu_particle].velZ.val[2]; 
	}
}

void velocity_filter1_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, int threads_num = 5)
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

			// The revisements of velocities in XYZ direction
			double tmp_velX_rev = 0.0,
				tmp_velY_rev = 0.0,
				tmp_velZ_rev = 0.0;

			double epsilon = .2;

			for (int j = 1; j < neighbors_list[i].size(); ++j)
			{
				int label_ij = neighbors_list[i][j];

				double dis_x, dis_y, dis_z;
				dis_x = particles[label_ij].coorX.val[0] - particles[calcu_particle].coorX.val[0];
				dis_y = particles[label_ij].coorY.val[0] - particles[calcu_particle].coorY.val[0];
				dis_z = particles[label_ij].coorZ.val[0] - particles[calcu_particle].coorZ.val[0];

				double coeff0 = epsilon * particles[label_ij].mass.val[0];
				coeff0 /= (.5 * (particles[label_ij].density.val[0] + particles[calcu_particle].density.val[0]));

				double tmp_x, tmp_y, tmp_z;
				tmp_x = coeff0 * (particles[label_ij].velX.val[0] - particles[calcu_particle].velX.val[0]) * kernel_function(dis_x, h);
				tmp_y = coeff0 * (particles[label_ij].velY.val[0] - particles[calcu_particle].velY.val[0]) * kernel_function(dis_y, h);
				tmp_z = coeff0 * (particles[label_ij].velZ.val[0] - particles[calcu_particle].velZ.val[0]) * kernel_function(dis_z, h);

				tmp_velX_rev += tmp_x;
				tmp_velY_rev += tmp_y;
				tmp_velZ_rev += tmp_z;
			}

			particles[calcu_particle].velX.val[2] = particles[calcu_particle].velX.val[0] + tmp_velX_rev;
			particles[calcu_particle].velY.val[2] = particles[calcu_particle].velY.val[0] + tmp_velY_rev;
			particles[calcu_particle].velZ.val[2] = particles[calcu_particle].velZ.val[0] + tmp_velZ_rev;
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

			particles[calcu_particle].velX.val[0] = particles[calcu_particle].velX.val[2];
			particles[calcu_particle].velY.val[0] = particles[calcu_particle].velY.val[2];
			particles[calcu_particle].velZ.val[0] = particles[calcu_particle].velZ.val[2];
		}
	}
}

void velocity_filter2(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles)
{
	for (int i = 0; i < neighbors_list.size(); ++i)
	{
		int calcu_particle = neighbors_list[i][0];

		double local_x, local_y, local_z;
		local_x = particles[calcu_particle].coorX.val[0];
		local_y = particles[calcu_particle].coorY.val[0];
		local_z = particles[calcu_particle].coorZ.val[0];
		double h = particles[calcu_particle].smthR.val[0];

		double _sum1_X = 0.0,
			_sum1_Y = 0.0,
			_sum1_Z = 0.0;

		double _sum2 = 0.0;

		for (int j = 0; j < neighbors_list[i].size(); ++j)
		{
			int label_ij = neighbors_list[i][j];

			if (particles[label_ij].id == 0 || particles[label_ij].id == 1 || particles[label_ij].id == 2 || particles[label_ij].id == 3) continue;

			double dis_x, dis_y, dis_z;
			dis_x = local_x - particles[label_ij].coorX.val[0];
			dis_y = local_y - particles[label_ij].coorY.val[0];
			dis_z = local_z - particles[label_ij].coorZ.val[0];

			double dis = pow(dis_x, 2.0) + pow(dis_y, 2.0) + pow(dis_z, 2.0); dis = pow(dis, .5);

			double kernel_val = kernel_function(dis, h);

			if (particles[label_ij].density.val[0] == 0.0) continue;

			// X direction
			_sum1_X += (particles[label_ij].velX.val[0]) * (particles[label_ij].mass.val[0] / particles[label_ij].density.val[0]) * kernel_val;

			// Y direction
			_sum1_Y += (particles[label_ij].velY.val[0]) * (particles[label_ij].mass.val[0] / particles[label_ij].density.val[0]) * kernel_val;

			// Z direction
			_sum1_Z += (particles[label_ij].velZ.val[0]) * (particles[label_ij].mass.val[0] / particles[label_ij].density.val[0]) * kernel_val;


			_sum2 += (particles[label_ij].mass.val[0] / particles[label_ij].density.val[0]) * kernel_val;
		}

		if (_sum2 == 0.0) continue;

		particles[calcu_particle].velX.val[2] = _sum1_X / _sum2;
		particles[calcu_particle].velY.val[2] = _sum1_Y / _sum2;
		particles[calcu_particle].velZ.val[2] = _sum1_Z / _sum2;
	}

	for (int i = 0; i < neighbors_list.size(); ++i)
	{
		int calcu_particle = neighbors_list[i][0];

		particles[calcu_particle].velX.val[0] = particles[calcu_particle].velX.val[2];
		particles[calcu_particle].velY.val[0] = particles[calcu_particle].velY.val[2];
		particles[calcu_particle].velZ.val[0] = particles[calcu_particle].velZ.val[2];
	}
}

void velocity_filter2_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, int threads_num = 5)
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

			double _sum1_X = 0.0,
				_sum1_Y = 0.0,
				_sum1_Z = 0.0;

			double _sum2 = 0.0;

			for (int j = 0; j < neighbors_list[i].size(); ++j)
			{
				int label_ij = neighbors_list[i][j];

				if (particles[label_ij].id == 0 || particles[label_ij].id == 1 || particles[label_ij].id == 2 || particles[label_ij].id == 3) continue;

				double dis_x, dis_y, dis_z;
				dis_x = local_x - particles[label_ij].coorX.val[0];
				dis_y = local_y - particles[label_ij].coorY.val[0];
				dis_z = local_z - particles[label_ij].coorZ.val[0];

				double dis = pow(dis_x, 2.0) + pow(dis_y, 2.0) + pow(dis_z, 2.0); dis = pow(dis, .5);

				double kernel_val = kernel_function(dis, h);

				if (particles[label_ij].density.val[0] == 0.0) continue;

				// X direction
				_sum1_X += (particles[label_ij].velX.val[0]) * (particles[label_ij].mass.val[0] / particles[label_ij].density.val[0]) * kernel_val;

				// Y direction
				_sum1_Y += (particles[label_ij].velY.val[0]) * (particles[label_ij].mass.val[0] / particles[label_ij].density.val[0]) * kernel_val;

				// Z direction
				_sum1_Z += (particles[label_ij].velZ.val[0]) * (particles[label_ij].mass.val[0] / particles[label_ij].density.val[0]) * kernel_val;


				_sum2 += (particles[label_ij].mass.val[0] / particles[label_ij].density.val[0]) * kernel_val;
			}

			if (_sum2 == 0.0) continue;

			particles[calcu_particle].velX.val[2] = _sum1_X / _sum2;
			particles[calcu_particle].velY.val[2] = _sum1_Y / _sum2;
			particles[calcu_particle].velZ.val[2] = _sum1_Z / _sum2;
		}
	}

	for (int i = 0; i < neighbors_list.size(); ++i)
	{
		int calcu_particle = neighbors_list[i][0];

		particles[calcu_particle].velX.val[0] = particles[calcu_particle].velX.val[2];
		particles[calcu_particle].velY.val[0] = particles[calcu_particle].velY.val[2];
		particles[calcu_particle].velZ.val[0] = particles[calcu_particle].velZ.val[2];
	}
}

void velocity_filter11(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, int threads_num = 8)
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

			double _sum1_X = 0.0,
				_sum1_Y = 0.0,
				_sum1_Z = 0.0;

			double _sum2 = 0.0;

			for (int j = 0; j < neighbors_list[i].size(); ++j)
			{
				int label_ij = neighbors_list[i][j];

				if (particles[label_ij].id == 0 || particles[label_ij].id == 1 || particles[label_ij].id == 2 || particles[label_ij].id == 3) continue;

				double dis_x, dis_y, dis_z;
				dis_x = local_x - particles[label_ij].coorX.val[0];
				dis_y = local_y - particles[label_ij].coorY.val[0];
				dis_z = local_z - particles[label_ij].coorZ.val[0];

				double kernel_val[3] = { kernel_function(dis_x, h), kernel_function(dis_y, h), kernel_function(dis_z, h) };

				if (particles[label_ij].density.val[0] == 0.0) continue;

				double tmp_volume = particles[label_ij].mass.val[0] / particles[label_ij].density.val[0];

				// X direction
				_sum1_X += (particles[label_ij].velX.val[0]) * tmp_volume * kernel_val[0];
				_sum1_X += (particles[label_ij].velX.val[0]) * tmp_volume * kernel_val[1];
				if (CASE_DIM == 1) _sum1_X += (particles[label_ij].velX.val[0]) * tmp_volume * kernel_val[2];

				// Y direction
				_sum1_Y += (particles[label_ij].velY.val[0]) * tmp_volume * kernel_val[0];
				_sum1_Y += (particles[label_ij].velY.val[0]) * tmp_volume * kernel_val[1];
				if (CASE_DIM == 1) _sum1_Y += (particles[label_ij].velY.val[0]) * tmp_volume * kernel_val[2];

				// Z direction
				_sum1_Z += (particles[label_ij].velZ.val[0]) * tmp_volume * kernel_val[0];
				_sum1_Z += (particles[label_ij].velZ.val[0]) * tmp_volume * kernel_val[1];
				_sum1_Z += (particles[label_ij].velZ.val[0]) * tmp_volume * kernel_val[2];


				_sum2 += tmp_volume * kernel_val[0];
				_sum2 += tmp_volume * kernel_val[1];
				if (CASE_DIM == 1) _sum2 += tmp_volume * kernel_val[2];
			}

			if (_sum2 == 0.0) continue;

			particles[calcu_particle].velX.val[2] = _sum1_X / _sum2;
			particles[calcu_particle].velY.val[2] = _sum1_Y / _sum2;
			particles[calcu_particle].velZ.val[2] = _sum1_Z / _sum2;
		}
	}

	for (int i = 0; i < neighbors_list.size(); ++i)
	{
		int calcu_particle = neighbors_list[i][0];

		particles[calcu_particle].velX.val[0] = particles[calcu_particle].velX.val[2];
		particles[calcu_particle].velY.val[0] = particles[calcu_particle].velY.val[2];
		particles[calcu_particle].velZ.val[0] = particles[calcu_particle].velZ.val[2];
	}
}