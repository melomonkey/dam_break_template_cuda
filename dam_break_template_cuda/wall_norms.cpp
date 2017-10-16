#include "struct_particle.h"
#include "configuration.h"
#include "global_variables.h"
#include "basic_func.h"
#include <vector>
using namespace std;

struct VEC
{
	double x; double y; double z; double length; 
}; 

void normal_vec_initialize2D()
{
	void generate_neighbors_list_wall(); 
	generate_neighbors_list_wall(); 

	for(int i = 0; i < neighbors_list_wall.size(); ++ i)
	{
		vector<int> eff_particle; 
		for(int j = 0; j < neighbors_list_wall[i].size(); ++ j)
		{
			if(particles[neighbors_list_wall[i][j]].id >= 4) continue; 
			eff_particle.push_back(neighbors_list_wall[i][j]); 
		}

		vector<VEC> norm_vecs; 
		for(int j = 1; j < eff_particle.size(); ++ j)
		{
			int loc_idx = eff_particle[0], 
				tgr_idx = eff_particle[j]; 

			VEC tmp; 
			tmp.x = particles[loc_idx].coorX.val[0] - particles[tgr_idx].coorX.val[0]; 
			tmp.y = particles[loc_idx].coorY.val[0] - particles[tgr_idx].coorY.val[0]; 
			tmp.z = particles[loc_idx].coorZ.val[0] - particles[tgr_idx].coorZ.val[0];
			tmp.length = tmp.x * tmp.x + tmp.y * tmp.y; 
			tmp.length = pow(tmp.length, .5); 

			// rotate 90 degrees
			double tmp_x, tmp_y; 
			tmp_x = tmp.y; tmp_y = -tmp.x; 
			tmp.x = tmp_x; tmp.y = tmp_y;

			norm_vecs.push_back(tmp); 
		}

		if(norm_vecs.size() == 0) continue; 

		for(int j = 1; j < norm_vecs.size(); ++ j)
		{
			double tmp = norm_vecs[0].x * norm_vecs[j].x + norm_vecs[0].y * norm_vecs[j].y; 
			if(tmp < 0) { norm_vecs[j].x *= -1; norm_vecs[j].y *= -1;  }
		}

		double length_sum = 0.0, 
			   length_avg = 0.0;
		for(int j = 0; j < norm_vecs.size(); ++ j)
		{
			length_sum += norm_vecs[j].length; 
		}
		length_avg = length_sum / norm_vecs.size(); 

		VEC norm; 
		norm.x = 0.0; norm.y = 0.0; norm.z = 0.0; norm.length = 0.0;  

		for(int j = 0; j < norm_vecs.size(); ++ j) 
		{
			double weight = (length_avg + (length_avg - norm_vecs[j].length)) / length_sum; 
			norm.x += weight * norm_vecs[j].x; 
			norm.y += weight * norm_vecs[j].y; 
		}

		double tmp_L  = norm.x * norm.x + norm.y * norm.y; 
		tmp_L = pow(tmp_L, .5); 

		norm.x /= tmp_L; norm.y /= tmp_L; 

		particles[neighbors_list_wall[i][0]].normal_vec[0] = norm.x; 
		particles[neighbors_list_wall[i][0]].normal_vec[1] = norm.y; 
		particles[neighbors_list_wall[i][0]].normal_vec[2] = 0.0; 
	}
}

void normal_vec_initialize3D()
{
	void generate_neighbors_list_wall(); 
	generate_neighbors_list_wall(); 

	for(int i = 0; i < neighbors_list_wall.size(); ++ i)
	{
		vector<int> eff_particle; 
		for(int j = 0; j < neighbors_list_wall[i].size(); ++ j)
		{
			if(particles[neighbors_list_wall[i][j]].id >= 4) continue; 
			eff_particle.push_back(neighbors_list_wall[i][j]); 
		}

		vector<VEC> dir_vecs; 
		for(int j = 1; j < eff_particle.size(); ++ j)
		{
			int loc_idx = eff_particle[0], 
				tgr_idx = eff_particle[j]; 

			VEC tmp; 
			tmp.x = particles[loc_idx].coorX.val[0] - particles[tgr_idx].coorX.val[0]; 
			tmp.y = particles[loc_idx].coorY.val[0] - particles[tgr_idx].coorY.val[0]; 
			tmp.z = particles[loc_idx].coorZ.val[0] - particles[tgr_idx].coorZ.val[0];
			tmp.length = tmp.x * tmp.x + tmp.y * tmp.y + tmp.z * tmp.z; 
			tmp.length = pow(tmp.length, .5); 


			dir_vecs.push_back(tmp); 
		}

		if(dir_vecs.size() == 1) continue; 

		if(dir_vecs.size() * 0.8 > 2)
		{
			int end = 0.8 * dir_vecs.size(); 

			vector<VEC> dir_vecs_; dir_vecs_.assign(dir_vecs.begin(), dir_vecs.begin() + end); 
		    dir_vecs = dir_vecs_; 
		}

		vector<VEC> norm_vecs;
		for(int j = 0; j < dir_vecs.size(); ++ j)
		{
			for(int k = j + 1; k < dir_vecs.size(); ++ k)
			{
				double v1[3], v2[3], ans[3]; 
				v1[0] = dir_vecs[j].x; v1[1] = dir_vecs[j].y; v1[2] = dir_vecs[j].z;
				v2[0] = dir_vecs[k].x; v2[1] = dir_vecs[k].y; v2[2] = dir_vecs[k].z;

				cross_prod_3_d(v1, v2, ans); 

				VEC tmp; 
				tmp.x = ans[0]; tmp.y = ans[1]; tmp.z = ans[2]; 
				tmp.length = tmp.x * tmp.x + tmp.y * tmp.y + tmp.z * tmp.z; 
			    tmp.length = pow(tmp.length, .5);  

				norm_vecs.push_back(tmp); 
			}
		}

		if(norm_vecs.size() == 0) continue; 

		double length_sum = 0.0, 
			   length_avg = 0.0;
		for(int j = 0; j < norm_vecs.size(); ++ j)
		{
			length_sum += norm_vecs[j].length; 
		}
		length_avg = length_sum / norm_vecs.size(); 

		VEC norm; 
		norm.x = 0.0; norm.y = 0.0; norm.z = 0.0; norm.length = 0.0;  

		for(int j = 0; j < norm_vecs.size(); ++ j) 
		{
			double weight = (length_avg + (length_avg - norm_vecs[j].length)) / length_sum; 
			norm.x += weight * norm_vecs[j].x; 
			norm.y += weight * norm_vecs[j].y; 
			norm.z += weight * norm_vecs[j].z; 
		}

		double tmp_L  = norm.x * norm.x + norm.y * norm.y + norm.z * norm.z; 
		tmp_L = pow(tmp_L, .5); 

		norm.x /= tmp_L; norm.y /= tmp_L; norm.z /= tmp_L; 

		particles[neighbors_list_wall[i][0]].normal_vec[0] = norm.x; 
		particles[neighbors_list_wall[i][0]].normal_vec[1] = norm.y; 
		particles[neighbors_list_wall[i][0]].normal_vec[2] = norm.z; 
	}
}