#include "struct_particle.h"
#include "global_variables.h"
#include <vector>
#include <omp.h>
#include "fast_search_neighbors\fast_search_neighbors.h"
using namespace std;

int grid_x, grid_y, grid_z;
float bins_factor = 1.3; 

void work_out_grids()
{
	float x_range, y_range, z_range;

	float mx_x = particles[0].coorX.val[0], mn_x = particles[0].coorX.val[0];
	float mx_y = particles[0].coorY.val[0], mn_y = particles[0].coorY.val[0];
	float mx_z = particles[0].coorZ.val[0], mn_z = particles[0].coorZ.val[0];

	float mx_radius = particles[0].smthR.val[0];

	float mean_radius = 0.0;
	int _cnt = 0; 

	for (int i = 0; i < particles.size(); ++i)
	{
		if (particles[i].id < 4) continue; 

		mean_radius += particles[i].smthR.val[0]; 

		if (particles[i].coorX.val[0] > mx_x) mx_x = particles[i].coorX.val[0];
		if (particles[i].coorX.val[0] < mn_x) mn_x = particles[i].coorX.val[0];

		if (particles[i].coorY.val[0] > mx_y) mx_y = particles[i].coorY.val[0];
		if (particles[i].coorY.val[0] < mn_y) mn_y = particles[i].coorY.val[0];

		if (particles[i].coorZ.val[0] > mx_z) mx_z = particles[i].coorZ.val[0];
		if (particles[i].coorZ.val[0] < mn_z) mn_z = particles[i].coorZ.val[0];

		_cnt++; 
	}

	mean_radius /= _cnt; 

	x_range = mx_x - mn_x;
	y_range = mx_y - mn_y;
	z_range = mx_z - mn_z;

	// adjustable factor
	float bin_size = bins_factor * mean_radius;

	grid_x = x_range / bin_size;
	grid_y = y_range / bin_size;
	grid_z = z_range / bin_size;

	if (grid_x < 1) grid_x = 1;
	if (grid_y < 1) grid_y = 1;
	if (grid_z < 1) grid_z = 1;
}

void neighbor_list_shrink(int threads_num = 5)
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

		for (int i = beg_idx; i < end_idx; ++i)
		{
			for (int j = 0; j < neighbors_list[i].size(); ++j)
			{
				int neigh_idx = neighbors_list[i][j];

				if (abs(particles[i].coorX.val[0] - particles[neigh_idx].coorX.val[0]) > 2.0 * particles[i].smthR.val[0])
				{
					neighbors_list[i].erase(neighbors_list[i].begin() + j);
					j--;
					continue;
				}
				if (abs(particles[i].coorY.val[0] - particles[neigh_idx].coorY.val[0]) > 2.0 * particles[i].smthR.val[0])
				{
					neighbors_list[i].erase(neighbors_list[i].begin() + j);
					j--;
					continue;
				}
				if (abs(particles[i].coorZ.val[0] - particles[neigh_idx].coorZ.val[0]) > 2.0 * particles[i].smthR.val[0])
				{
					neighbors_list[i].erase(neighbors_list[i].begin() + j);
					j--;
					continue;
				}
			}
		}
	}
}

void generate_neighbors_list()
{
	work_out_grids();

	vector<float> val_x, val_y, val_z;
	val_x.resize(particles.size());
	val_y.resize(particles.size());
	val_z.resize(particles.size());

	for (int i = 0; i < particles.size(); ++i)
	{
		val_x[i] = particles[i].coorX.val[0];
		val_y[i] = particles[i].coorY.val[0];
		val_z[i] = particles[i].coorZ.val[0];
	}

	fast_search_neighbors FTS;

	vector<vector<vector<vector<int>>>> list_3D;
	vector<vector<int>> label_3D;
	vector<vector<int>> neigh_list_3D;
	vector<vector<int>> neigh_list_3D_another;

	FTS.split_grid_3D_label(val_x, val_y, val_z, list_3D, label_3D, grid_x, grid_y, grid_z);

	FTS.search_neighbor_3D_multi_threads(list_3D, label_3D, neigh_list_3D, 1, 10);

	/* * * * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * */

	neighbors_list.clear();

	neighbors_list = neigh_list_3D;

	neighbor_list_shrink(10);

	for (int i = 0; i < neighbors_list.size(); ++i) neighbors_list[i].insert(neighbors_list[i].begin(), i); 


	for (int i = 0; i < neighbors_list.size(); ++i)
	{
		for (int j = 0; j < neighbors_list[i].size(); ++j)
			neighbors_list[i][j] ++; 
	}
}