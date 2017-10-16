#pragma once

#ifndef FAST_SEARCH_NEIGHBORS_H_ 
#define FAST_SEARCH_NEIGHBORS_H_ 

#include <vector>
#include <random>
#include <omp.h>
using namespace std;

class fast_search_neighbors
{
public:

	void split_grid_3D_label(vector<float> val_x, vector<float> val_y, vector<float> val_z, vector<vector<vector<vector<int>>>> & grid_list, vector<vector<int>> & label, int split_bins_x, int split_bins_y, int split_bins_z);
	
	void search_neighbor_3D(vector<vector<vector<vector<int>>>> grid_list, vector<vector<int>> label, vector<vector<int>> & neigh_list, int search_offset);

	void search_neighbor_3D_multi_threads(vector<vector<vector<vector<int>>>> grid_list, vector<vector<int>> label, vector<vector<int>> & neigh_list, int search_offset, int threads_num);

	/* * * * * * * * * * * * * * * * * * */

	void split_grid_2D_label(vector<float> val_x, vector<float> val_y, vector<vector<vector<int>>> & grid_list, vector<vector<int>> & label, int split_bins_x, int split_bins_y);

	void search_neighbor_2D(vector<vector<vector<int>>> grid_list, vector<vector<int>> label, vector<vector<int>> & neigh_list, int search_offset);

	void search_neighbor_2D_multi_threads(vector<vector<vector<int>>> grid_list, vector<vector<int>> label, vector<vector<int>> & neigh_list, int search_offset, int threads_num);

	/* * * * * * * * * * * * * * * * * * */

	void split_grid_1D_label(vector<float> val, vector<vector<int>> & list, vector<int> & label, int split_bins);

	void search_neighbor_1D(vector<vector<int>> grid_list, vector<int> label, vector<vector<int>> & neigh_list, int search_offset);

	void search_neighbor_1D_multi_threads(vector<vector<int>> grid_list, vector<int> label, vector<vector<int>> & neigh_list, int search_offset, int threads_num);
};

/* * * USAGE * * */

/*













*/

#endif
