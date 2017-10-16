#include <vector>
#include "fast_search_neighbors.h"
using namespace std;

void fast_search_neighbors::split_grid_3D_label(vector<float> val_x, vector<float> val_y, vector<float> val_z, vector<vector<vector<vector<int>>>> & grid_list, vector<vector<int>> & label, int split_bins_x = 8, int split_bins_y = 8, int split_bins_z = 8)
{
	float mx_X, mn_X;
	mx_X = val_x[0];
	mn_X = val_x[0];

	for (int i = 0; i < val_x.size(); ++i)
	{
		if (val_x[i] > mx_X) mx_X = val_x[i];
		if (val_x[i] < mn_X) mn_X = val_x[i];
	}

	float interval_X = (mx_X - mn_X) / split_bins_x;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * */

	float mx_Y, mn_Y;
	mx_Y = val_y[0];
	mn_Y = val_y[0];

	for (int i = 0; i < val_y.size(); ++i)
	{
		if (val_y[i] > mx_Y) mx_Y = val_y[i];
		if (val_y[i] < mn_Y) mn_Y = val_y[i];
	}

	float interval_Y = (mx_Y - mn_Y) / split_bins_y;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * */

	float mx_Z, mn_Z;
	mx_Z = val_z[0];
	mn_Z = val_z[0];

	for (int i = 0; i < val_z.size(); ++i)
	{
		if (val_z[i] > mx_Z) mx_Z = val_z[i];
		if (val_z[i] < mn_Z) mn_Z = val_z[i];
	}

	float interval_Z = (mx_Z - mn_Z) / split_bins_z;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * */

	grid_list.resize(split_bins_x);
	for (int i = 0; i < split_bins_x; ++i) grid_list[i].resize(split_bins_y);

	for (int i = 0; i < split_bins_x; ++i)
	{
		for (int j = 0; j < split_bins_y; ++j) grid_list[i][j].resize(split_bins_z);
	}

	for (int i = 0; i < split_bins_x; ++i)
	{
		for (int j = 0; j < split_bins_y; ++j)
		{
			for (int k = 0; k < split_bins_z; ++k)
			{
				grid_list[i][j][k].reserve(val_x.size() / (split_bins_x * split_bins_y * split_bins_z));
			}
		}
	}

	label.resize(val_x.size());
	for (int i = 0; i < label.size(); ++i)
	{
		vector<int> tmp; tmp.resize(3);
		label[i] = tmp;
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * */

	for (int i = 0; i < val_x.size(); ++i)
	{
		int loc_X = (val_x[i] - mn_X) / (interval_X * 1.0);

		if (loc_X > split_bins_x - 1) loc_X = split_bins_x - 1;
		if (loc_X < 0) loc_X = 0;

		/* * * * * * * * * * * * * * * * * */

		int loc_Y = (val_y[i] - mn_Y) / (interval_Y * 1.0);

		if (loc_Y > split_bins_y - 1) loc_Y = split_bins_y - 1;
		if (loc_Y < 0) loc_Y = 0;

		/* * * * * * * * * * * * * * * * * */

		int loc_Z = (val_z[i] - mn_Z) / (interval_Z * 1.0);

		if (loc_Z > split_bins_z - 1) loc_Z = split_bins_z - 1;
		if (loc_Z < 0) loc_Z = 0;

		/* * * * * * * * * * * * * * * * * */

		grid_list[loc_X][loc_Y][loc_Z].push_back(i);

		label[i][0] = loc_X; label[i][1] = loc_Y; label[i][2] = loc_Z;
	}
}

void fast_search_neighbors::search_neighbor_3D(vector<vector<vector<vector<int>>>> grid_list, vector<vector<int>> label, vector<vector<int>> & neigh_list, int search_offset = 1)
{
	if (search_offset < 0) search_offset = 0;

	int x_up_bound = grid_list.size() - 1;
	int y_up_bound = grid_list[0].size() - 1;
	int z_up_bound = grid_list[0][0].size() - 1;

	neigh_list.resize(label.size());

	for (int i = 0; i < label.size(); ++i)
	{
		int idx_x, idx_y, idx_z;
		idx_x = label[i][0];
		idx_y = label[i][1];
		idx_z = label[i][2];

		int beg_x, end_x, beg_y, end_y, beg_z, end_z;
		beg_x = idx_x - search_offset; if (beg_x < 0) beg_x = 0;
		end_x = idx_x + search_offset; if (end_x > x_up_bound) end_x = x_up_bound;

		beg_y = idx_y - search_offset; if (beg_y < 0) beg_y = 0;
		end_y = idx_y + search_offset; if (end_y > y_up_bound) end_y = y_up_bound;

		beg_z = idx_z - search_offset; if (beg_z < 0) beg_z = 0;
		end_z = idx_z + search_offset; if (end_z > z_up_bound) end_z = z_up_bound;

		for (int ii = beg_x; ii <= end_x; ++ii)
		{
			for (int jj = beg_y; jj <= end_y; ++jj)
			{
				for (int kk = beg_z; kk <= end_z; ++kk)
				{
					neigh_list[i].insert(neigh_list[i].end(), grid_list[ii][jj][kk].begin(), grid_list[ii][jj][kk].end());
				}
			}
		}
	}
}

void fast_search_neighbors::search_neighbor_3D_multi_threads(vector<vector<vector<vector<int>>>> grid_list, vector<vector<int>> label, vector<vector<int>> & neigh_list, int search_offset = 1, int threads_num = 5)
{
	if (search_offset < 0) search_offset = 0;

	int x_up_bound = grid_list.size() - 1;
	int y_up_bound = grid_list[0].size() - 1;
	int z_up_bound = grid_list[0][0].size() - 1;

	neigh_list.resize(label.size());

	int _size = label.size();
	int _thread_cnt = _size / threads_num + 1;

#pragma omp parallel for
	for (int iii = 0; iii < threads_num; ++iii)
	{
		int _beg, _end;
		_beg = iii * _thread_cnt;
		_end = (iii + 1) * _thread_cnt;

		if (_end > label.size() - 1) _end = label.size() - 1;

		if (search_offset < 0) search_offset = 0;

		int x_up_bound = grid_list.size() - 1;
		int y_up_bound = grid_list[0].size() - 1;
		int z_up_bound = grid_list[0][0].size() - 1;

		for (int i = _beg; i <= _end; ++i)
		{
			int idx_x, idx_y, idx_z;
			idx_x = label[i][0];
			idx_y = label[i][1];
			idx_z = label[i][2];

			int beg_x, end_x, beg_y, end_y, beg_z, end_z;
			beg_x = idx_x - search_offset; if (beg_x < 0) beg_x = 0;
			end_x = idx_x + search_offset; if (end_x > x_up_bound) end_x = x_up_bound;

			beg_y = idx_y - search_offset; if (beg_y < 0) beg_y = 0;
			end_y = idx_y + search_offset; if (end_y > y_up_bound) end_y = y_up_bound;

			beg_z = idx_z - search_offset; if (beg_z < 0) beg_z = 0;
			end_z = idx_z + search_offset; if (end_z > z_up_bound) end_z = z_up_bound;

			for (int ii = beg_x; ii <= end_x; ++ii)
			{
				for (int jj = beg_y; jj <= end_y; ++jj)
				{
					for (int kk = beg_z; kk <= end_z; ++kk)
					{
						neigh_list[i].insert(neigh_list[i].end(), grid_list[ii][jj][kk].begin(), grid_list[ii][jj][kk].end());
					}
				}
			}
		}
	}
}

/* * * * * * * * * * * * * * * * * * */

void fast_search_neighbors::split_grid_2D_label(vector<float> val_x, vector<float> val_y, vector<vector<vector<int>>> & grid_list, vector<vector<int>> & label, int split_bins_x = 8, int split_bins_y = 8)
{
	float mx_X, mn_X;
	mx_X = val_x[0];
	mn_X = val_x[0];

	for (int i = 0; i < val_x.size(); ++i)
	{
		if (val_x[i] > mx_X) mx_X = val_x[i];
		if (val_x[i] < mn_X) mn_X = val_x[i];
	}

	float interval_X = (mx_X - mn_X) / split_bins_x;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * */

	float mx_Y, mn_Y;
	mx_Y = val_y[0];
	mn_Y = val_y[0];

	for (int i = 0; i < val_y.size(); ++i)
	{
		if (val_y[i] > mx_Y) mx_Y = val_y[i];
		if (val_y[i] < mn_Y) mn_Y = val_y[i];
	}

	float interval_Y = (mx_Y - mn_Y) / split_bins_y;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * */

	grid_list.resize(split_bins_x);
	for (int i = 0; i < split_bins_x; ++i) grid_list[i].resize(split_bins_y);
	for (int i = 0; i < split_bins_x; ++i)
	{
		for (int j = 0; j < split_bins_y; ++j) grid_list[i][j].reserve(val_x.size() / (split_bins_x * split_bins_y));
	}

	label.resize(val_x.size());
	for (int i = 0; i < label.size(); ++i)
	{
		vector<int> tmp; tmp.resize(2);
		label[i] = tmp;
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * */

	for (int i = 0; i < val_x.size(); ++i)
	{
		int loc_X = (val_x[i] - mn_X) / (interval_X * 1.0);

		if (loc_X > split_bins_x - 1) loc_X = split_bins_x - 1;
		if (loc_X < 0) loc_X = 0;

		/* * * * * * * * * * * * * * * * * */

		int loc_Y = (val_y[i] - mn_Y) / (interval_Y * 1.0);

		if (loc_Y > split_bins_y - 1) loc_Y = split_bins_y - 1;
		if (loc_Y < 0) loc_Y = 0;

		/* * * * * * * * * * * * * * * * * */

		grid_list[loc_X][loc_Y].push_back(i);

		label[i][0] = loc_X; label[i][1] = loc_Y;
	}
}

void fast_search_neighbors::search_neighbor_2D(vector<vector<vector<int>>> grid_list, vector<vector<int>> label, vector<vector<int>> & neigh_list, int search_offset = 1)
{
	if (search_offset < 0) search_offset = 0;

	int x_up_bound = grid_list.size() - 1;
	int y_up_bound = grid_list[0].size() - 1;

	neigh_list.resize(label.size());

	for (int i = 0; i < label.size(); ++i)
	{
		int idx_x, idx_y;
		idx_x = label[i][0];
		idx_y = label[i][1];

		int beg_x, end_x, beg_y, end_y, beg_z, end_z;
		beg_x = idx_x - search_offset; if (beg_x < 0) beg_x = 0;
		end_x = idx_x + search_offset; if (end_x > x_up_bound) end_x = x_up_bound;

		beg_y = idx_y - search_offset; if (beg_y < 0) beg_y = 0;
		end_y = idx_y + search_offset; if (end_y > y_up_bound) end_y = y_up_bound;

		for (int ii = beg_x; ii <= end_x; ++ii)
		{
			for (int jj = beg_y; jj <= end_y; ++jj)
			{
				neigh_list[i].insert(neigh_list[i].end(), grid_list[ii][jj].begin(), grid_list[ii][jj].end());
			}
		}
	}
}

void fast_search_neighbors::search_neighbor_2D_multi_threads(vector<vector<vector<int>>> grid_list, vector<vector<int>> label, vector<vector<int>> & neigh_list, int search_offset = 1, int threads_num = 5)
{
	if (search_offset < 0) search_offset = 0;

	int x_up_bound = grid_list.size() - 1;
	int y_up_bound = grid_list[0].size() - 1;

	neigh_list.resize(label.size());

	int _size = label.size();
	int _thread_cnt = _size / threads_num + 1;

#pragma omp parallel for
	for (int iii = 0; iii < threads_num; ++iii)
	{
		int _beg, _end;
		_beg = iii * _thread_cnt;
		_end = (iii + 1) * _thread_cnt;

		if (_end > label.size() - 1) _end = label.size() - 1;

		if (search_offset < 0) search_offset = 0;

		int x_up_bound = grid_list.size() - 1;
		int y_up_bound = grid_list[0].size() - 1;

		for (int i = _beg; i <= _end; ++i)
		{
			int idx_x, idx_y;
			idx_x = label[i][0];
			idx_y = label[i][1];

			int beg_x, end_x, beg_y, end_y, beg_z, end_z;
			beg_x = idx_x - search_offset; if (beg_x < 0) beg_x = 0;
			end_x = idx_x + search_offset; if (end_x > x_up_bound) end_x = x_up_bound;

			beg_y = idx_y - search_offset; if (beg_y < 0) beg_y = 0;
			end_y = idx_y + search_offset; if (end_y > y_up_bound) end_y = y_up_bound;

			for (int ii = beg_x; ii <= end_x; ++ii)
			{
				for (int jj = beg_y; jj <= end_y; ++jj)
				{
					neigh_list[i].insert(neigh_list[i].end(), grid_list[ii][jj].begin(), grid_list[ii][jj].end());
				}
			}
		}
	}
}

/* * * * * * * * * * * * * * * * * * */

void fast_search_neighbors::split_grid_1D_label(vector<float> val, vector<vector<int>> & list, vector<int> & label, int split_bins = 8)
{
	float mx, mn;
	mx = val[0];
	mn = val[0];

	for (int i = 0; i < val.size(); ++i)
	{
		if (val[i] > mx) mx = val[i];
		if (val[i] < mn) mn = val[i];
	}

	float interval = (mx - mn) / split_bins;

	list.resize(split_bins);
	for (int i = 0; i < split_bins; ++i) list[i].reserve(val.size() / split_bins);

	label.resize(val.size());

	for (int i = 0; i < val.size(); ++i)
	{
		int loc = (val[i] - mn) / (interval * 1.0);

		if (loc > split_bins - 1) loc = split_bins - 1;
		if (loc < 0) loc = 0;

		list[loc].push_back(i);

		label[i] = loc;
	}
}

void fast_search_neighbors::search_neighbor_1D(vector<vector<int>> grid_list, vector<int> label, vector<vector<int>> & neigh_list, int search_offset = 1)
{
	if (search_offset < 0) search_offset = 0;

	int x_up_bound = grid_list.size() - 1;

	neigh_list.resize(label.size());

	for (int i = 0; i < label.size(); ++i)
	{
		int idx_x;
		idx_x = label[i];

		int beg_x, end_x, beg_y, end_y, beg_z, end_z;
		beg_x = idx_x - search_offset; if (beg_x < 0) beg_x = 0;
		end_x = idx_x + search_offset; if (end_x > x_up_bound) end_x = x_up_bound;

		for (int ii = beg_x; ii <= end_x; ++ii)
		{
			neigh_list[i].insert(neigh_list[i].end(), grid_list[ii].begin(), grid_list[ii].end());
		}
	}
}

void fast_search_neighbors::search_neighbor_1D_multi_threads(vector<vector<int>> grid_list, vector<int> label, vector<vector<int>> & neigh_list, int search_offset = 1, int threads_num = 5)
{
	if (search_offset < 0) search_offset = 0;

	int x_up_bound = grid_list.size() - 1;

	neigh_list.resize(label.size());

	int _size = label.size();
	int _thread_cnt = _size / threads_num + 1;

#pragma omp parallel for
	for (int iii = 0; iii < threads_num; ++iii)
	{
		int _beg, _end;
		_beg = iii * _thread_cnt;
		_end = (iii + 1) * _thread_cnt;

		if (_end > label.size() - 1) _end = label.size() - 1;

		if (search_offset < 0) search_offset = 0;

		int x_up_bound = grid_list.size() - 1;

		for (int i = _beg; i <= _end; ++i)
		{
			int idx_x;
			idx_x = label[i];

			int beg_x, end_x, beg_y, end_y, beg_z, end_z;
			beg_x = idx_x - search_offset; if (beg_x < 0) beg_x = 0;
			end_x = idx_x + search_offset; if (end_x > x_up_bound) end_x = x_up_bound;

			for (int ii = beg_x; ii <= end_x; ++ii)
			{
				neigh_list[i].insert(neigh_list[i].end(), grid_list[ii].begin(), grid_list[ii].end());
			}
		}
	}
}
