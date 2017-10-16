#ifndef global_variables_h_
#define global_variables_h_

#include <vector>
#include "struct_particle.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
using std::vector;

extern vector<PARTICLE> particles;
extern unsigned int _size_particles;
extern PARTICLE * _particles_buff;
extern PARTICLE * _particles_cuda;
extern vector<vector<double>> coordination;
extern double * _coordination_buff;
extern double * _coordination_cuda;

extern vector<vector<int>> neighbors_list;
extern unsigned int _neigh_list_row, _neigh_list_col;
extern unsigned int * _neighbors_list_buff;
extern unsigned int * _neighbors_list_cuda;

extern vector<vector<int>> neighbors_list_wall;

extern double cur_calculation_time;

extern cudaError_t cudaStatus; 

#endif