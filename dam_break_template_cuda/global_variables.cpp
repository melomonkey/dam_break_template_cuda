#include <vector>
#include "struct_particle.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
using namespace std;

vector<PARTICLE> particles;
unsigned int _size_particles;
PARTICLE * _particles_buff;
PARTICLE * _particles_cuda;
vector<vector<double>> coordination;
double * _coordination_buff; 
double * _coordination_cuda; 

vector<vector<int>> neighbors_list;
unsigned int _neigh_list_row, _neigh_list_col;
unsigned int * _neighbors_list_buff;
unsigned int * _neighbors_list_cuda;

vector<vector<int>> neighbors_list_wall;

double cur_calculation_time; 

cudaError_t cudaStatus;
