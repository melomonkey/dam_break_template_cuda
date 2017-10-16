#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "struct_particle.h"
#include <stdio.h>
#include "global_variables.h"

void _neighbors_list_IN_gpu()
{
	unsigned int mx_idx = neighbors_list[0].size();

	for (int i = 0; i < neighbors_list.size(); ++i)
	{
		if (neighbors_list[i].size() > mx_idx) mx_idx = neighbors_list[i].size();
	}

	_neigh_list_row = neighbors_list.size();
	_neigh_list_col = mx_idx;

	_neighbors_list_buff = new unsigned int[_neigh_list_row * _neigh_list_col]; 
	for (int i = 0; i < _neigh_list_row * _neigh_list_col; ++i) _neighbors_list_buff[i] = 0; 

	for (int i = 0; i < neighbors_list.size(); ++i)
	{
		for (int j = 0; j < neighbors_list[i].size(); ++j)
		{
			_neighbors_list_buff[j + i * _neigh_list_col] = neighbors_list[i][j]; 
		}
	}

	cudaFree(_neighbors_list_cuda);
	cudaStatus = cudaMalloc((void **)&_neighbors_list_cuda, sizeof(unsigned int) * _neigh_list_row * _neigh_list_col);
	cudaStatus = cudaMemcpy(_neighbors_list_cuda, _neighbors_list_buff, sizeof(unsigned int) * _neigh_list_row * _neigh_list_col, cudaMemcpyHostToDevice);
	
	free(_neighbors_list_buff); 
}

void _coordination_list_OUT_gpu()
{
	
	_coordination_buff = new double[3 * _size_particles];

	cudaStatus = cudaMemcpy(_coordination_buff, _coordination_cuda, sizeof(double) * 3 * _size_particles, cudaMemcpyDeviceToHost);

	coordination.clear(); 
	coordination.reserve(particles.size()); 
	for (int i = 0; i < particles.size(); ++i)
	{
		vector<double> tmp(3, 0.0);
		tmp[0] = _coordination_buff[0 + i * 3];
		tmp[1] = _coordination_buff[1 + i * 3];
		tmp[2] = _coordination_buff[2 + i * 3];

		coordination.push_back(tmp);
	}

	free(_coordination_buff);
}

void _initialize_lists()
{
	_size_particles = particles.size(); 

	_particles_buff = new PARTICLE[_size_particles]; 

	cudaStatus = cudaMalloc((void **)&_particles_cuda, sizeof(PARTICLE) * _size_particles);

	cudaStatus = cudaMalloc((void **)&_coordination_cuda, sizeof(double) * 3 * _size_particles);

	cudaStatus = cudaMalloc((void **)&_neighbors_list_cuda, sizeof(unsigned int) * 2);

	coordination.clear();
	coordination.reserve(particles.size());
	for (int i = 0; i < particles.size(); ++i)
	{
		vector<double> tmp(3, 0.0);
		tmp[0] = particles[i].coorX.val[0];
		tmp[1] = particles[i].coorY.val[0];
		tmp[2] = particles[i].coorZ.val[0];

		coordination.push_back(tmp);
	}
}

void _particles_IN_gpu()
{
	for (int i = 0; i < particles.size(); ++i) _particles_buff[i] = particles[i];

	cudaStatus = cudaMalloc((void **)&_particles_cuda, sizeof(PARTICLE) * _size_particles);
	cudaStatus = cudaMemcpy(_particles_cuda, _particles_buff, sizeof(PARTICLE) * _size_particles, cudaMemcpyHostToDevice);

	
	/*__global__ void update_coordination_list_gpu(unsigned int _size, double * _coordination_list, PARTICLE * _particles);
	dim3 grid(1, 1, 1);
	dim3 thrd(1, 1, 1);

	update_coordination_list_gpu <<<grid, thrd>>> (_size_particles, _coordination_cuda, _particles_cuda);*/
}

void particles_OUT_gpu()
{
	cudaStatus = cudaMemcpy(_particles_buff, _particles_cuda, sizeof(PARTICLE) * _size_particles, cudaMemcpyDeviceToHost);

	for (int i = 0; i < particles.size(); ++i) particles[i] = _particles_buff[i]; 
}

__global__ void _update_coordination_list_gpu(unsigned int _size, double * _coordination_list_cuda, PARTICLE * _particles)
{
	for (int i = 0; i < _size; ++i)
	{
		_coordination_list_cuda[0 + i * 3] = _particles[i].coorX.val[0];
		_coordination_list_cuda[1 + i * 3] = _particles[i].coorY.val[0];
		_coordination_list_cuda[2 + i * 3] = _particles[i].coorZ.val[0];
	}
}

void update_coordination_list_gpu()
{
	dim3 grid(1, 1, 1);
	dim3 thrd(1, 1, 1);

	_update_coordination_list_gpu << <grid, thrd >> > (_size_particles, _coordination_cuda, _particles_cuda);
}