#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__device__ double dot_prod_3_d_gpu(double * v1, double * v2)
{
	double tmp = 0.0;
	for (int i = 0; i < 3; ++i) tmp += v1[i] * v2[i];
	return tmp;
}

__device__ float dot_prod_3_f_gpu(float * v1, float * v2)
{
	float tmp = 0.0;
	for (int i = 0; i < 3; ++i) tmp += v1[i] * v2[i];
	return tmp;
}

__device__ int dot_prod_3_i_gpu(int * v1, int * v2)
{
	int tmp = 0.0;
	for (int i = 0; i < 3; ++i) tmp += v1[i] * v2[i];
	return tmp;
}