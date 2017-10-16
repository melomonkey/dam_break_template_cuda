#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__device__ double modulus_3_d_gpu(double * v)
{
	double tmp = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	tmp = pow(tmp, .5);
	return tmp;
}

__device__ float modulus_3_f_gpu(float * v)
{
	float tmp = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	tmp = powf(tmp, .5);
	return tmp;
}

__device__ int modulus_3_i_gpu(int * v)
{
	float tmp = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	tmp = powf(tmp, .5);
	return tmp;
}