#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__device__ double constrain_d_gpu(double val, double min, double max)
{
	double tmp = val;

	if (tmp <= min) tmp = min; if (tmp >= max) tmp = max;

	return tmp;
}

__device__ float constrain_f_gpu(float val, float min, float max)
{
	float tmp = val;

	if (tmp <= min) tmp = min; if (tmp >= max) tmp = max;

	return tmp;
}

__device__ int constrain_i_gpu(int val, int min, int max)
{
	int tmp = val;

	if (tmp <= min) tmp = min; if (tmp >= max) tmp = max;

	return tmp;
}