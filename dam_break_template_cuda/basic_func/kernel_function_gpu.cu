#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <math.h>
#include "..\configuration.h"

__device__ double kernel_function_gpu(double dis, double h)
{
	double R = dis / h;
	double alpha;

	if (CASE_DIM == 0) alpha = (1.0 / (h * h)) * pow(PI, -1.0);
	if (CASE_DIM == 1) alpha = (1.0 / (h * h * h)) * pow(PI, -1.5);

	double ans = alpha * exp(-(R * R));
	return ans;
}

__device__ double kernel_function_1dev_gpu(double dis, double h)
{
	double R = dis / h;
	double alpha;

	if (CASE_DIM == 0) alpha = (1.0 / (h * h)) * pow(3.14, -1.0);
	if (CASE_DIM == 1) alpha = (1.0 / (h * h * h)) * pow(3.14, -1.5);

	double ans = alpha * exp(-(R * R)) * (-2.0 * R);

	return ans / h;
}

__device__ double kernel_function_2dev_gpu(double dis, double h)
{
	double R = dis / h;
	double alpha;

	if (CASE_DIM == 0) alpha = (1.0 / (h * h)) * pow(PI, -1.0);
	if (CASE_DIM == 1) alpha = (1.0 / (h * h * h)) * pow(PI, -1.5);

	double ans = -2.0 * alpha * (1 - 2 * R * R) * exp(-(R * R));

	return ans / h / h; // local particle's position coordination's derivative, so include the R's derivative, twice.
}