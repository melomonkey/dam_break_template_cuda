
#ifndef basic_func_h_
#define basic_func_h_

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

double kernel_function(double dis, double h);
double kernel_function_1dev(double dis, double h);
double kernel_function_2dev(double dis, double h);

/* * * * * * * * * * * * * * * * * * * * * */

double constrain_d(double val, double min, double max);
float constrain_f(float val, float min, float max);
int constrain_i(int val, int min, int max);

/* * * * * * * * * * * * * * * * * * * * * */

double dot_prod_3_d(double * v1, double * v2);
float dot_prod_3_f(float * v1, float * v2);
int dot_prod_3_i(int * v1, int * v2);

/* * * * * * * * * * * * * * * * * * * * * */

void cross_prod_3_d(double * v1, double * v2, double * ans);
void cross_prod_3_f(float * v1, float * v2, float * ans);
void cross_prod_3_i(int * v1, int * v2, int * ans);

/* * * * * * * * * * * * * * * * * * * * * */

double modulus_3_d(double * v);
float modulus_3_f(float * v);
int modulus_3_i(int * v);

/* * * * * * * * * * * * * * * * * * * * * */

extern __device__ double kernel_function_gpu(double dis, double h);
extern __device__ double kernel_function_1dev_gpu(double dis, double h);
extern __device__ double kernel_function_2dev_gpu(double dis, double h);

/* * * * * * * * * * * * * * * * * * * * * */

extern __device__ double constrain_d_gpu(double val, double min, double max);
extern __device__ float constrain_f_gpu(float val, float min, float max);
extern __device__ int constrain_i_gpu(int val, int min, int max);

/* * * * * * * * * * * * * * * * * * * * * */

extern __device__ double dot_prod_3_d_gpu(double * v1, double * v2);
extern __device__ float dot_prod_3_f_gpu(float * v1, float * v2);
extern __device__ int dot_prod_3_i_gpu(int * v1, int * v2);

/* * * * * * * * * * * * * * * * * * * * * */

extern __device__ void cross_prod_3_d_gpu(double * v1, double * v2, double * ans);
extern __device__ void cross_prod_3_f_gpu(float * v1, float * v2, float * ans);
extern __device__ void cross_prod_3_i_gpu(int * v1, int * v2, int * ans);

/* * * * * * * * * * * * * * * * * * * * * */

extern __device__ double modulus_3_d_gpu(double * v);
extern __device__ float modulus_3_f_gpu(float * v);
extern __device__ int modulus_3_i_gpu(int * v);

/* * * * * * * * * * * * * * * * * * * * * */



#endif 