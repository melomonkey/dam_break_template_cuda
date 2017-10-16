#include <math.h>

double modulus_3_d(double * v)
{
	double tmp = v[0] * v[0] + v[1] * v[1] + v[2] * v[2]; 
	tmp = pow(tmp, .5); 
	return tmp; 
}

float modulus_3_f(float * v)
{
	float tmp = v[0] * v[0] + v[1] * v[1] + v[2] * v[2]; 
	tmp = powf(tmp, .5); 
	return tmp; 
}

int modulus_3_i(int * v)
{
	float tmp = v[0] * v[0] + v[1] * v[1] + v[2] * v[2]; 
	tmp = powf(tmp, .5); 
	return tmp; 
}

