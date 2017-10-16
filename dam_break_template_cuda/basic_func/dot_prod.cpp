double dot_prod_3_d(double * v1, double * v2)
{
	double tmp = 0.0; 
	for(int i = 0; i < 3; ++ i) tmp += v1[i] * v2[i]; 
	return tmp; 
}

float dot_prod_3_f(float * v1, float * v2)
{
	float tmp = 0.0; 
	for(int i = 0; i < 3; ++ i) tmp += v1[i] * v2[i]; 
	return tmp; 
}

int dot_prod_3_i(int * v1, int * v2)
{
	int tmp = 0.0; 
	for(int i = 0; i < 3; ++ i) tmp += v1[i] * v2[i]; 
	return tmp; 
}

