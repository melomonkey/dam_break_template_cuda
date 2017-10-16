double constrain_d(double val, double min, double max)
{
	double tmp = val; 

	if(tmp <= min) tmp = min; if(tmp >= max) tmp = max; 

	return tmp; 
}

float constrain_f(float val, float min, float max)
{
	float tmp = val; 

	if(tmp <= min) tmp = min; if(tmp >= max) tmp = max; 

	return tmp; 
}

int constrain_i(int val, int min, int max)
{
	int tmp = val; 

	if(tmp <= min) tmp = min; if(tmp >= max) tmp = max; 

	return tmp; 
}


