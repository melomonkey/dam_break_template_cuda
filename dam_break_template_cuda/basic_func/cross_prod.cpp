void cross_prod_3_d(double * v1, double * v2, double * ans)
{
	ans[0] = v1[1] * v2[2] - v1[2] * v2[1]; 
	ans[1] = v1[0] * v2[2] - v1[2] * v2[0];  
	ans[2] = v1[0] * v2[1] - v1[1] * v2[0];  
}

void cross_prod_3_f(float * v1, float * v2, float * ans)
{
	ans[0] = v1[1] * v2[2] - v1[2] * v2[1]; 
	ans[1] = v1[0] * v2[2] - v1[2] * v2[0];  
	ans[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void cross_prod_3_i(int * v1, int * v2, int * ans)
{
	ans[0] = v1[1] * v2[2] - v1[2] * v2[1]; 
	ans[1] = v1[0] * v2[2] - v1[2] * v2[0];  
	ans[2] = v1[0] * v2[1] - v1[1] * v2[0];
}
