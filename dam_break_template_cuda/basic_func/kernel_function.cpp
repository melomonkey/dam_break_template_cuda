#include <math.h>
#include "..\configuration.h"

#define PI 3.1415927

/* * * * * * * * Gaussian Kernel * * * * * * * * * * */
double kernel_function(double dis, double h)
{
	double R = dis / h; 
	double alpha;

	if(CASE_DIM == 0) alpha = (1.0 / (h * h)) * pow(PI, -1.0);
	if(CASE_DIM == 1) alpha = (1.0 / (h * h * h)) * pow(PI, -1.5);

	double ans = alpha * exp(-(R * R));  
	return ans; 
}

double kernel_function_1dev(double dis, double h)
{
	double R = dis / h;
	double alpha;

	if(CASE_DIM == 0) alpha = (1.0 / (h * h)) * pow(PI, -1.0);
	if(CASE_DIM == 1) alpha = (1.0 / (h * h * h)) * pow(PI, -1.5);

	double ans = alpha * exp(-(R * R)) * (-2.0 * R);
	
	return ans / h;
}

double kernel_function_2dev(double dis, double h)
{
	double R = dis / h;
	double alpha;

	if(CASE_DIM == 0) alpha = (1.0 / (h * h)) * pow(PI, -1.0);
	if(CASE_DIM == 1) alpha = (1.0 / (h * h * h)) * pow(PI, -1.5);

	double ans = -2.0 * alpha * (1 - 2 * R * R) * exp(-(R * R));

	return ans / h / h; // local particle's position coordination's derivative, so include the R's derivative, twice.
}





/* * * * * * * * 2D Quadratic * * * * * * * * * * */
//double kernel_function(double dis, double h)
//{
//	double R = dis / h; 
//	double alpha = 2.0 / (PI * pow(h, 2.0)); 
//
//	double sign; 
//	sign = (R < 0.0) ? -1.0 : 1.0; 
//
//	R = abs(R); 
//
//	double ans = alpha * ((3.0 / 16.0) * R * R - (3.0 / 4.0) * R + (3.0 / 4));
//
//	return ans; 
//}

//double kernel_function_1dev(double dis, double h)
//{
//	double R = dis / h;
//	double alpha = 2.0 / (PI * pow(h, 2.0));
//
//	double sign; 
//	sign = (R < 0.0) ? -1.0 : 1.0; 
//
//	R = abs(R); 
//
//	double ans = alpha * ((3.0 / 8.0) * R - (3.0 / 4.0));
//
//	ans *= sign; 
//
//	return ans;
//}

/* * * * * * * * 3D Quadratic * * * * * * * * * * */
//double kernel_function(double dis, double h)
//{
//	double R = dis / h; 
//	double alpha = 5.0 / (4.0 * PI * pow(h, 3.0)); 
//
//	double sign; 
//	sign = (R < 0.0) ? -1.0 : 1.0; 
//
//	R = abs(R);
//
//	double ans = alpha * ((3.0 / 16.0) * R * R - (3.0 / 4.0) * R + (3.0 / 4));
//
//	return ans; 
//}

//double kernel_function_1dev(double dis, double h)
//{
//	double R = dis / h;
//	double alpha = 5.0 / (4.0 * PI * pow(h, 3.0));
//
//	double sign; 
//	sign = (R < 0.0) ? -1.0 : 1.0; 
//
//	R = abs(R);
//
//	double ans = alpha * ((3.0 / 8.0) * R - (3.0 / 4.0));
//
//	ans *= sign; 
//
//	return ans;
//}



/* * * * * * * * 2D Cubic Spline * * * * * * * * * * */
//double kernel_function(double dis, double h)
//{
//	double R = dis / h; 
//	double alpha = 10.0 / (7.0 * PI * pow(h, 2.0)); 
//
//	double sign; 
//	sign = (R < 0.0) ? -1.0 : 1.0; 
//
//	R = abs(R);
//
//	double ans = 0.0; 
//
//	if(R >= 0.0 && R <= 1.0) ans = 1 - (3.0 / 2.0) * pow(R, 2.0) + (3.0 / 4.0) * pow(R, 3.0); 
//
//	if(R >= 1 && R <= 2) ans = (1.0 / 4.0) * pow(2.0 - R, 3.0); 
//
//	ans *= alpha; 
//
//	return ans; 
//}

//double kernel_function_1dev(double dis, double h)
//{
//	double R = dis / h;
//	double alpha = 10.0 / (7.0 * PI * pow(h, 2.0));
//
//	double sign; 
//	sign = (R < 0.0) ? -1.0 : 1.0; 
//
//	R = abs(R);
//
//	double ans = 0.0; 
//
//	if(R >= 0.0 && R <= 1.0) ans = 3.0 * R + (9.0 / 4.0) * pow(R, 2.0); 
//
//	if(R >= 1 && R <= 2) ans = (-3.0 / 4.0) * pow(2.0 - R, 2.0); 
//
//	ans *= alpha; 
//
//	ans *= sign; 
//
//	return ans; 
//}

/* * * * * * * * 3D Cubic Spline * * * * * * * * * * */
//double kernel_function(double dis, double h)
//{
//	double R = dis / h; 
//	double alpha = 1.0 / (PI * pow(h, 3.0)); 
//
//	double sign; 
//	sign = (R < 0.0) ? -1.0 : 1.0; 
//
//	R = abs(R);
//
//	double ans = 0.0; 
//
//	if(R >= 0.0 && R <= 1.0) ans = 1 - (3.0 / 2.0) * pow(R, 2.0) + (3.0 / 4.0) * pow(R, 3.0); 
//
//	if(R >= 1 && R <= 2) ans = (1.0 / 4.0) * pow(2.0 - R, 3.0); 
//
//	ans *= alpha; 
//
//	return ans;  
//}

//double kernel_function_1dev(double dis, double h)
//{
//	double R = dis / h;
//	double alpha = 1.0 / (PI * pow(h, 3.0));
//
//	double sign; 
//	sign = (R < 0.0) ? -1.0 : 1.0; 
//
//	R = abs(R);
//
//	double ans = 0.0; 
//
//	if(R >= 0.0 && R <= 1.0) ans = 3.0 * R + (9.0 / 4.0) * pow(R, 2.0); 
//
//	if(R >= 1 && R <= 2) ans = (-3.0 / 4.0) * pow(2.0 - R, 2.0); 
//
//	ans *= alpha; 
//
//	ans *= sign; 
//
//	return ans; 
//}



/* * * * * * * * 2D Lucy * * * * * * * * * * */
//double kernel_function(double dis, double h)
//{
//	double R = dis / h; 
//	double alpha = 7.0 / (4.0 * PI * pow(h, 2.0)); 
//
//	double sign; 
//	sign = (R < 0.0) ? -1.0 : 1.0; 
//
//	R = abs(R);
//
//	double ans = 0.0; 
//
//	if(R >= 0.0 && R <= 1.0) ans = (1 + 3 * R) * pow((1 - R), 3.0); 
//
//	ans *= alpha; 
//
//	return ans; 
//}

//double kernel_function_1dev(double dis, double h)
//{
//	double R = dis / h;
//	double alpha = 10.0 / (7.0 * PI * pow(h, 2.0));
//
//	double sign; 
//	sign = (R < 0.0) ? -1.0 : 1.0; 
//
//	R = abs(R);
//
//	double ans = 0.0; 
//
//	if(R >= 0.0 && R <= 1.0) ans = 3.0 * pow((1 - R), 3.0) + (1 + 3 * R) * 3.0 * (-1.0) * pow((1 - R), 2.0); 
//
//	ans *= alpha; 
//
//	ans *= sign; 
//
//	return ans; 
//}

/* * * * * * * * 3D Lucy * * * * * * * * * * */
//double kernel_function(double dis, double h)
//{
//	double R = dis / h; 
//	double alpha = 1.0 / (PI * pow(h, 3.0)); 
//
//	double sign; 
//	sign = (R < 0.0) ? -1.0 : 1.0; 
//
//	R = abs(R);
//
//	double ans = 0.0; 
//
//	if(R >= 0.0 && R <= 1.0) if(R >= 0.0 && R <= 1.0) ans = (1 + 3 * R) * pow((1 - R), 3.0);
//
//	ans *= alpha; 
//
//	return ans;  
//}

//double kernel_function_1dev(double dis, double h)
//{
//	double R = dis / h;
//	double alpha = 1.0 / (PI * pow(h, 3.0));
//
//	double sign; 
//	sign = (R < 0.0) ? -1.0 : 1.0; 
//
//	R = abs(R);
//
//	double ans = 0.0; 
//
//	if(R >= 0.0 && R <= 1.0) ans = 3.0 * pow((1 - R), 3.0) + (1 + 3 * R) * 3.0 * (-1.0) * pow((1 - R), 2.0); 
//
//	ans *= alpha;
//
//	ans *= sign; 
//
//	return ans; 
//}
