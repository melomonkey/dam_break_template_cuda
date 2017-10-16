#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\global_variables.h"
#include "..\update_calculation.h"
#include <ctime>
using namespace std;

void renew_leapfrog(double dt)
{
	{
		density_update_leapfrog_schm1a(neighbors_list, &particles, .5 * dt); 
		printf("Density Update...Done"); printf("\n");
	}

	{
		pressure_update_schm1(neighbors_list, &particles); 
		printf("Pressure Update...Done"); printf("\n");
	}

	{
		acc_pressure_term_update_schm1(neighbors_list, &particles);
		acc_boundaries_update_schm1(neighbors_list, &particles); 
	    printf("Acceleration Update...Done"); printf("\n");
	}

	{
		velocity_update_leap_frog_schm1a(neighbors_list, &particles, .5 * dt); 
		printf("Velocity Update...Done"); printf("\n");
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * */

	{
		position_time_update_schm1(&particles, dt); 
	    printf("Time and Position Update...Done"); printf("\n");
	}

	{
		void generate_neighbors_list(); 

		generate_neighbors_list(); 
		
	    printf("Generating the Neighbor List...Done"); printf("\n"); 
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * */

	{
		density_update_leapfrog_schm1b(neighbors_list, &particles, .5 * dt); 
		printf("Density Update...Done"); printf("\n");
	}

	{
		pressure_update_schm1(neighbors_list, &particles); 
		printf("Pressure Update...Done"); printf("\n");
	}

	{
		velocity_update_leap_frog_schm1b(neighbors_list, &particles, .5 * dt); 
		printf("Velocity Update...Done"); printf("\n");
	}

	printf("%12.11f s Calculation Done.\n\n\n", particles[0].time.val[0]); 
}

void renew_leapfrog1(double dt)
{
	{
		density_update_leapfrog_schm1a(neighbors_list, &particles, .5 * dt); 
		density_filter1(neighbors_list, &particles); 
		printf("Density Update...Done"); printf("\n");
	}

	{
		pressure_update_schm1(neighbors_list, &particles); 
		printf("Pressure Update...Done"); printf("\n");
	}

	{
		acc_pressure_term_update_schm1(neighbors_list, &particles);
		acc_boundaries_update_schm1(neighbors_list, &particles); 
	    printf("Acceleration Update...Done"); printf("\n");
	}

	{
		velocity_update_leap_frog_schm1a(neighbors_list, &particles, .5 * dt); 
		velocity_filter3(neighbors_list, &particles);
		printf("Velocity Update...Done"); printf("\n");
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * */

	{
		position_time_update_schm1(&particles, dt); 
	    printf("Time and Position Update...Done"); printf("\n");
	}

	{
		void generate_neighbors_list(); 

		generate_neighbors_list(); 
		
	    printf("Generating the Neighbor List...Done"); printf("\n"); 
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * */

	{
		density_update_leapfrog_schm1b(neighbors_list, &particles, .5 * dt); 
		density_filter1(neighbors_list, &particles); 
		printf("Density Update...Done"); printf("\n");
	}

	{
		pressure_update_schm1(neighbors_list, &particles); 
		printf("Pressure Update...Done"); printf("\n");
	}

	{
		velocity_update_leap_frog_schm1b(neighbors_list, &particles, .5 * dt);
		velocity_filter3(neighbors_list, &particles); 
		printf("Velocity Update...Done"); printf("\n");
	}

	printf("%12.11f s Calculation Done.\n\n\n", particles[0].time.val[0]); 
}