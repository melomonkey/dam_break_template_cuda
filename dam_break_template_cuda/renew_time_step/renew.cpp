#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\global_variables.h"
#include "..\update_calculation.h"
#include <ctime>
using namespace std;

void renew(double dt)
{
	{
		generate_neighbors_list();

		void _neighbors_list_IN_gpu();
		_neighbors_list_IN_gpu();

		printf("Generating the Neighbor List...Done"); printf("\n");
	}

	{
		density_update_schm1_cuda(_neigh_list_row, _neigh_list_col, _neighbors_list_cuda, _particles_cuda, dt, 600);
		density_filter_schm1_cuda(_neigh_list_row, _neigh_list_col, _neighbors_list_cuda, _particles_cuda, 600);
		printf("Density Update...Done"); printf("\n");
	}

	{
		pressure_update_schm1_cuda(_neigh_list_row, _neigh_list_col, _neighbors_list_cuda, _particles_cuda, 600);
		printf("Pressure Update...Done"); printf("\n");
	}

	{
		acc_pressure_term_update_schm1_cuda(_neigh_list_row, _neigh_list_col, _neighbors_list_cuda, _particles_cuda, dt, 600);
		printf("Acceleration Update...Done"); printf("\n");
	}

	{
		velocity_update_schm1_cuda(_neigh_list_row, _neigh_list_col, _neighbors_list_cuda, _particles_cuda, dt, 600);
		velocity_filter_schm1_cuda(_neigh_list_row, _neigh_list_col, _neighbors_list_cuda, _particles_cuda, 600);
		printf("Velocity Update...Done"); printf("\n");
	}

	{
		position_time_update_cuda(_neigh_list_row, _neigh_list_col, _neighbors_list_cuda, _particles_cuda, dt, 600);
		printf("Time and Position Update...Done"); printf("\n");
	}

	{
		support_radius_update_schm1_cuda(_neigh_list_row, _neigh_list_col, _neighbors_list_cuda, _particles_cuda, dt, 600);
		printf("Radius Update...Done"); printf("\n");
	}

	void _coordination_list_OUT_gpu();
	void update_coordination_list_gpu();

	update_coordination_list_gpu();
	_coordination_list_OUT_gpu();

	cur_calculation_time += dt;

	printf("%12.11f s Calculation Done.\n\n\n", cur_calculation_time);
}

void renew1(double dt)
{
	{
		generate_neighbors_list();

		void _neighbors_list_IN_gpu(); 
		_neighbors_list_IN_gpu(); 

		printf("Generating the Neighbor List...Done"); printf("\n");
	}

	{
		density_update_schm1_cuda(_neigh_list_row, _neigh_list_col, _neighbors_list_cuda, _particles_cuda, dt, 1200);
		printf("Density Update...Done"); printf("\n");
	}
	
	{
		pressure_update_schm1_cuda(_neigh_list_row, _neigh_list_col, _neighbors_list_cuda, _particles_cuda, 1200);
		printf("Pressure Update...Done"); printf("\n");
	}

	{
		acc_pressure_term_update_schm1_cuda(_neigh_list_row, _neigh_list_col, _neighbors_list_cuda, _particles_cuda, dt, 1200);
		printf("Acceleration Update...Done"); printf("\n");
	}

	{
		velocity_update_schm1_cuda(_neigh_list_row, _neigh_list_col, _neighbors_list_cuda, _particles_cuda, dt, 1200);
		printf("Velocity Update...Done"); printf("\n");
	}

	{
		position_time_update_cuda(_neigh_list_row, _neigh_list_col, _neighbors_list_cuda, _particles_cuda, dt, 1200);
		printf("Time and Position Update...Done"); printf("\n");
	}

	{
		support_radius_update_schm1_cuda(_neigh_list_row, _neigh_list_col, _neighbors_list_cuda, _particles_cuda, dt); 
		printf("Radius Update...Done"); printf("\n");
	}

	void _coordination_list_OUT_gpu();
	void update_coordination_list_gpu();

	update_coordination_list_gpu();
	_coordination_list_OUT_gpu();

	cur_calculation_time += dt; 

	printf("%12.11f s Calculation Done.\n\n\n", cur_calculation_time);
}

//void renew(double dt)
//{
//	{
//		generate_neighbors_list();
//		printf("Generating the Neighbor List...Done"); printf("\n");
//	}
//
//	{
//		density_update_schm1_multi_threads(neighbors_list, particles, dt, 9);
//		density_filter1_multi_threads(neighbors_list, particles, 9);
//		printf("Density Update...Done"); printf("\n");
//	}
//
//	{
//		pressure_update_schm1_multi_threads(neighbors_list, particles, 9);
//		printf("Pressure Update...Done"); printf("\n");
//	}
//
//	{
//		acc_pressure_term_update_schm1_multi_threads(neighbors_list, particles, 9);
//		printf("Acceleration Update...Done"); printf("\n");
//	}
//
//	{
//		velocity_update_schm1_multi_threads(neighbors_list, particles, dt, 9);
//		velocity_filter1_multi_threads(neighbors_list, particles, 9);
//		printf("Velocity Update...Done"); printf("\n");
//	}
//
//	{
//		position_time_update_schm1(particles, dt);
//		printf("Time and Position Update...Done"); printf("\n");
//	}
//
//	{
//		support_radius_update_schm2_multi_threads(neighbors_list, particles, dt);
//		printf("Radius Update...Done"); printf("\n");
//	}
//
//	printf("%12.11f s Calculation Done.\n\n\n", particles[0].time.val[0]);
//}
//
//void renew1(double dt)
//{
//	{
//		generate_neighbors_list();
//		printf("Generating the Neighbor List...Done"); printf("\n");
//	}
//
//	{
//		density_update_schm1_multi_threads(neighbors_list, particles, dt, 9);
//		printf("Density Update...Done"); printf("\n");
//	}
//
//	{
//		pressure_update_schm1_multi_threads(neighbors_list, particles, 9);
//		printf("Pressure Update...Done"); printf("\n");
//	}
//
//	{
//		acc_pressure_term_update_schm1_multi_threads(neighbors_list, particles, 9);
//		printf("Acceleration Update...Done"); printf("\n");
//	}
//
//	{
//		velocity_update_schm1_multi_threads(neighbors_list, particles, dt, 9);
//		printf("Velocity Update...Done"); printf("\n");
//	}
//
//	{
//		position_time_update_schm1(particles, dt);
//		printf("Time and Position Update...Done"); printf("\n");
//	}
//
//	{
//		support_radius_update_schm2_multi_threads(neighbors_list, particles, dt);
//		printf("Radius Update...Done"); printf("\n");
//	}
//
//	printf("%12.11f s Calculation Done.\n\n\n", particles[0].time.val[0]);
//}