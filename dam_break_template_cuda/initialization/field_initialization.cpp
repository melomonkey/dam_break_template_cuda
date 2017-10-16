#include <iostream>
#include <vector>
#include <fstream>
#include "..\struct_particle.h"
#include "..\configuration.h"
#include "..\global_variables.h"
#include "..\update_calculation.h"
#include "..\basic_func\basic_func.h"

using namespace std;

void initial_particles_mass(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles)
{
	for (int ii = 0; ii < 15; ++ii)
	{
		for (int i = 0; i < neighbors_list.size(); ++i) particles[neighbors_list[i][0]].mass.val[1] = 0.0;

		float total_J = 0.0; 

		for (int i = 0; i < neighbors_list.size(); ++i)
		{
			int calcu_particle = neighbors_list[i][0];

			double local_x, local_y, local_z;
			local_x = particles[calcu_particle].coorX.val[0];
			local_y = particles[calcu_particle].coorY.val[0];
			local_z = particles[calcu_particle].coorZ.val[0];
			double h = particles[calcu_particle].smthR.val[0];

			double tmp_density = 0.0;
			for (int j = 1; j < neighbors_list[i].size(); ++j)
			{
				int label_ij = neighbors_list[i][j];

				if (particles[label_ij].id == 0 || particles[label_ij].id == 1 || particles[label_ij].id == 2 || particles[label_ij].id == 3) continue;

				double tmp_x, tmp_y, tmp_z;
				tmp_x = particles[label_ij].coorX.val[0];
				tmp_y = particles[label_ij].coorY.val[0];
				tmp_z = particles[label_ij].coorZ.val[0];

				double dis_x, dis_y, dis_z;
				dis_x = local_x - tmp_x;
				dis_y = local_y - tmp_y;
				dis_z = local_z - tmp_z;

				tmp_density += particles[label_ij].mass.val[0] * kernel_function(dis_x, h);

				double aa = kernel_function(dis_x, h);
				tmp_density += particles[label_ij].mass.val[0] * kernel_function(dis_y, h);
				if (CASE_DIM == 1) tmp_density += particles[label_ij].mass.val[0] * kernel_function(dis_z, h);
			}

			double erro = tmp_density - particles[calcu_particle].density.val[0];

			total_J += erro * erro; 

			for (int j = 1; j < neighbors_list[i].size(); ++j)
			{
				int label_ij = neighbors_list[i][j];

				if (particles[label_ij].id == 0 || particles[label_ij].id == 1 || particles[label_ij].id == 2 || particles[label_ij].id == 3) continue;

				double tmp_x, tmp_y, tmp_z;
				tmp_x = particles[label_ij].coorX.val[0];
				tmp_y = particles[label_ij].coorY.val[0];
				tmp_z = particles[label_ij].coorZ.val[0];

				double dis_x, dis_y, dis_z;
				dis_x = local_x - tmp_x;
				dis_y = local_y - tmp_y;
				dis_z = local_z - tmp_z;

				double rlx_factor = 1e-3 * (particles[calcu_particle].mass.val[0] / erro) * particles[calcu_particle].mass.val[0] / particles[calcu_particle].density.val[0];

				particles[label_ij].mass.val[1] += -erro * rlx_factor * kernel_function(dis_x, h);
				particles[label_ij].mass.val[1] += -erro * rlx_factor * kernel_function(dis_y, h);
				if (CASE_DIM == 1) particles[label_ij].mass.val[1] += -erro * rlx_factor * kernel_function(dis_z, h);
			}
		}

		total_J /= neighbors_list.size();

		for (int i = 0; i < neighbors_list.size(); ++i) particles[neighbors_list[i][0]].mass.val[0] += particles[neighbors_list[i][0]].mass.val[1];
	}
	
}



double investigate_interpolation_mass()
{
	double obv_box[6] = { 0.4, 0.6, 0.4, 0.6, 0.4, 0.6 };

	double tmp_sum_obv_val = 0.0;
	int obv_num = 0;
	for (int i = 0; i < particles.size(); ++i)
	{
		if (particles[i].coorX.val[0] >= obv_box[0] && particles[i].coorX.val[0] <= obv_box[1])
		{
			if (particles[i].coorY.val[0] >= obv_box[2] && particles[i].coorY.val[0] <= obv_box[3])
			{
				if (particles[i].coorZ.val[0] >= obv_box[4] && particles[i].coorZ.val[0] <= obv_box[5])
				{
					tmp_sum_obv_val += particles[i].density.val[0];
					obv_num++;
				}
			}
		}
	}

	return tmp_sum_obv_val /= obv_num;
}

void drop3d_initialization()
{
	void generate_neighbors_list();
	generate_neighbors_list();

	while (1)
	{
		density_update_schm2(neighbors_list, particles);
		double abbre = investigate_interpolation_mass();

		if (abs(abbre - REF_DENSITY) < 1e-3) break;

		double mass = particles[0].mass.val[0];
		if (abbre - REF_DENSITY > 0) mass -= 0.01 * mass;
		if (abbre - REF_DENSITY < 0) mass += 0.01 * mass;

		for (int i = 0; i < particles.size(); ++i)
			particles[i].mass.val[0] = mass;
	}
}

void drop_initialization()
{
}

void heat_transfer_initialization()
{
	for (int i = 0; i < particles.size(); ++i)
	{
		if (particles[i].coorZ.val[0] >= 0.1)
		{
			particles[i].temperature.val[0] = 100.0;
		}

		if (abs(particles[i].coorX.val[0]) >= 0.1 - 0.0045 || abs(particles[i].coorY.val[0]) >= 0.1 - 0.0045 || abs(particles[i].coorZ.val[0]) >= 0.1 - 0.0045) particles[i].id = 0;
		if (abs(particles[i].coorX.val[0]) <= 0.0045 || abs(particles[i].coorY.val[0]) <= 0.0045 || abs(particles[i].coorZ.val[0]) <= 0.0045) particles[i].id = 0;
	}
}

double investigation_balanced_mass()
{
	double obv_box[6] = { 0.2, 0.3, 0.3, 0.4, 0.1, 0.2 };

	double tmp_sum_obv_val = 0.0;
	int obv_num = 0;
	for (int i = 0; i < particles.size(); ++i)
	{
		if (particles[i].coorX.val[0] >= obv_box[0] && particles[i].coorX.val[0] <= obv_box[1])
		{
			if (particles[i].coorY.val[0] >= obv_box[2] && particles[i].coorY.val[0] <= obv_box[3])
			{
				if (particles[i].coorZ.val[0] >= obv_box[4] && particles[i].coorZ.val[0] <= obv_box[5])
				{
					tmp_sum_obv_val += particles[i].accZ.val[0];
					obv_num++;
				}
			}
		}
	}

	return tmp_sum_obv_val /= obv_num;
}

double investigation_balanced_mass_2D()
{
	double obv_box[4] = { 0.2, 0.4, 0.1, 0.15 };

	double tmp_sum_obv_val = 0.0;
	int obv_num = 0;
	for (int i = 0; i < particles.size(); ++i)
	{
		if (particles[i].coorX.val[0] >= obv_box[0] && particles[i].coorX.val[0] <= obv_box[1])
		{
			if (particles[i].coorY.val[0] >= obv_box[2] && particles[i].coorY.val[0] <= obv_box[3])
			{
				tmp_sum_obv_val += particles[i].accY.val[0];
				obv_num++;
			}
		}
	}

	return tmp_sum_obv_val /= obv_num;
}

void dambreak3d_field_initialization()
{
	for (int i = 0; i < particles.size(); ++i)
	{
		if (particles[i].id == 4)
		{
			particles[i].pressure.val[0] = REF_DENSITY * GRAVITY * (0.3 - particles[i].coorY.val[0]);
			particles[i].pressure.val[1] = particles[i].pressure.val[0];
			particles[i].pressure.val[2] = particles[i].pressure.val[0];
		}
	}

	for (int i = 0; i < particles.size(); ++i)
	{
		if (particles[i].id == 4)
		{
			particles[i].density.val[0] = REF_DENSITY * pow((particles[i].pressure.val[0] / COMPRS_LIM_CONST_B) + 1.0, (1.0 / 7.0));
			particles[i].density.val[1] = particles[i].density.val[0];
			particles[i].density.val[2] = particles[i].density.val[0];
		}
	}

	//support_radius_update_schm1(particles);

	/*void generate_neighbors_list();
	generate_neighbors_list();

	while(1)
	{
		acc_pressure_term_update_schm1(neighbors_list, &particles);
		double abbre = investigation_balanced_mass();

		if(abs(abbre) < 1e-5) break;

		double mass = particles[0].mass.val[0];
		if(abbre > 0) mass -= 0.01 * mass;
		if(abbre < 0) mass += 0.01 * mass;

		for(int i = 0; i < particles.size(); ++ i)
			particles[i].mass.val[0] = mass;
	}*/
}

void dambreak2d_field_initialization()
{
	for (int i = 0; i < particles.size(); ++i)
	{
		if (particles[i].id == 4)
		{
			particles[i].pressure.val[0] = REF_DENSITY * GRAVITY * (0.25 - particles[i].coorY.val[0]);
			particles[i].pressure.val[1] = particles[i].pressure.val[0];
			particles[i].pressure.val[2] = particles[i].pressure.val[0];
		}
	}

	for (int i = 0; i < particles.size(); ++i)
	{
		if (particles[i].id == 4)
		{
			particles[i].density.val[0] = REF_DENSITY * pow((particles[i].pressure.val[0] / COMPRS_LIM_CONST_B) + 1.0, (1.0 / 7.0));
			particles[i].density.val[1] = particles[i].density.val[0];
			particles[i].density.val[2] = particles[i].density.val[0];
		}
	}

	/*void generate_neighbors_list();
	generate_neighbors_list();

	while(1)
	{
		acc_pressure_term_update_schm1(neighbors_list, &particles);
		double abbre = investigation_balanced_mass_2D();

		if(abs(abbre) < 1e-5) break;

		double mass = particles[0].mass.val[0];
		if(abbre > 0) mass -= 0.01 * mass;
		if(abbre < 0) mass += 0.01 * mass;

		for(int i = 0; i < particles.size(); ++ i)
			particles[i].mass.val[0] = mass;
	}*/
}

