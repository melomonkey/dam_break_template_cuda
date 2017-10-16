#include "output.h"

void output_wall()
{
	string file_name = "OUTPUT_RESULTS\\PLOT_TECPLOT.plt";

	fstream fout(file_name, ios::out | ios::app);

	int count = 0; 
	for(int i = 0; i < particles.size(); ++ i)
	{
		if(particles[i].id == 0) count++; 
	}

	if(count == 0) return; 

	fout <<"VARIABLES = \"X\" \"Y\" \"Z\" \"VEL_X\" \"VEL_Y\" \"VEL_Z\" \"ACC_X\" \"ACC_Y\" \"ACC_Z\" \"PRESSURE\" \"DENSITY\" \"VISCO\" \"RADIUS\" \"TYPE\""<< endl;
	fout <<"ZONE I = " << count << ", F = POINT" << ", SolutionTime = "<< particles[0].time.val[0] <<endl; 

	for(int i = 0; i < particles.size(); ++ i)
	{
		if(particles[i].id == 0)
		{
			fout << particles[i].coorX.val[0] << " " << particles[i].coorY.val[0] << " " << particles[i].coorZ.val[0] << " "
				 << particles[i].velX.val[0] << " " << particles[i].velY.val[0] << " " << particles[i].velZ.val[0] << " " 
				 << particles[i].normal_vec[0] << " " << particles[i].normal_vec[1]<< " " << particles[i].normal_vec[2] << " " 
				 << particles[i].pressure.val[0] << " " << particles[i].density.val[0] << " " << particles[i].viscosity.val[0] << " " << particles[i].smthR.val[0] << " "
				 << particles[i].id << " "
			     << endl; 
		}
	}
}


void output_res_bin()
{
	ofstream fout("OUTPUT_RESULTS\\result_data_bin.bin", std::ios::binary | ios::app);

	int output_particles_num = 0;
	for(int i = 0; i < particles.size(); ++ i) output_particles_num ++; 

	int var_num = 18; 

	fout.write((char *)&output_particles_num, sizeof(int)); 
	fout.write((char *)&var_num, sizeof(int)); 
	
	for(int i = 0; i < particles.size(); ++ i)
	{
		double tmp_id = (double)particles[i].id; 
		fout.write((char *)&tmp_id, sizeof(double)); // #0

		fout.write((char *)&particles[i].mass.val[0], sizeof(double)); // #1
		fout.write((char *)&particles[i].density.val[0], sizeof(double)); // #2
		fout.write((char *)&particles[i].viscosity.val[0], sizeof(double)); // #3
		fout.write((char *)&particles[i].pressure.val[0], sizeof(double)); // #4
		fout.write((char *)&particles[i].temperature.val[0], sizeof(double)); // #5
		fout.write((char *)&particles[i].energy.val[0], sizeof(double)); // #6
		fout.write((char *)&particles[i].velX.val[0], sizeof(double)); // #7
		fout.write((char *)&particles[i].velY.val[0], sizeof(double)); // #8
		fout.write((char *)&particles[i].velZ.val[0], sizeof(double)); // #9
		fout.write((char *)&particles[i].accX.val[0], sizeof(double)); // #10
		fout.write((char *)&particles[i].accY.val[0], sizeof(double)); // #11
		fout.write((char *)&particles[i].accZ.val[0], sizeof(double)); // #12
		fout.write((char *)&particles[i].coorX.val[0], sizeof(double)); // #13
		fout.write((char *)&particles[i].coorY.val[0], sizeof(double)); // #14
		fout.write((char *)&particles[i].coorZ.val[0], sizeof(double)); // #15
		fout.write((char *)&particles[i].smthR.val[0], sizeof(double)); // #16
		fout.write((char *)&particles[i].time.val[0], sizeof(double)); // #17
	}

	fout.close(); 
}

void output_saved()
{
	string file_name = "OUTPUT_RESULTS\\_SAVED.txt";

	fstream fout(file_name, ios::out);

	fout << particles.size() << endl; 
	for(int i = 0; i < particles.size(); ++ i)
	{
		fout << particles[i].id << " " <<
			    particles[i].mass.val[0] << " " << particles[i].density.val[0] << " " << particles[i].viscosity.val[0] << " " << particles[i].pressure.val[0] << " " << 
				particles[i].temperature.val[0] << " " << particles[i].energy.val[0] << " " <<
			    particles[i].velX.val[0] << " " << particles[i].velY.val[0] << " " << particles[i].velZ.val[0] << " " <<
				particles[i].accX.val[0] << " " << particles[i].accY.val[0] << " " << particles[i].accZ.val[0] << " " <<
				particles[i].coorX.val[0] << " " << particles[i].coorY.val[0] << " " << particles[i].coorZ.val[0] << " " <<
				particles[i].smthR.val[0] << " " << particles[i].time.val[0] << " " << endl;				  
	}
}



void output_tecplot()
{
	string file_name = "OUTPUT_RESULTS\\PLOT_TECPLOT.plt";

	fstream fout(file_name, ios::out | ios::app);

	int count = 0; 
	for(int i = 0; i < particles.size(); ++ i)
	{
		if(particles[i].id == 4) count++; 
	}

	if(count == 0) return; 

	fout <<"VARIABLES = \"X\" \"Y\" \"Z\" \"VEL_X\" \"VEL_Y\" \"VEL_Z\" \"ACC_X\" \"ACC_Y\" \"ACC_Z\" \"PRESSURE\" \"DENSITY\" \"VISCO\" \"RADIUS\" \"TYPE\""<< endl;
	fout <<"ZONE I = " << count << ", F = POINT" << ", SolutionTime = "<< particles[0].time.val[0] <<endl; 

	for(int i = 0; i < particles.size(); ++ i)
	{
		if(particles[i].id == 4)
		{
			fout << particles[i].coorX.val[0] << " " << particles[i].coorY.val[0] << " " << particles[i].coorZ.val[0] << " "
				 << particles[i].velX.val[0] << " " << particles[i].velY.val[0] << " " << particles[i].velZ.val[0] << " " 
				 << particles[i].accX.val[0] << " " << particles[i].accY.val[0] << " " << particles[i].accZ.val[0] << " " 
				 << particles[i].pressure.val[0] << " " << particles[i].density.val[0] << " " << particles[i].viscosity.val[0] << " " << particles[i].smthR.val[0] << " "
				 << particles[i].id << " "
			     << endl; 
		}
	}

	count = 0; 
	for(int i = 0; i < particles.size(); ++ i)
	{
		if(particles[i].id == 0) count++; 
	}

	if(count == 0) return; 

	fout <<"VARIABLES = \"X\" \"Y\" \"Z\" \"VEL_X\" \"VEL_Y\" \"VEL_Z\" \"ACC_X\" \"ACC_Y\" \"ACC_Z\" \"PRESSURE\" \"DENSITY\" \"VISCO\" \"RADIUS\" \"TYPE\""<< endl;
	fout <<"ZONE I = " << count << ", F = POINT" << ", SolutionTime = "<< particles[0].time.val[0] <<endl; 

	for(int i = 0; i < particles.size(); ++ i)
	{
		if(particles[i].id == 0)
		{
			fout << particles[i].coorX.val[0] << " " << particles[i].coorY.val[0] << " " << particles[i].coorZ.val[0] << " "
				 << particles[i].velX.val[0] << " " << particles[i].velY.val[0] << " " << particles[i].velZ.val[0] << " " 
				 << particles[i].accX.val[0] << " " << particles[i].accY.val[0] << " " << particles[i].accZ.val[0] << " " 
				 << particles[i].pressure.val[0] << " " << particles[i].density.val[0] << " " << particles[i].viscosity.val[0] << " " << particles[i].smthR.val[0] << " "
				 << particles[i].id << " "
			     << endl; 
		}
	}
}