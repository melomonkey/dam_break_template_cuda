#include "input.h"
/* Initial input variables do not contain energy */

void input_initial_water()
{
	fstream fin; 
	fin.open("INPUT_MODELS\\test2d.txt", ios::in);

	vector<string> inf; 

	while(!fin.eof())
	{
		string tmp_line; 
	    getline(fin, tmp_line); 

		inf.push_back(tmp_line); 
	}

	inf.pop_back();

	vector<vector<double>> coordination; 
	coordination.resize(inf.size()); 
	for(int i = 0; i < inf.size(); ++ i)
	{
		vector<double> tmp; 
		tmp.resize(3);

		float a, b, c; 	  
		sscanf_s(inf[i].c_str(), "%f %f %f", &a, &b, &c);

		tmp[0] = a; tmp[1] = b; tmp[2] = c; 

		coordination[i] = tmp; 
	}
	
	for(int i = 0; i < coordination.size(); ++ i) 
	{
		PARTICLE tmp; 

		double scale = 1.0; 

		tmp.id = 4; 

		tmp.coorX.val[0] = coordination[i][0] * scale; tmp.coorX.val[1] = tmp.coorX.val[0]; tmp.coorX.val[2] = tmp.coorX.val[0]; 
		tmp.coorY.val[0] = coordination[i][1] * scale; tmp.coorY.val[1] = tmp.coorY.val[0]; tmp.coorY.val[2] = tmp.coorY.val[0];
		tmp.coorZ.val[0] = coordination[i][2] * scale; tmp.coorZ.val[1] = tmp.coorZ.val[0]; tmp.coorZ.val[2] = tmp.coorZ.val[0];

		tmp.mass.val[0] = 0.0082; tmp.mass.val[1] = tmp.mass.val[0]; tmp.mass.val[2] = tmp.mass.val[0];
		tmp.density.val[0] = REF_DENSITY; tmp.density.val[1] = tmp.density.val[0]; tmp.density.val[2] = tmp.density.val[0];
		tmp.viscosity.val[0] = 1e-6; tmp.viscosity.val[1] = tmp.viscosity.val[0]; tmp.viscosity.val[2] = tmp.viscosity.val[0];

		tmp.pressure.val[0] = 0.0; tmp.pressure.val[1] = tmp.pressure.val[0]; tmp.pressure.val[2] = tmp.pressure.val[0];
		tmp.temperature.val[0] = 0.0; tmp.temperature.val[1] = tmp.temperature.val[0]; tmp.temperature.val[2] = tmp.temperature.val[0];
		tmp.energy.val[0] = 0.0; tmp.energy.val[1] = tmp.energy.val[0]; tmp.energy.val[2] = tmp.energy.val[0];

		tmp.accX.val[0] = 0.0; tmp.accX.val[1] = tmp.accX.val[0]; tmp.accX.val[2] = tmp.accX.val[0];
		tmp.accY.val[0] = 0.0; tmp.accY.val[1] = tmp.accY.val[0]; tmp.accY.val[2] = tmp.accY.val[0];
		tmp.accZ.val[0] = 0.0; tmp.accZ.val[1] = tmp.accZ.val[0]; tmp.accZ.val[2] = tmp.accZ.val[0];

		tmp.velX.val[0] = 0.0; tmp.velX.val[1] = tmp.velX.val[0]; tmp.velX.val[2] = tmp.velX.val[0];
		tmp.velY.val[0] = 0.0; tmp.velY.val[1] = tmp.velY.val[0]; tmp.velY.val[2] = tmp.velY.val[0];
		tmp.velZ.val[0] = 0.0; tmp.velZ.val[1] = tmp.velZ.val[0]; tmp.velZ.val[2] = tmp.velZ.val[0];

		tmp.smthR.val[0] = 3.333e-3 * 1.5; tmp.smthR.val[1] = tmp.smthR.val[0]; tmp.smthR.val[2] = tmp.smthR.val[0];

		tmp.time.val[0] = 0.0; tmp.time.val[1] = tmp.time.val[0]; tmp.time.val[2] = tmp.time.val[0];

		particles.push_back(tmp); 
	}

	fin.close(); 
}

//void input_initial_box()
//{
//	fstream fin; 
//	fin.open("INPUT_MODELS\\box0.02.txt", ios::in);
//
//	vector<string> inf; 
//
//	while(!fin.eof())
//	{
//		string tmp_line; 
//	    getline(fin, tmp_line); 
//
//		inf.push_back(tmp_line); 
//	}
//
//	inf.pop_back();
//
//	vector<vector<double>> coordination; 
//	coordination.resize(inf.size()); 
//	for(int i = 0; i < inf.size(); ++ i)
//	{
//		vector<double> tmp; 
//		tmp.resize(3);
//
//		float a, b, c; 	  
//		sscanf_s(inf[i].c_str(), "%f %f %f", &a, &b, &c);
//
//		tmp[0] = a; tmp[1] = b; tmp[2] = c; 
//
//		coordination[i] = tmp; 
//	}
//	
//	for(int i = 0; i < coordination.size(); ++ i) 
//	{
//		PARTICLE tmp; 
//
//		double scale = 1.0; 
//
//		tmp.id = 0; 
//
//		tmp.coorX.val[0] = coordination[i][0] * scale; 
//		tmp.coorY.val[0] = coordination[i][1] * scale; 
//		tmp.coorZ.val[0] = coordination[i][2] * scale; 
//
//		tmp.mass.val[0] = 0.0; tmp.mass.val[1] = tmp.mass.val[0]; tmp.mass.val[2] = tmp.mass.val[0];
//		tmp.density.val[0] = REF_DENSITY; tmp.density.val[1] = tmp.density.val[0]; tmp.density.val[2] = tmp.density.val[0];
//		tmp.viscosity.val[0] = 0.0; tmp.viscosity.val[1] = tmp.viscosity.val[0]; tmp.viscosity.val[2] = tmp.viscosity.val[0];
//
//		tmp.pressure.val[0] = 0.0; tmp.pressure.val[1] = tmp.pressure.val[0]; tmp.pressure.val[2] = tmp.pressure.val[0];
//		tmp.temperature.val[0] = 0.0; tmp.temperature.val[1] = tmp.temperature.val[0]; tmp.temperature.val[2] = tmp.temperature.val[0];
//		tmp.energy.val[0] = 0.0; tmp.energy.val[1] = tmp.energy.val[0]; tmp.energy.val[2] = tmp.energy.val[0];
//
//		tmp.accX.val[0] = 0.0; tmp.accX.val[1] = tmp.accX.val[0]; tmp.accX.val[2] = tmp.accX.val[0];
//		tmp.accY.val[0] = 0.0; tmp.accY.val[1] = tmp.accY.val[0]; tmp.accY.val[2] = tmp.accY.val[0];
//		tmp.accZ.val[0] = 0.0; tmp.accZ.val[1] = tmp.accZ.val[0]; tmp.accZ.val[2] = tmp.accZ.val[0];
//
//		tmp.velX.val[0] = 0.0; tmp.velX.val[1] = tmp.velX.val[0]; tmp.velX.val[2] = tmp.velX.val[0];
//		tmp.velY.val[0] = 0.0; tmp.velY.val[1] = tmp.velY.val[0]; tmp.velY.val[2] = tmp.velY.val[0];
//		tmp.velZ.val[0] = 0.0; tmp.velZ.val[1] = tmp.velZ.val[0]; tmp.velZ.val[2] = tmp.velZ.val[0];
//
//		tmp.smthR.val[0] = 0.02 * 1.5; tmp.smthR.val[1] = tmp.smthR.val[0]; tmp.smthR.val[2] = tmp.smthR.val[0];
//
//		tmp.normal_vec[0] = 0.0; tmp.normal_vec[1] = 0.0; tmp.normal_vec[2] = 0.0;
//
//		particles.push_back(tmp); 
//	}
//
//	fin.close(); 
//}

void input_saved()
{
	fstream fin("OUTPUT_RESULTS\\_SAVED.txt", ios::in);

	int tmp_num; 
	fin >> tmp_num; 

	for(int i = 0; i < tmp_num; ++ i)
	{
		PARTICLE tmp;
		fin >> tmp.id;
		fin >> tmp.mass.val[0]; fin >> tmp.density.val[0]; fin >> tmp.viscosity.val[0]; fin >> tmp.pressure.val[0]; 
		fin >> tmp.temperature.val[0]; fin >> tmp.energy.val[0];
		fin >> tmp.velX.val[0]; fin >> tmp.velY.val[0]; fin >> tmp.velZ.val[0];
		fin >> tmp.accX.val[0]; fin >> tmp.accY.val[0]; fin >> tmp.accZ.val[0]; 
		fin >> tmp.coorX.val[0]; fin >> tmp.coorY.val[0]; fin >> tmp.coorZ.val[0]; 
		fin >> tmp.smthR.val[0]; fin >> tmp.time.val[0]; 
		 
		particles.push_back(tmp); 
	}
}