#ifndef update_calculation_h_
#define update_calculation_h_

#include "struct_particle.h"
#include <vector>
using namespace std;

void generate_neighbors_list(); 

/* * * * * * * * * * * * * * * * * * * * * */

void density_filter1(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles);
void density_filter1_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, int threads_num = 5);
void density_filter2(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles);
void density_filter2_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, int threads_num = 5);

void density_filter_schm1_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * _neigh_list_cuda, PARTICLE * _particles_cuda, int blcks = 50); 

void density_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, double dt);
void density_update_schm1_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, double dt, int threads_num = 5);
void density_update_schm2(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles);
void density_update_schm2_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, int threads_num = 5);

void density_update_schm1_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * _neigh_list_cuda, PARTICLE * _particles_cuda, double dt, int blcks = 50);

void density_update_leapfrog_schm1a(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles, double dt);
void density_update_leapfrog_schm1b(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles, double dt);
/* * * * * * * * * * * * * * * * * * * * * */

void pressure_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles);
void pressure_update_schm1_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, int threads_num = 5);
void pressure_update_schm2(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles);

void pressure_update_schm1_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * _neigh_list_cuda, PARTICLE * _particles_cuda, int blcks = 50);

/* * * * * * * * * * * * * * * * * * * * * */

void temperature_filter1(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles);
void temperature_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles, double dt);
void temperature_update_schm2(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles, double dt);

/* * * * * * * * * * * * * * * * * * * * * */

void energy_pressure_term_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles, double dt);
void energy_pressure_term_update_schm2(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles, double dt);

/* * * * * * * * * * * * * * * * * * * * * */
void acc_pressure_term_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles);
void acc_pressure_term_update_schm1_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, int threads_num = 5);
void acc_pressure_term_update_schm2(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, double dt);
void acc_pressure_term_update_schm2_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, double dt, int threads_num = 5);

void acc_pressure_term_update_schm1_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * _neigh_list_cuda, PARTICLE * _particles_cuda, double dt, int blcks = 50); 

void acc_boundaries_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles);
void acc_viscosity_term_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles);
void acc_tension_term_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles);

/* * * * * * * * * * * * * * * * * * * * * */

void velocity_filter1(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles);
void velocity_filter1_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, int threads_num = 5);
void velocity_filter2(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles);
void velocity_filter2_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, int threads_num = 5);

void velocity_filter_schm1_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * _neigh_list_cuda, PARTICLE * _particles_cuda, int blcks = 50);

void velocity_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, double dt);
void velocity_update_schm1_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, double dt, int threads_num = 5);

void velocity_update_schm1_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * _neigh_list_cuda, PARTICLE * _particles_cuda, double dt, int blcks = 50); 

void velocity_update_leap_frog_schm1a(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles, double dt);
void velocity_update_leap_frog_schm1b(vector<vector<int>> neighbors_list, vector<PARTICLE> * particles, double dt);

/* * * * * * * * * * * * * * * * * * * * * */

void position_time_update_schm1(vector<PARTICLE> & particles, double dt);

void position_time_update_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * _neigh_list_cuda, PARTICLE * _particles_cuda, double dt, int blcks = 50);

/* * * * * * * * * * * * * * * * * * * * * */

void support_radius_update_schm1(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles);
void support_radius_update_schm2_multi_threads(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, double dt, int threads_num = 5);
void support_radius_update_schm2(vector<vector<int>> neighbors_list, vector<PARTICLE> & particles, double dt);

void support_radius_update_schm1_cuda(unsigned int neigh_list_length, unsigned int neigh_list_width, unsigned int * _neigh_list_cuda, PARTICLE * _particles_cuda, double dt, int blcks = 50);

/* * * * * * * * * * * * * * * * * * * * * */

#endif