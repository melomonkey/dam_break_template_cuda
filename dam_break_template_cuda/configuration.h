#ifndef configuration_h_
#define configuration_h_

#ifndef PI
#define PI 3.1415927
#endif


#ifndef GRAVITY
#define GRAVITY 9.81
#endif

/* * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef REF_PRESSURE
#define REF_PRESSURE 0.0
#endif

#ifndef REF_DENSITY
#define REF_DENSITY 1000.0
#endif

#ifndef COMPRS_LIM_CONST_K
#define COMPRS_LIM_CONST_K 600.0
#endif

#ifndef COMPRS_LIM_CONST_B
#define COMPRS_LIM_CONST_B 92040
#endif

#ifndef INI_TEN_DIR_LEN
#define INI_TEN_DIR_LEN 4e-4
#endif

#ifndef INI_TENSION_ACC
#define INI_TENSION_ACC 1.4e-3
#endif

#ifndef WALL_BOUNCE_CONST
#define WALL_BOUNCE_CONST 80
#endif

/* * * * * * * * * * * * * * * * * * * * * * * * */

#define DENS_CHANGE_LIM 3.0
#define DENS_LIM_U 1010.0
#define DENS_LIM_L 1000.0


#define PRESSURE_LIM_U 7000
#define PRESSURE_LIM_L 0.0


#define ACC_X_LIM_U 500
#define ACC_X_LIM_L -500

#define ACC_Y_LIM_U 500
#define ACC_Y_LIM_L -500

#define ACC_Z_LIM_U 500
#define ACC_Z_LIM_L -500


#define VEL_X_LIM_U 1e30
#define VEL_X_LIM_L -1e30

#define VEL_Y_LIM_U 1e30
#define VEL_Y_LIM_L -1e30

#define VEL_Z_LIM_U 1e30
#define VEL_Z_LIM_L -1e30


#ifndef BOX_CONTAIN // 0 -- No box; 1 -- Box; 
#define BOX_CONTAIN 1
#endif

#ifndef CASE_DIM // 2D -- 0; 3D -- 1;
#define CASE_DIM 0
#endif


#define BOX_X_MIN 0.0
#define BOX_X_MAX 1.6
#define BOX_Y_MIN 0.0
#define BOX_Y_MAX 0.6
#define BOX_Z_MIN 0.0
#define BOX_Z_MAX 1.0

/* * * * * * * * * * * * * * * * * * * * * * * * */

#endif