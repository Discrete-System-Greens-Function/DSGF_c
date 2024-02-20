// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Functions used in the Discrete System Green's Function 
// Developed by RETL Lab, Department of Mechanical Engineering, The University of Utah, USA.
// LAST UPDATE: MAY 31, 2022
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifndef __functions_dsgf_h //https://stackoverflow.com/questions/28361391/calling-define-from-another-file
#define __functions_dsgf_h

#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <sys/resource.h>
#include <sys/stat.h>
// ########## dsgf_functions.h ##########

// Function for wave vector in free space calculation:
double k_0_function(double omega, double epsilon_0, double mu_0);


// Function for wave vector in background reference medium calculation:
double k_function(double omega, double epsilon_ref, double epsilon_0, double mu_0);

// Function for the thermal volumes calculation:
double vol_sphere(double radius, double pi);

// Definition used to import the discretization positions
typedef struct node {
	float x;
	float y;
	float z;
} subvol;


// Definition for Green's tensors
struct Matrix3x3 {
	double complex xx, xy, xz;
	double complex yx, yy, yz;
	double complex zx, zy, zz;
};

struct Eye3x3 {
	double  xx, xy, xz;
	double  yx, yy, yz;
	double  zx, zy, zz;
};


// Function a for free space green's function calculation:
double a_j_function(double delta_V_vector, double pi);

// Function for the mean energy of an electromagnetic state calculation:
double theta_function(double omega, double T, double h_bar, double k_b);

// Function for the calculation of the derivative of the mean energy of an electromagnetic state by the temperature
double  dtheta_dT_function(double omega,double T_calc, double h_bar, double k_b);

// How to measure memory usage inside my program? (getrusage): https://youtu.be/Os5cK0H8EOA
long get_mem_usage();

// Code's output files definitions
char matrices_folder[100];
char frequency_folder[100];
char spectral_transmissivity_folder[100];
char spectral_conductance_folder[100];
char dirPathpos_processing_summary_FileName[260];

// ########## end dsgf_functions.h ##########

#endif // __functions_dsgf_h
