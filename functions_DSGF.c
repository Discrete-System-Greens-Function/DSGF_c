// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Definition of some functions in DSGF framework
// Developed by RETL Lab, Department of Mechanical Engineering, The University of Utah, USA.
// LAST UPDATE: JUNE 20, 2022  
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "functions_DSGF.h"

struct stat st = {0};

double k_0_function(double omega, double epsilon_0, double mu_0){
	return omega*sqrt(epsilon_0*mu_0);
}

double k_function(double omega, double epsilon_ref, double epsilon_0, double mu_0){
	return omega*sqrt(epsilon_ref*epsilon_0*mu_0);	
}

double vol_sphere(double radius, double pi){
	return (4./3)*pi*pow(radius, 3);
}

double a_j_function(double delta_V_vector, double pi){
	return pow( (3. * delta_V_vector)/(4*pi), 1./3);
}

double theta_function(double omega, double T, double h_bar, double k_b){
	return (h_bar*omega)/(exp(h_bar*omega/(k_b*T))-1);
}

double  dtheta_dT_function(double omega,double T_calc, double h_bar, double k_b){
	double psi = h_bar*(omega)/(k_b*T_calc);     // Input for exponential dtheta_dT function[dimensionless]
	double dtheta_dT = (k_b*pow(psi,2)*exp(psi))/(pow(exp(psi)-1.,2)); // [J/K]
	return dtheta_dT;           // return statement
}

long get_mem_usage() 
{
	struct rusage myusage;
	getrusage(RUSAGE_SELF,&myusage); // getrusage() https://man7.org/linux/man-pages/man2/getrusage.2.html
	return myusage.ru_maxrss;           // return statement

}
