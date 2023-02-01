#include "material/SiO2.h"
#include "math.h"

double complex calculate_epsilon_SiO2(double q, double omega, double h_bar){

	// Precision of this code is defined on the dielectric function, with 5 digits of precision.
	// Dielectric function for SiO2
	// Defined in Czapla and Narayanaswamy, JQSRT, 2019. DOI: https://doi.org/10.1016/j.jqsrt.2019.01.020 
	double epsilon_inf = 2.03843;
	double h_bar_omega_0[] = {0.05624, 0.09952, 0.13355};    // [eV]
	double S[] = {0.93752, 0.05050, 0.60642};                // [-]
	double Gamma[] = {0.09906, 0.05511,0.05246};
	double omega_0[3]; 
	double epsilon_term1[3];
	double epsilon_term2[3];
	double complex epsilon_num;
	double epsilon_denom; 
	double complex epsilon = epsilon_inf + 0.*I ;
	for(int i_epsilon = 0; i_epsilon < 3; i_epsilon++)
	{
		omega_0[i_epsilon] = h_bar_omega_0[i_epsilon]*q/h_bar;
		epsilon_term1[i_epsilon] = pow(omega/omega_0[i_epsilon],2);
		epsilon_term2[i_epsilon] = Gamma[i_epsilon]*omega/omega_0[i_epsilon];
		epsilon_num = 0. + 0.*I;
		epsilon_num = S[i_epsilon]*(1-epsilon_term1[i_epsilon]) + S[i_epsilon]*epsilon_term2[i_epsilon]*I;
		epsilon_denom=0.;
		epsilon_denom = pow(1. - epsilon_term1[i_epsilon],2) + pow(epsilon_term2[i_epsilon],2); 
		epsilon += epsilon_num/epsilon_denom ;  
	}
	return epsilon;
}
