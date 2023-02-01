#include "material/SiN.h"
#include <math.h>

double complex calculate_epsilon_SiN(double omega, double pi){

	// Dielectric function for SiN -  Lorentz oscillator model:
	// Defined in Cataldo et al, Optics letters 2012. DOI: https://doi.org/10.1364/OL.37.004200
	int M=5;
	double complex epsilon_j[] = {7.582+ 0*I, 6.754+0.3759*I, 6.601+ 0.0041*I, 5.430+0.1179*I, 4.601+0.2073*I, 4.562 + 0.0124*I};
	double complex epsilon_inf = epsilon_j[5];
	double omega_T[] = {13.913, 15.053, 24.521, 26.440, 31.724};
	double Gamma[] = {5.810, 6.436, 2.751, 3.482, 5.948};
	double alpha[] = {0.0001, 0.3427, 0.0006, 0.0002, 0.0080};
	double omega_T_converted, Gamma_converted;
	double complex delta_epsilon, num_1, denom_1,Gamma_prime,num_2, denom_2;

	double complex summation = 0. + 0.*I;
	for(int i_epsilon = 0; i_epsilon < M; i_epsilon++)
	{				
		omega_T_converted = omega_T[i_epsilon]*(2*pi)*(1.e12);
		Gamma_converted =Gamma[i_epsilon]*(2*pi)*(1.e12);
		num_1 = pow(omega_T_converted,2)- pow(omega,2);
		denom_1 = omega*Gamma_converted;

		Gamma_prime = Gamma_converted*cexp(-alpha[i_epsilon]*pow(num_1/denom_1,2));

		delta_epsilon = epsilon_j[i_epsilon] - epsilon_j[i_epsilon+1]; 

		num_2 = delta_epsilon*pow(omega_T_converted,2);
		denom_2 = pow(omega_T_converted,2) - pow(omega,2) - omega*Gamma_prime*I;
		summation += num_2/denom_2;

	}
	double complex epsilon = epsilon_inf + summation;
	return epsilon;
}
