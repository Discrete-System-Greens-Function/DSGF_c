#include "material/SiC.h"
#include "math.h"

double complex calculate_epsilon_SiC(double omega){


	// Dielectric function for SiC -  Lorentz oscillator model:
	// Defined in Francoeur et al., PRB, 2011. DOI: https://doi.org/10.1103/PhysRevB.84.075436           
	double epsilon_inf = 6.7;              // [-]
	double Gamma = 8.966e11; 		//[1/s]  Drude relaxation time
	double omega_TO = 1.494e14;		//[rad/s]		
	double omega_LO = 1.825e14;		//[rad/s]
	double complex epsilon = epsilon_inf*(pow(omega,2)-pow(omega_LO,2)+I*Gamma*omega)/(pow(omega,2)-pow(omega_TO,2)+I*Gamma*omega) ; // omega_TO omega_LO Gamma epsilon_inf
	return epsilon;
}
