#include "material/u_SiC.h"
#include "math.h"

double complex calculate_epsilon_u_SiC(double omega, double pi, double c_0){


	// Dielectric function for u-SiC -  Lorentz oscillator model:
	// Defined in St-Gelais et al., Nature nanotechnology, 2016. DOI: https://doi.org/10.1038/nnano.2016.20    
	double to_rad_s = 100*2*pi*c_0;       
	double epsilon_inf = 8;              // [-]
	double Gamma = 20*to_rad_s; 		//[1/s]  Drude relaxation time
	double omega_TO = 789*to_rad_s;		//[rad/s]		
	double omega_LO = 956*to_rad_s;		//[rad/s]
	double complex epsilon = epsilon_inf*(pow(omega,2)-pow(omega_LO,2)+I*Gamma*omega)/(pow(omega,2)-pow(omega_TO,2)+I*Gamma*omega) ; // omega_TO omega_LO Gamma epsilon_inf

	return epsilon;
}
