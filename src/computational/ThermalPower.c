#include "computational/ThermalPower.h"
#include <complex.h>
#include <math.h>
#include "functions_DSGF.h"

void spectral_post_processing(int tot_sub_vol, int i_omega, int const_N_omega, double k_0, double h_bar, double k_b, double complex epsilon, double omega_value, double T_vector[],double delta_V_vector[], double const_N_subvolumes_per_object, double pi,double complex G_sys[3*tot_sub_vol][3*tot_sub_vol], double *sum_trans_coeff, double Q_subvol[tot_sub_vol][const_N_omega]){



	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		double complex Q_omega_subvol = 0;
		for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
		{
			double complex G_element = 0;
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					// Extract one Green's function component: G_sys[ig_0][jg_0][i_subG_0]
					double complex transpose_G_sys = G_sys[ig_0_2d][jg_0_2d];
					double complex G_sys_cross = conj(transpose_G_sys);
					G_element+=G_sys[ig_0_2d][jg_0_2d]*G_sys_cross; 

				}    
			}
			// Transmissivity coefficient matrix tau(omega) for comparison with Czapla Mathematica output [dimensionless]
			double complex trans_coeff = 4.*pow(k_0,4)*delta_V_vector[ig_0]*delta_V_vector[jg_0]*cimag(epsilon)*cimag(epsilon)*G_element; 
			double inner_sum;
			// Thermal power dissipated calculation, based on the matlab code (Using Tervo's Eq. 26)
			if (ig_0 != jg_0) 
			{
				inner_sum = (theta_function(omega_value, T_vector[jg_0], h_bar, k_b) - theta_function(omega_value, T_vector[ig_0], h_bar, k_b)) * trans_coeff;
			}
			else {
				inner_sum = 0;
			}

			//Trans_bulk: Transmission coefficient between two bulk objects
			// This function calculates the transmission coefficient between bulk objects given the transmission coefficient between every pair of dipoles for a given frequency.
			if(ig_0 < const_N_subvolumes_per_object && jg_0 >= const_N_subvolumes_per_object)// bulk 1 -> 2
			{
				*sum_trans_coeff += trans_coeff;
			} 

			Q_omega_subvol += (1 / (2 * pi)) * inner_sum; // calculates the thermal power dissipated per subvolume

		} 
		Q_subvol[ig_0][i_omega] = Q_omega_subvol;
	}       
}
