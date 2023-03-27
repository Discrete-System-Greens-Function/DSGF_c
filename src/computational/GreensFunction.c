#include "computational/GreensFunction.h"
#include <math.h>
#include "functions_DSGF.h"
#include <stdlib.h>

void setup_G_0_matrices(int tot_sub_vol, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double R[][3]){


	double (*r)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*r));
	double (*abs_r_ij)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*abs_r_ij));
	double (*unit_r_ij)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*unit_r_ij));
	double complex (*unit_conj_r_ij)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*unit_conj_r_ij));

	for (int i=0; i<tot_sub_vol; i++)
	{   
		for (int j=0; j<tot_sub_vol; j++)
		{
			for (int i_alpha=0; i_alpha<3; i_alpha++)
			{
				r[i][j][i_alpha]= R[i][i_alpha] - R[j][i_alpha] ; //general case 
				abs_r_ij[i][j][i_alpha] = fabs(r[i][j][i_alpha]);
				modulo_r_i_j[i][j] += pow(abs_r_ij[i][j][i_alpha],2);
			}
			modulo_r_i_j[i][j] = sqrt(modulo_r_i_j[i][j]);
		}
	}

	for (int i_i = 0; i_i < tot_sub_vol; i_i++)
	{ 
		for (int i_j = 0; i_j < tot_sub_vol; i_j++)
		{ 
			for (int i_alpha = 0; i_alpha < 3; i_alpha++)
			{  
				unit_r_ij[i_i][i_j][i_alpha] = r[i_i][i_j][i_alpha]/modulo_r_i_j[i_i][i_j]; // ˆr -- unit distance 
				double transpose = r[i_i][i_j][i_alpha]/modulo_r_i_j[i_i][i_j]; // https://www.programiz.com/c-programming/examples/matrix-transpose        
				unit_conj_r_ij[i_i][i_j][i_alpha] = conj(transpose);// r^+ Conjugate transpose unit distance 
			}
		}
	}

	for (int i_i = 0; i_i < tot_sub_vol; i_i++)
	{ 
		for (int i_j = 0; i_j < tot_sub_vol; i_j++)
		{ 
			for (int i_alpha = 0; i_alpha < 3; i_alpha++)
			{  
				for (int j_alpha = 0; j_alpha < 3; j_alpha++)
				{ 
					r_i_j_outer_r_i_j[i_i][i_j][i_alpha][j_alpha] = unit_r_ij[i_i][i_j][i_alpha]*unit_conj_r_ij[i_i][i_j][j_alpha];

				}
			}    
		}
	}

	free(r);
	free(abs_r_ij);
	free(unit_r_ij);
	free(unit_conj_r_ij);
}

void get_G0_A_matrices(int tot_sub_vol, double complex G_0[tot_sub_vol][tot_sub_vol][3][3], double complex A[tot_sub_vol][tot_sub_vol][3][3], double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double complex alpha_0[tot_sub_vol], double delta_V_vector[tot_sub_vol]){

	// ################### MATRICES STRUCTURE LOOPS ###########################
	// 3N X 3N Matrices structure loops for G^0 and A:
	double denom_1, denom_2 ; // used in G^0_ij function
	double complex const_1, const_2, const_3;

	double eyeG_0[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};


	// eq. 25 from Lindsay's paper 

	//for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++)//
	for (int jg_0 = 0; jg_0 < tot_sub_vol-1; jg_0++) //tot_sub_vol
	{
		//for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++)//
		for (int ig_0 = jg_0+1; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
		{
			if (ig_0!=jg_0)
			{
			const_1 = cexp(k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0]*I)/(4.*pi*modulo_r_i_j[ig_0][jg_0]); 
			denom_1 = epsilon_ref*pow(k_0*modulo_r_i_j[ig_0][jg_0],2);
			denom_2 = k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0];
			const_2 = (1. - 1./denom_1 + 1.*I/denom_2 ) ;
			const_3 = (1. - 3./denom_1 + 3.*I/denom_2) ;
			/*
			//split of G_0:

			if (jg_0<=tot_sub_vol/2 && ig_0>tot_sub_vol/2) //subvolumes in different objects
			{
				const_2 = (1. - 1./denom_1 + 1.*I/denom_2 ) ; // total 
				const_3 = (1. - 3./denom_1 + 3.*I/denom_2) ;  // total 
				//Goal: compute only the propagating wave contribution of DSGF
				//const_2 = (1. ) ; //propagating only
				//const_3 = (1. ) ; //propagating only
				//Goal: compute only the evanescent wave contribution of DSGF
				//const_2 = (- 1./denom_1 + 1.*I/denom_2 ) ; //evanescent only
				//const_3 = (- 3./denom_1 + 3.*I/denom_2) ; //evanescent only
			}
			else 
			{
				const_2 = (1. - 1./denom_1 + 1.*I/denom_2 ) ;
				const_3 = (1. - 3./denom_1 + 3.*I/denom_2) ; 
			}
			*/

			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					G_0[ig_0][jg_0][i_subG_0][j_subG_0] = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
					G_0[jg_0][ig_0][i_subG_0][j_subG_0] = G_0[ig_0][jg_0][i_subG_0][j_subG_0];
					A[ig_0][jg_0][i_subG_0][j_subG_0] = 0 - pow(k_0,2)*alpha_0[ig_0]*G_0[ig_0][jg_0][i_subG_0][j_subG_0]; 
					A[jg_0][ig_0][i_subG_0][j_subG_0] =A[ig_0][jg_0][i_subG_0][j_subG_0];
				}    
			}
			}
		}    
	} //end j_subG_0
	// eq. 26 from Lindsay's paper: 


	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		double a_j = a_j_function(delta_V_vector[ig_0], pi);
		double part1ii = 1./(3.*delta_V_vector[ig_0]*epsilon_ref*pow(k_0,2));
		double complex part2ii = a_j*k_0*sqrt(epsilon_ref)*I; // com i term 
		double complex part2iiexp = cexp(0. + a_j*k_0*sqrt(epsilon_ref)*I); 
		// part3ii is inside brackets
		double complex part3ii = part2iiexp*(1-part2ii) - 1. ;
		for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
		{
			for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
			{ 
				if (ig_0==jg_0) // if i=j:
				{
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
					{
						G_0[ig_0][jg_0][i_subG_0][j_subG_0] = eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.); 
						A[ig_0][jg_0][i_subG_0][j_subG_0] = eyeG_0[i_subG_0][j_subG_0] - pow(k_0,2)*alpha_0[ig_0]*G_0[ig_0][jg_0][i_subG_0][j_subG_0]; 
					}
				} 
			}   //end jg_0   
		} //end i_subG_0     
	} //end ig_0 
}
