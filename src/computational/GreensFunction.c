#include "computational/GreensFunction.h"
#include <math.h>

void setup_G_0_matrices(int tot_sub_vol, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double R[][3]){


	double r[tot_sub_vol][tot_sub_vol][3];
	double abs_r_ij[tot_sub_vol][tot_sub_vol][3];
	double unit_r_ij[tot_sub_vol][tot_sub_vol][3];
	double complex unit_conj_r_ij[tot_sub_vol][tot_sub_vol][3];

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

}
