#include "computational/solvers/direct_solver.h"
#include <lapacke.h>
#include "computational/GreensFunction.h"
#include "functions_DSGF.h" // Definitions of functions used in DSGF

void A_b_direct_populator(int tot_sub_vol, double complex A[tot_sub_vol][tot_sub_vol][3][3], double complex G_0[tot_sub_vol][tot_sub_vol][3][3], double complex A_direct[3*3*tot_sub_vol*tot_sub_vol], double complex b_direct[3*3*tot_sub_vol*tot_sub_vol]){

	int ipack=0; 
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
		{
			for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
			{
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
				{
					A_direct[ipack] = A[ig_0][jg_0][i_subG_0][j_subG_0];
					b_direct[ipack]= G_0[ig_0][jg_0][i_subG_0][j_subG_0];
					ipack++;
				}    
			}
		}        
	}   

}

void A_direct_populator(int tot_sub_vol, double complex A[tot_sub_vol][tot_sub_vol][3][3], double complex A_direct[3*3*tot_sub_vol*tot_sub_vol]){

	int ipack=0; 
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
		{
			for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
			{
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
				{
					A_direct[ipack] = A[ig_0][jg_0][i_subG_0][j_subG_0];
					ipack++;
				}    
			}
		}        
	}   

}

void b_direct_populator(int tot_sub_vol,  double complex G_0[tot_sub_vol][tot_sub_vol][3][3], double complex b_direct[3*3*tot_sub_vol*tot_sub_vol]){

	int ipack=0; 
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
		{
			for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
			{
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
				{
					b_direct[ipack]= G_0[ig_0][jg_0][i_subG_0][j_subG_0];
					ipack++;
				}    
			}
		}        
	}   

}


void populate_G_sys(int tot_sub_vol, double complex b_direct[3*3*tot_sub_vol*tot_sub_vol], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol]){

	int gpack=0;
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
		{
			int ig_0_2d = (3*ig_0 + i_subG_0);
			for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
			{
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
				{
					int jg_0_2d = (3*jg_0 + j_subG_0);
					G_sys[ig_0_2d][jg_0_2d] = b_direct[gpack];
					gpack++;
				}    
			}
		}        
	} 
}

/*
void direct_solver(int tot_sub_vol, double complex A[tot_sub_vol][tot_sub_vol][3][3], double complex G_0[tot_sub_vol][tot_sub_vol][3][3], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol]){

	double complex (*A_direct) = calloc(3*3*tot_sub_vol*tot_sub_vol, sizeof(*A_direct));
	double complex (*b_direct) = calloc(3*3*tot_sub_vol*tot_sub_vol, sizeof(*b_direct));

	A_b_direct_populator(tot_sub_vol, A, G_0, A_direct, b_direct);

	// F08ANF (ZGELS) solves linear least-squares problems using a QR or LQ factorization of A
	//Description of ZGELS: https://extras.csc.fi/math/nag/mark21/pdf/F08/f08anf.pdf
	int info = LAPACKE_zgels(LAPACK_ROW_MAJOR,'N',3*tot_sub_vol,3*tot_sub_vol,3*tot_sub_vol,A_direct,3*tot_sub_vol,b_direct,3*tot_sub_vol); 

	populate_G_sys(tot_sub_vol, b_direct, G_sys);


	free(A_direct);
	free(b_direct);	
}
*/

//void direct_solver(int tot_sub_vol, double complex A_direct[3*3*tot_sub_vol*tot_sub_vol],double complex b_direct[3*3*tot_sub_vol*tot_sub_vol], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol]){
//void direct_solver(int tot_sub_vol, double complex G_sys[3*tot_sub_vol][3*tot_sub_vol],double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char wave_type, double complex alpha_0[tot_sub_vol]){
void direct_solver(int tot_sub_vol, double complex G_sys[3*tot_sub_vol][3*tot_sub_vol],double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char wave_type, double complex alpha_0[tot_sub_vol]){
	double complex (*G_0)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*G_0)); 		
	double complex (*A)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*A));
	
	//get_G0_A_matrices(tot_sub_vol, G_0, A, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, alpha_0, delta_V_vector, wave_type); //this function is only applicable in the case of uniform discretization.

	//get_G0_matrix(tot_sub_vol, G_0, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, wave_type);
	get_G0_matrix_memory(tot_sub_vol, G_0, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, wave_type);
	
	get_A_matrix(tot_sub_vol, G_0, A, k_0, alpha_0); // function applicable for uniform and non-uniform discretization

	double complex (*A_direct) = calloc(3*3*tot_sub_vol*tot_sub_vol, sizeof(*A_direct));
	A_direct_populator(tot_sub_vol, A, A_direct);
	free(A);
			
	double complex (*b_direct) = calloc(3*3*tot_sub_vol*tot_sub_vol, sizeof(*b_direct));
	b_direct_populator(tot_sub_vol, G_0, b_direct);
	free(G_0);

	// F08ANF (ZGELS) solves linear least-squares problems using a QR or LQ factorization of A
	//Description of ZGELS: https://extras.csc.fi/math/nag/mark21/pdf/F08/f08anf.pdf
	int info = LAPACKE_zgels(LAPACK_ROW_MAJOR,'N',3*tot_sub_vol,3*tot_sub_vol,3*tot_sub_vol,A_direct,3*tot_sub_vol,b_direct,3*tot_sub_vol); 
	
	free(A_direct);

	populate_G_sys(tot_sub_vol, b_direct, G_sys);
	
	free(b_direct);

}
