#include "computational/solvers/direct_solver.h"
#include <lapacke.h>

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


void direct_solver(int tot_sub_vol, double complex A[tot_sub_vol][tot_sub_vol][3][3], double complex G_0[tot_sub_vol][tot_sub_vol][3][3], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol]){


	double complex (*A_direct) = calloc(3*3*tot_sub_vol*tot_sub_vol, sizeof(*A_direct));
	double complex (*b_direct) = calloc(3*3*tot_sub_vol*tot_sub_vol, sizeof(*b_direct));

	A_b_direct_populator(tot_sub_vol, A, G_0, A_direct, b_direct);


	// F08ANF (ZGELS) solves linear least-squares problems using a QR or LQ factorization of A
	int info = LAPACKE_zgels(LAPACK_ROW_MAJOR,'N',3*tot_sub_vol,3*tot_sub_vol,3*tot_sub_vol,A_direct,3*tot_sub_vol,b_direct,3*tot_sub_vol); 

	populate_G_sys(tot_sub_vol, b_direct, G_sys);

	free(A_direct);
	free(b_direct);	
}
