#include "computational/solvers/direct_solver.h"
#include <lapacke.h> // for desktop
// #include <mkl_lapacke.h> // for CHPC
#include "computational/GreensFunction.h"
#include "functions_DSGF.h" // Definitions of functions used in DSGF

void A_direct_populator_2D(int tot_sub_vol, double complex A[3*tot_sub_vol][3*tot_sub_vol], double complex A_direct[3*3*tot_sub_vol*tot_sub_vol],char multithread){

	#pragma omp parallel for if (multithread == 'Y')
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
		{
			int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
			for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
			{
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					int index = j_subG_0+jg_0*3+i_subG_0*3*tot_sub_vol+ig_0*3*tot_sub_vol*3;
					A_direct[index] = A[ig_0_2d][jg_0_2d];
				}    
			}
		}        
	}   

}


void b_direct_populator_2D(int tot_sub_vol,  double complex G_0[3*tot_sub_vol][3*tot_sub_vol], double complex b_direct[3*3*tot_sub_vol*tot_sub_vol], char multithread){

	#pragma omp parallel for if (multithread == 'Y')
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
		{
			int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
			for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
			{
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					int index = j_subG_0+jg_0*3+i_subG_0*3*tot_sub_vol+ig_0*3*tot_sub_vol*3;
					//printf("%d - ",index);
					b_direct[index]= G_0[ig_0_2d][jg_0_2d];
					//b_direct[ipack]= G_0[ig_0_2d][jg_0_2d];
					//ipack++;
				}    
			}
		}        
	}   
	//printf("\n");
}


void populate_G_sys(int tot_sub_vol, double complex b_direct[3*3*tot_sub_vol*tot_sub_vol], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol], char multithread){

	#pragma omp parallel for if (multithread == 'Y')
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
					int index = j_subG_0+jg_0*3+i_subG_0*3*tot_sub_vol+ig_0*3*tot_sub_vol*3;
					G_sys[ig_0_2d][jg_0_2d] = b_direct[index];
				}    
			}
		}        
	} 
}

void direct_solver_memory(int tot_sub_vol, double complex A_direct[3*tot_sub_vol*3*tot_sub_vol], double complex b_direct[3*tot_sub_vol*3*tot_sub_vol]){
	
	// F08ANF (ZGELS) solves linear least-squares problems using a QR or LQ factorization of A. Description of ZGELS: https://extras.csc.fi/math/nag/mark21/pdf/F08/f08anf.pdf
	int info = LAPACKE_zgels(LAPACK_ROW_MAJOR,'N',3*tot_sub_vol,3*tot_sub_vol,3*tot_sub_vol,A_direct,3*tot_sub_vol,b_direct,3*tot_sub_vol); 

}