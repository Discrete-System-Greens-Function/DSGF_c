#include "computational/solvers/direct_solver.h"


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
