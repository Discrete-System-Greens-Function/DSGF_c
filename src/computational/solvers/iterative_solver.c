#include "computational/solvers/iterative_solver.h"

void matrix_reshape(int inner_size, int outer_size, double complex matrix_2d_1[][3*outer_size], double complex matrix_4d_1[outer_size][outer_size][inner_size][inner_size]){

	for (int major_row = 0; major_row < outer_size; major_row++)
	{
		for (int major_y = 0; major_y < outer_size; major_y++)
		{ 
			for(int minor_x = 0; minor_x < inner_size; minor_x++)
			{
				for(int minor_y = 0; minor_y < inner_size; minor_y++)
				{
					int new_x = 3*major_row + minor_x;
					int new_y = 3*major_y + minor_y;
					matrix_2d_1[new_x][new_y]=matrix_4d_1[major_row][major_y][minor_x][minor_y]; 
				}
			}
		}
	}

}

void A2d_solver(double complex epsilon, int mm, int tot_sub_vol, double eyeA_2d[3][3], double delta_V, double complex G_sys_old[][3*tot_sub_vol], double complex A_2d[3][3], double k){

	for(int minor_x = 0; minor_x < 3; minor_x++){
		int new_x = 3*mm + minor_x;
		for (int minor_y = 0; minor_y < 3; minor_y++){
			int new_y = 3*mm + minor_y;
			A_2d[minor_x][minor_y] = eyeA_2d[minor_x][minor_y] - pow(k,2)*delta_V*epsilon*G_sys_old[new_x][new_y];
		}
	}
}

void A_b_iterative_populator(int tot_sub_vol, double complex *A_iterative, double complex *b_iterative, double complex A_2d[3][3], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol],int mm, int jg_0){


	int ipack=0;

	for (int row = 0; row < 3; row++)
	{
		int new_x = (3*mm + row);
		for(int col = 0; col < 3; col++) // 3D coordinate positions 
		{
			int new_y = (3*jg_0 + col);
			A_iterative[ipack] = A_2d[row][col]; //A_2d[mm_2d][mm_2d]; //A[mm][mm][i_subG_0][j_subG_0];
			b_iterative[ipack] = G_sys_old[new_x][new_y]; //G_sys_old[mm][jg_0][i_subG_0][j_subG_0];
			ipack++;
		}    
	}
}


void G_sys_new_populator(int tot_sub_vol, int mm, int jg_0, double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], double complex b_iterative[]){


	int gpack=0;     	

	for (int mm_sub = 0; mm_sub < 3; mm_sub++) // 3D coordinate positions
	{
		int mm_2d = (3*mm + mm_sub);
		for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
		{
			int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices

			G_sys_new[mm_2d][jg_0_2d] = b_iterative[gpack]; // stores G^1_11

			gpack++;
		}  
	}
}


void offdiagonal_solver(int tot_sub_vol, int mm, double k, double complex alpha_0[], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol]){


	// Next, solve all systems of equations for ii not equal to mm		
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{ 
		for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
		{
			int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
			for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //lower triangular matrix
			{
				if(ig_0 != mm)
				{
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
					{
						int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
						double complex G_sys_prod = 0.;
						for(int m_sub = 0;  m_sub < 3;  m_sub++)//loop for matricial multiplication
						{
							int mm_2d = (3*mm + m_sub);
							G_sys_prod += G_sys_old[ig_0_2d][mm_2d]*G_sys_new[mm_2d][jg_0_2d];
						}

						G_sys_new[ig_0_2d][jg_0_2d] = G_sys_old[ig_0_2d][jg_0_2d] + pow(k,2)*alpha_0[mm]*G_sys_prod;

					} // j_subG_0            				
				} // if(ig_0 != mm) 
			} // jg_0   	
		}// i_subG_0   	                	
	} // ig_0    

}
