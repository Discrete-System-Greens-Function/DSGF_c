#include "iterative_solver.h"

void matrix_reshape(int inner_size, int outer_size, double complex matrix_2d_1[][3*outer_size], double matrix_2d_2[][3*outer_size], double complex matrix_4d_1[outer_size][outer_size][inner_size][inner_size], double matrix_4d_2[outer_size][outer_size][inner_size][inner_size]){

	index_map *indexing_map = init_index_map();

	for (int major_row = 0; major_row < outer_size; major_row++)
	{
		for (int major_y = 0; major_y < outer_size; major_y++)
		{ 
			for(int minor_x = 0; minor_x < inner_size; minor_x++)
			{
				for(int minor_y = 0; minor_y < inner_size; minor_y++)
				{
					four_d_two_d_mapping(indexing_map, 3, major_row, minor_x, major_y, minor_y);
					matrix_2d_1[indexing_map->new_x][indexing_map->new_y]=matrix_4d_1[major_row][major_y][minor_x][minor_y]; 
					matrix_2d_2[indexing_map->new_x][indexing_map->new_y] = matrix_4d_2[major_row][major_y][minor_x][minor_y];
				}
			}
		}
	}

	free(indexing_map);

}

void A2d_solver(double complex epsilon, int mm, int tot_sub_vol, double eyeA_2d[3][3], double delta_V, double complex G_sys_old[][3*tot_sub_vol], double complex A_2d[3][3], double k){

	index_map *new_index_map = init_index_map();

	for(int minor_x = 0; minor_x < 3; minor_x++){
		int new_x = 3*mm + minor_x;
		for (int minor_y = 0; minor_y < 3; minor_y++){
			int new_y = 3*mm + minor_y;
			//four_d_two_d_mapping(new_index_map, 3, mm, minor_x, mm, minor_y);
			//Amm[mm_sub][mm_sub_n] = eye_iter[mm_sub][mm_sub_n] - pow(k,2)*epsilon_s[mm]*delta_V_vector[mm]*G_sys_old[mm_2d][mm_2d_n]; 
			A_2d[minor_x][minor_y] = eyeA_2d[minor_x][minor_y] - pow(k,2)*delta_V*epsilon*G_sys_old[new_x][new_y]; //modification...see if it works
		}
	}

	free(new_index_map);
}

void A1_lapack_B1_lapack_populator(int tot_sub_vol, double complex *A1lapack, double complex *b1lapack, double complex A_2d[3][3], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol],int mm, int jg_0){

	
	index_map *indexing_map_2D = init_index_map();
	index_map *indexing_map_4D = init_index_map();
	
	int ipack=0;

	for (int row = 0; row < 3; row++)
	{
		int new_x = (3*mm + row);
		for(int col = 0; col < 3; col++) // 3D coordinate positions 
		{
			int new_y = (3*jg_0 + col);
			A1lapack[ipack] = A_2d[row][col]; //A_2d[mm_2d][mm_2d]; //A[mm][mm][i_subG_0][j_subG_0];
			b1lapack[ipack] = G_sys_old[new_x][new_y]; //G_sys_old[mm][jg_0][i_subG_0][j_subG_0];
			ipack = ipack + 1;
			//four_d_two_d_mapping(indexing_map_4D, 3, mm, row, jg_0, col);
			//two_d_one_d_mapping(indexing_map_2D, 3, row, col, false);
	
			//A1lapack[indexing_map_2D->new_x] = A_2d[row][col];
			//b1lapack[indexing_map_2D->new_x] = G_sys_old[indexing_map_4D->new_x][indexing_map_4D->new_y];
		}    
	}
	
	free(indexing_map_2D);
	free(indexing_map_4D);

}


