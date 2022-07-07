#include "iterative_solver.h"

void matrix_reshape(int inner_size, int outer_size, double complex matrix_2d_1[][3*outer_size], double matrix_2d_2[][3*outer_size], double complex matrix_4d_1[outer_size][outer_size][inner_size][inner_size], double matrix_4d_2[outer_size][outer_size][inner_size][inner_size]){
	
	index_map *indexing_map = init_index_map();

	for (int major_x = 0; major_x < outer_size; major_x++)
	{
		for (int major_y = 0; major_y < outer_size; major_y++)
		{ 
			for(int minor_x = 0; minor_x < inner_size; minor_x++)
			{
				for(int minor_y = 0; minor_y < inner_size; minor_y++)
				{
					four_d_two_d_mapping(indexing_map, 3, major_x, minor_x, major_y, minor_y);
					matrix_2d_1[indexing_map->new_x][indexing_map->new_y]=matrix_4d_1[major_x][major_y][minor_x][minor_y]; 
					matrix_2d_2[indexing_map->new_x][indexing_map->new_y] = matrix_4d_2[major_x][major_y][minor_x][minor_y];
				}
			}
		}
	}

	free(indexing_map);

}

void A2d_solver(double complex epsilon, int mm, int tot_sub_vol, double eyeA_2d[][3*tot_sub_vol], double delta_V, double complex G_sys_old[][3*tot_sub_vol], double complex A_2d[3][3], double k){
	
	index_map *new_index_map = init_index_map();

	for(int minor_x = 0; minor_x < 3; minor_x++){
		for (int minor_y = 0; minor_y < 3; minor_y++){
			four_d_two_d_mapping(new_index_map, 3, mm, minor_x, mm, minor_y);

			A_2d[minor_x][minor_y] = eyeA_2d[new_index_map->new_x][new_index_map->new_y] - pow(k,2)*delta_V*epsilon*G_sys_old[new_index_map->new_x][new_index_map->new_y]; //modification...see if it works
		}
	}

	free(new_index_map);
}
