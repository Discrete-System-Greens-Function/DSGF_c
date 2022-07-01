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
					four_d_two_d_mapping(indexing_map, major_x, minor_x, major_y, minor_y);
					matrix_2d_1[indexing_map->new_x][indexing_map->new_y]=matrix_4d_1[major_x][major_y][minor_x][minor_y]; 
					matrix_2d_2[indexing_map->new_x][indexing_map->new_y] = matrix_4d_2[major_x][major_y][minor_x][minor_y];
				}
			}
		}
	}

	free(indexing_map);

}
