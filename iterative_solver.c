#include "iterative_solver.h"

typedef struct {
	int new_x;
	int new_y;
} mapping;

mapping *init_mapping(){
	mapping *new_map = malloc(sizeof(mapping));
	new_map->new_x = -1;
	new_map->new_y = -1;
	return new_map;
}

void get_new_point(mapping *mapping_struct, int major_x, int minor_x, int major_y, int minor_y){
	mapping_struct->new_x = 3*major_x + minor_x;
	mapping_struct->new_y = 3*major_y + minor_y;
}

void matrix_reshape(int inner_size, int outer_size, double complex matrix_2d_1[][3*outer_size], double matrix_2d_2[][3*outer_size], double complex matrix_4d_1[outer_size][outer_size][inner_size][inner_size], double matrix_4d_2[outer_size][outer_size][inner_size][inner_size]){
	
	mapping *new_map = init_mapping();

	for (int major_x = 0; major_x < outer_size; major_x++)
	{
		for (int major_y = 0; major_y < outer_size; major_y++)
		{ 
			for(int minor_x = 0; minor_x < inner_size; minor_x++)
			{
				for(int minor_y = 0; minor_y < inner_size; minor_y++)
				{
					get_new_point(new_map, major_x, minor_x, major_y, minor_y);
					matrix_2d_1[new_map->new_x][new_map->new_y]=matrix_4d_1[major_x][major_y][minor_x][minor_y]; 
					matrix_2d_2[new_map->new_x][new_map->new_y] = matrix_4d_2[major_x][major_y][minor_x][minor_y];
				}
			}
		}
	}

	free(new_map);

}
