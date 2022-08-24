#include "indexing_util.h"

index_map *init_index_map(){
	index_map *new_index_map = malloc(sizeof(index_map));
	new_index_map->new_x = -1;
	new_index_map->new_y = -1;
	return new_index_map;
}

void four_d_two_d_mapping(index_map *index_map_struct, int offset, int major_x, int minor_x, int major_y, int minor_y){
	
	index_map_struct->new_x = offset*major_x + minor_x;
	index_map_struct->new_y = offset*major_y + minor_y;

}

void two_d_one_d_mapping(index_map *indexing_map_struct, int offset, int row, int col, bool row_leading){

	if (row_leading){
		indexing_map_struct->new_x = offset * row + col;
	}else {
		indexing_map_struct->new_x = offset * col + row;
	}	
}
