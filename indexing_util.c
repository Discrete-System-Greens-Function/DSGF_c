#include "indexing_util.h"

index_map *init_index_map(){
	index_map *new_index_map = malloc(sizeof(index_map));
	new_index_map->new_x = -1;
	new_index_map->new_y = -1;
	return new_index_map;
}

void four_d_two_d_mapping(index_map *index_map_struct, int major_x, int minor_x, int major_y, int minor_y){
	
	index_map_struct->new_x = 3*major_x + minor_x;
	index_map_struct->new_y = 3*major_y + minor_y;

}

void two_d_one_d_mapping(index_map *indexing_map_struct, int x, int y, bool x_leading){

	
}
