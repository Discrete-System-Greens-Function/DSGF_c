#ifndef __indexing_util_h
#define __indexing_util_h

#include <stdlib.h>
#include <stdbool.h>

typedef struct {
	int new_x;
	int new_y;
} index_map;

index_map *init_index_map();

/*
 * this function maps a 4D coordinate to a 2D coordinate (only works for square matrices!!)
 * 
 * the offset input dictates the dimensions for the new matrix
 * for example: if the offset is 3 then the 2D matrix should be of dimensions [3*outer_dim]
 * TODO: this description needs to be finished
 * new_x = 3*major_x + minor_x
 * new_y = 3*major_y + minor_y
 *
 *
 */
void four_d_two_d_mapping(index_map *index_map_struct, int major_x, int minor_x, int major_y, int minor_y);

/*
 * this function maps a 2D coordinate to a 1D coordinate
 *
 * it has an option for x_leading -> if this is enabled the x coordinate is used for the offset
 * TODO: this description needs to be finished
 */
void two_d_one_d_mapping(index_map *indexing_map_struct, int x, int y, bool x_leading);

#endif
