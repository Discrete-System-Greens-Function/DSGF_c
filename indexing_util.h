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
 * for this function to work the 4D matrix has to be [dim1][dim1][dim2][dim2]
 * this means the 2D matrix would be [dim1*dim2][dim1*dim2]
 *
 * the offset input is the same as dim2
 * for example: if the offset is 3 then the 4D matrix should be of dimensions [dim1][dim1][3][3]
 * 		this means the 2D matrix would be [3*dim1][3*dim1]
 *
 * new_x = offset*major_x + minor_x
 * new_y = offset*major_y + minor_y
 * 
 * the function sets the new_x and new_y properties of the index_map struct that gets passed into it
 *
 * INPUTS:
 * 	- *index_map -> struct which holds the new coordinates
 * 	- int offset -> same as dim2 in the 4D matrix
 * 	- int major_x, int major_y, int minor_x, int minor_y -> coordinates for the point in the 4D matrix [major_x][major_y][minor_x][minor_y]
 */
void four_d_two_d_mapping(index_map *index_map_struct, int offset, int major_x, int minor_x, int major_y, int minor_y);

/*
 * this function maps a 2D coordinate to a 1D coordinate
 * 
 * the 2D matrix can be [dim1][dim2] where dim1 != dim2
 * then !d matrix being mapped to is [dim1*dim2]
 *
 * it has an option for x_leading -> true: the rows get pasted first - offset needs to be dim2
 * 				     false: the columns get pasted first - offset needs to be dim1 
 *
 * INPUTS:
 * 	- *index_map -> struct which holds the mapped coordinates (only new_x will be modified)
 * 	- int offset -> dim1 or dim2 from the 2D matrix
 * 	- int x, int y -> coordinate of the point in the 2D matrix
 */
void two_d_one_d_mapping(index_map *indexing_map_struct, int offset, int x, int y, bool x_leading);

#endif
