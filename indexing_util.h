#ifndef __indexing_util_h
#define __indexing_util_h

#include <stdlib.h>

typedef struct {
	int new_x;
	int new_y;
} index_map;

index_map *init_index_map();

void four_d_two_d_mapping(index_map *index_map_struct, int major_x, int minor_x, int major_y, int minor_y);

#endif
