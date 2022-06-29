#ifndef __iterative_soler_h
#define __iterative_soler_h

#include <stdlib.h>
#include <complex.h>

void matrix_reshape(int inner_size, int outer_size, double complex matrix_2d_1[][3*outer_size], double matrix_2d_2[][3*outer_size], double complex matrix_4d_1[outer_size][outer_size][inner_size][inner_size], double matrix_4d_2[outer_size][outer_size][inner_size][inner_size]);

#endif
