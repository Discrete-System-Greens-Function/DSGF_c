#ifndef __iterative_soler_h
#define __iterative_soler_h

#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "indexing_util.h"
#include "debugging/debugging_utils.h"

void matrix_reshape(int inner_size, int outer_size, double complex matrix_2d_1[][3*outer_size], double matrix_2d_2[][3*outer_size], double complex matrix_4d_1[outer_size][outer_size][inner_size][inner_size], double matrix_4d_2[outer_size][outer_size][inner_size][inner_size]);


void A2d_solver(double complex epsilon, int mm, int tot_sub_vol, double eyeA_2d[][3*tot_sub_vol], double delta_V, double complex G_sys_old[][3* tot_sub_vol], double complex A_2d[][3*tot_sub_vol], double k);

void A1_lapack_B1_lapack_populator(int tot_sub_vol ,double complex *A1lapack, double complex *b1lapack, double complex A_2d[3][3], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], int mm, int jg_0);

#endif
