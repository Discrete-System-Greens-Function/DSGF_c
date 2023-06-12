#ifndef __iterative_soler_h
#define __iterative_soler_h

#include <stdlib.h>
#include <complex.h>
#include <math.h>

void matrix_reshape(int inner_size, int outer_size, double complex matrix_2d_1[][3*outer_size], double complex matrix_4d_1[outer_size][outer_size][inner_size][inner_size]);

void A2d_solver(double complex epsilon, int mm, int tot_sub_vol, double delta_V, double complex G_sys_old[][3* tot_sub_vol], double complex A_2d[][3*tot_sub_vol], double k);

void A_b_iterative_populator(int tot_sub_vol ,double complex *A_iterative, double complex *b_iterative, double complex A_2d[3][3], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], int mm, int jg_0);

void G_sys_new_populator(int tot_sub_vol, int mm, int jg_0, double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], double complex b_iterative[]);

void offdiagonal_solver(int tot_sub_vol, int mm, double k, double complex alpha_0[], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol]);

void remaining_pertubations(int tot_sub_vol, int mm, double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], double complex A_2d[3][3]);

void core_solver(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[],double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol]);

//void iterative_solver(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol]);

void iterative_solver(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi,double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3],char wave_type);


#endif
