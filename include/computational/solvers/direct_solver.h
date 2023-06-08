#ifndef DIRECT_SOLVER_H
#define DIRECT_SOLVER_H

#include <complex.h>

void A_b_direct_populator(int tot_sub_vol, double complex A[tot_sub_vol][tot_sub_vol][3][3], double complex G_0[tot_sub_vol][tot_sub_vol][3][3], double complex A_direct[3*3*tot_sub_vol*tot_sub_vol], double complex b_direct[3*3*tot_sub_vol*tot_sub_vol]);

void A_direct_populator(int tot_sub_vol, double complex A[tot_sub_vol][tot_sub_vol][3][3], double complex A_direct[3*3*tot_sub_vol*tot_sub_vol]);

void b_direct_populator(int tot_sub_vol, double complex G_0[tot_sub_vol][tot_sub_vol][3][3], double complex b_direct[3*3*tot_sub_vol*tot_sub_vol]);

void populate_G_sys(int tot_sub_vol, double complex b_direct[3*3*tot_sub_vol*tot_sub_vol], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol]);

//void direct_solver(int tot_sub_vol, double complex A[tot_sub_vol][tot_sub_vol][3][3], double complex G_0[tot_sub_vol][tot_sub_vol][3][3], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol]);
// void direct_solver(int tot_sub_vol, double complex A_direct[3*3*tot_sub_vol*tot_sub_vol], double complex b_direct[3*3*tot_sub_vol*tot_sub_vol], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol]);
void direct_solver(int tot_sub_vol, double complex G_sys[3*tot_sub_vol][3*tot_sub_vol],double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char wave_type, double complex alpha_0[tot_sub_vol]);


#endif /* DIRECT_SOLVER_H */
