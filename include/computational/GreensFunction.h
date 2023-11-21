#ifndef GREENSFUNCTION_H
#define GREENSFUNCTION_H

#include <complex.h>

void set_up_get_G0_1D(int tot_sub_vol, double complex G_0[3*tot_sub_vol*3*tot_sub_vol], double k_0, double pi, double epsilon_ref, double delta_V_vector[tot_sub_vol], double R[][3]);
void get_A_matrix_1D(int tot_sub_vol, double complex G_0[3*tot_sub_vol*3*tot_sub_vol], double complex A[3*tot_sub_vol*3*tot_sub_vol], double k_0, double complex alpha_0[tot_sub_vol]);


void set_up_get_G0_2D(int tot_sub_vol, double complex G_0[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi, double epsilon_ref, double delta_V_vector[tot_sub_vol], double R[][3]);
//void set_up_get_G0_2D(int tot_sub_vol, double complex G_0[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi, double epsilon_ref, double delta_V_vector[tot_sub_vol],char multithread, double R[][3]);

void get_A_matrix_2D(int tot_sub_vol, double complex G_0[3*tot_sub_vol][3*tot_sub_vol], double complex A[3*tot_sub_vol][3*tot_sub_vol], double k_0, double complex alpha_0[tot_sub_vol]);
//void get_A_matrix_2D(int tot_sub_vol, double complex G_0[3*tot_sub_vol][3*tot_sub_vol], double complex A[3*tot_sub_vol][3*tot_sub_vol], double k_0, double complex alpha_0[tot_sub_vol], char multithread);

void set_up_get_G_old(int tot_sub_vol, double complex G_old[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi, double epsilon_ref, double delta_V_vector[tot_sub_vol], double R[][3]);
//void set_up_get_G_old(int tot_sub_vol, double complex G_old[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi, double epsilon_ref, double delta_V_vector[tot_sub_vol],char multithread, double R[][3]);

// old functions:

void setup_G_0_matrices(int tot_sub_vol, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double R[][3], char multithread);

void set_up_get_G_old_file(int tot_sub_vol, double k_0, double pi, double epsilon_ref, double delta_V_vector[tot_sub_vol],char multithread, char* G_old_file_name,double R[][3]);

void get_G0_A_matrices(int tot_sub_vol, double complex G_0[tot_sub_vol][tot_sub_vol][3][3], double complex A[tot_sub_vol][tot_sub_vol][3][3], double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double complex alpha_0[tot_sub_vol], double delta_V_vector[tot_sub_vol], char wave_type);

void get_G0_matrix(int tot_sub_vol, double complex G_0[tot_sub_vol][tot_sub_vol][3][3], double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol], char wave_type);

void get_G0_matrix_memory(int tot_sub_vol, double complex G_0[tot_sub_vol][tot_sub_vol][3][3], double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char wave_type); 

void get_G0_triangular_matrix(int tot_sub_vol, int size, double complex upperTriangularMatrix[size], double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char wave_type);

void get_G_old_matrix_memory(int tot_sub_vol, double complex G_old[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char wave_type);

void get_A_matrix(int tot_sub_vol, double complex G_0[tot_sub_vol][tot_sub_vol][3][3], double complex A[tot_sub_vol][tot_sub_vol][3][3], double k_0, double complex alpha_0[tot_sub_vol]);

void get_G_old_struct_matrix_memory(int tot_sub_vol, double complex G_old[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char multithread);

void get_G_old_struct_matrix_memory_file(int tot_sub_vol, double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char multithread, char* G_old_file_name);

void get_G0_matrix_memory_2D(int tot_sub_vol, double complex G_0[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char multithread);

void get_alpha_prod_for_A_2D(int tot_sub_vol, double complex G_0[3*tot_sub_vol][3*tot_sub_vol], double k_0, double complex alpha_0[tot_sub_vol], double complex prod[3*tot_sub_vol][3*tot_sub_vol]);

void fill_A_matrix_2D(int tot_sub_vol, double complex A[3*tot_sub_vol][3*tot_sub_vol], double complex prod[3*tot_sub_vol][3*tot_sub_vol]);

#endif /* GREENSFUNCTION_H */
