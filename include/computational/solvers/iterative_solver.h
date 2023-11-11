#ifndef __iterative_solver_h
#define __iterative_solver_h

#include <stdlib.h>
#include <complex.h>
#include <math.h>

void matrix_reshape(int inner_size, int outer_size, double complex matrix_2d_1[][3*outer_size], double complex matrix_4d_1[outer_size][outer_size][inner_size][inner_size]);

void A2d_solver(double complex epsilon, int mm, int tot_sub_vol, double delta_V, double complex G_sys_old[][3* tot_sub_vol], double complex A_2d[][3*tot_sub_vol], double k);
void A2d_solver_memory(double complex epsilon, int mm, int tot_sub_vol, double delta_V, int size, double complex upperTriangularMatrix[size], double complex A_2d[][3*tot_sub_vol], double k);

void A_b_iterative_populator(int tot_sub_vol ,double complex *A_iterative, double complex *b_iterative, double complex A_2d[3][3], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], int mm, int jg_0);
void A_b_iterative_populator_memory(int tot_sub_vol, double complex *A_iterative, double complex *b_iterative, double complex A_2d[3][3], int size, double complex G_0_TriangularMatrix[size],int mm, int jg_0);

void G_sys_new_populator(int tot_sub_vol, int mm, int jg_0, double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], double complex b_iterative[]);
void G_sys_new_populator_memory(int tot_sub_vol, int mm, int jg_0, int size, double complex G_sys_TriangularMatrix[size], double complex b_iterative[]);

void offdiagonal_solver(int tot_sub_vol, int mm, double k, double complex alpha_0[], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], char multithread);
void offdiagonal_solver_memory(int tot_sub_vol, int mm, double k, double complex alpha_0[], int size, double complex G_0_TriangularMatrix[size], double complex G_sys_TriangularMatrix[size]);

void remaining_pertubations(int tot_sub_vol, int mm, double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], double complex A_2d[3][3], char multithread);
void remaining_pertubations_memory(int tot_sub_vol, int mm, int size, double complex G_0_TriangularMatrix[size], double complex G_sys_TriangularMatrix[size], double complex A_2d[3][3]);

void core_solver(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[],double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], char multithread);
void core_solver_memory(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[], int size, double complex G_sys_TriangularMatrix[size], double complex G_0_TriangularMatrix[size]);
void core_solver_store(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[],double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], char* G_sys_file_name);

void iterative_solver(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi,char multithread, double R[][3]);
//void iterative_solver(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol]);
//void iterative_solver(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi,double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3],char multithread);
void iterative_solver_memory(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[], int size, double complex G_sys_TriangularMatrix[size], double k_0, double pi,double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3],char wave_type);
void iterative_solver_store(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol],double k_0, double pi,double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3],char wave_type, char* G_sys_file_name);

void iterative_solver_file_handler(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[],double k_0, double pi,char multithread, char* G_old_file_name, char* G_sys_file_name, double R[][3]);
//void iterative_solver_file_handler(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[],double k_0, double pi,double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3],char multithread, char* G_old_file_name, char* G_sys_file_name);

void extractSubmatrix( double complex *largerMatrix, int largerRows, int largerCols, double complex *smallerMatrix, int smallerRows, int smallerCols, int startRow, int startCol);

void insertSubmatrix( double complex *largerMatrix, int largerRows, int largerCols, double complex *smallerMatrix, int smallerRows, int smallerCols, int startRow, int startCol);

#endif
