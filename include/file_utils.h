#ifndef __file_utils_h
#define __file_utils_h

#include <complex.h>
#include <functions_DSGF.h>

void read_user_control(char *geometry,char *material, char *solution, char *single_spectrum_analysis, char *save_spectral_conductance, char *save_spectral_transmissivity, char *save_power_dissipated, int *number_bulk_objects, int *number_omega, int *number_subvolumes_per_object, char* wave_type, char *non_uniform_subvolumes);

void read_geometry_sphere(double *d, double *radius, double *T1, double *T2);

void read_geometry_thin_films(double *d, double *Lx, double *Ly, double *Lz, double *T1, double *T2);

void read_calculation_temperatures(int N_Tcalc, double Tcalc_vector[]);

void create_folder(char folder_name[]);

// this returns the final folder for all the variables and where the matrices will be stored
char* set_up_results(char material[], char geometry[], int tot_sub_vol, double d);

void write_to_csv_double_imag(char file_name[], int rows, int cols, double complex matrix[rows][cols]); 

void write_to_csv_double_matrix(char file_name[], int rows, int cols, double matrix[rows][cols]);

void write_to_csv_double_array(char file_name[], int length, double array[]);

void populate_subvol_struct(char file_name[], int array_length, subvol shape[array_length]);

void create_pos_processing(char file_name[], char material[], double initial, double end, double time_spent, double T_calc_vector[], double Total_conductance[], int N_Tcalc);

#endif
