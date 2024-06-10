#ifndef __file_utils_h
#define __file_utils_h

#include <complex.h>
#include <functions_DSGF.h>

void read_user_control(char *geometry,char *material, char *solution, char *single_spectrum_analysis, int *number_subvolumes_per_object_1, int *number_subvolumes_per_object_2, int *number_omega, char *multithread, double *epsilon_ref, char *frequency_set, char *save_spectral_conductance, char *save_total_conductance, char *save_power_dissipated_spectral_subvolumes,char *save_power_dissipated_total_subvolumes, char *save_power_dissipated_spectral_bulk, char *save_power_dissipated_total_bulk, char *save_power_density_total_subvolumes, char *save_spectral_transmissivity);
//void read_user_control(char *geometry,char *material, char *solution, char *single_spectrum_analysis, int *number_subvolumes_per_object_1, int *number_subvolumes_per_object_2, int *number_omega, char *multithread, double *epsilon_ref, char *uniform_spectra, char *save_spectral_conductance, char *save_total_conductance, char *save_power_dissipated_spectral_subvolumes,char *save_power_dissipated_total_subvolumes, char *save_power_dissipated_spectral_bulk, char *save_power_dissipated_total_bulk, char *save_power_density_total_subvolumes, char *save_spectral_transmissivity);

void read_geometry_sample(char *geometry_1, double *radius1, char *geometry_2, double *radius2, double *d, double *T1, double *T2);

void read_geometry_sphere(double *d, double *radius, double *T1, double *T2);

void read_geometry_thin_films(double *d, double *Lx, double *Ly, double *Lz, double *T1, double *T2);

void read_geometry_user_defined(double *d, char file_name_ud[], double *T1, double *T2);

void read_calculation_temperatures(int N_Tcalc, double Tcalc_vector[]);

void create_folder(char folder_name[]);

// this returns the final folder for all the variables and where the matrices will be stored
char* set_up_results(char material[], char geometry[], int tot_sub_vol, double d);

void write_to_csv_double_imag(char file_name[], int rows, int cols, double complex matrix[rows][cols]); 

void write_to_csv_double_matrix(char file_name[], int rows, int cols, double matrix[rows][cols]);

void write_to_csv_double_array(char file_name[], int length, double array[]);

void populate_subvol_struct(char file_name[], int array_length, subvol shape[array_length]);

void populate_subvol_struct_object_2(char file_name[], int array_length, subvol shape[array_length], int subvol_per_object);

void populate_subvol_delta_v(char file_name[], int array_length, double delta_V_vector[array_length]);

void create_pos_processing(char file_name[], char material[], double initial, double end, double time_spent, double T_calc_vector[], double Total_conductance[], int N_Tcalc);

void write_bin(int tot_sub_vol, double complex G_array[][3*tot_sub_vol], char* file_name);

void write_csv(int tot_sub_vol, double complex G_array[][3*tot_sub_vol], char* file_name);

void read_bin(int tot_sub_vol, double complex G_array[][3*tot_sub_vol], char* file_name);

void read_csv(int tot_sub_vol, double complex G_array[][3*tot_sub_vol], char* file_name);

void read_populator_bin(int tot_sub_vol, double complex G_array[3*tot_sub_vol*3*tot_sub_vol], char* file_name);



#endif
