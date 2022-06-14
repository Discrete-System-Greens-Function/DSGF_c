#ifndef __file_utils_h
#define __file_utils_h

void read_user_inputs(char *material, char *geometry, char *discretization_thin_film, double *d, double *radius, double *Lx, double *Ly, double *Lz, double *T1, double *T2, char *solution, char *single_spectrum_analysis, char *save_A_matrix, char *save_G0_matrix, char *save_SGF_matrix, char *save_spectral_conductance, char *save_spectral_transmissivity, char *save_power_dissipated);

void read_calculation_temperatures(int N_Tcalc, double Tcalc_vector[]);

int read_int_from_file(char file_name[]);

#endif
