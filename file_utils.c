#include "file_utils.h"
#include "stdio.h"

void read_user_inputs(char *material, char *geometry, char *discretization_thin_film, double *d, double *radius, double *Lx, double *Ly, double *Lz, double *T1, double *T2, char *solution, char *single_spectrum_analysis, char *save_A_matrix, char *save_G0_matrix, char *save_SGF_matrix, char *save_spectral_conductance, char *save_spectral_transmissivity, char *save_power_dissipated){
	
    char dirPathUserInputs[260] = "user_inputs.txt";
    FILE *import_inputs = fopen(dirPathUserInputs, "r"); 
    char buffer[256]; 
    fscanf(import_inputs,"%s = %s\n",buffer, material);
    fscanf(import_inputs,"%s = %s\n",buffer, geometry);
    fscanf(import_inputs,"%s = %s\n",buffer, discretization_thin_film);
    fscanf(import_inputs,"%s = %lf\n",buffer, d);
    fscanf(import_inputs,"%s = %lf\n",buffer, radius);
    fscanf(import_inputs,"%s = %lf\n",buffer, Lx);
    fscanf(import_inputs,"%s = %lf\n",buffer, Ly);
    fscanf(import_inputs,"%s = %lf\n",buffer, Lz);
    fscanf(import_inputs,"%s = %lf\n",buffer, T1);
    fscanf(import_inputs,"%s = %lf\n",buffer, T2);
    fscanf(import_inputs,"%s = %c\n",buffer, solution);
    fscanf(import_inputs,"%s = %c\n",buffer, single_spectrum_analysis);
    fscanf(import_inputs,"%s = %c\n",buffer, save_A_matrix);
    fscanf(import_inputs,"%s = %c\n",buffer, save_G0_matrix);
    fscanf(import_inputs,"%s = %c\n",buffer, save_SGF_matrix);
    fscanf(import_inputs,"%s = %c\n",buffer, save_spectral_conductance);
    fscanf(import_inputs,"%s = %c",buffer, save_spectral_transmissivity); 
    fscanf(import_inputs,"%s = %c\n",buffer, save_power_dissipated);	
    fclose(import_inputs);
}

void read_calculation_temperatures(int N_Tcalc, double Tcalc_vector[]){

    char dirPathT_calc[] = "T_calc.txt";
    FILE *import_T_calc = fopen(dirPathT_calc, "r"); 
    for (int i_T_calc=0; i_T_calc<N_Tcalc;i_T_calc++) //tot_sub_vol
    {
	    fscanf(import_T_calc,"%lf\n", &Tcalc_vector[i_T_calc]);
    }	
    fclose(import_T_calc);
}

int read_int_from_file(char *file_name){

	FILE *int_file = fopen(file_name, "r"); 

	int temp;

	fscanf(int_file,"%d", &temp);
	fclose(int_file);

	const int const_int = temp;
	return const_int;
}
