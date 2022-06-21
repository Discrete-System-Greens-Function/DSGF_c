// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Read user inputs and create folder name in DSGF framework
// Developed by RETL Lab, Department of Mechanical Engineering, The University of Utah, USA.
// LAST UPDATE: JUNE 20, 2022  
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "file_utils.h"

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

void create_folder(char folder_name[]){
	
	struct stat st = {0};

	if (stat(folder_name, &st) == -1){
		if (0 != mkdir(folder_name, 0700)){
			printf("when executing mkdir(\"%s\")\n", folder_name);
			perror("mkdir");
			exit(1);
		}
	}
}

void get_time_date_string(char *time_date){
	time_t rawtime;
	struct tm *timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(time_date, 80, "%y_%m_%d_%H_%M_%S", timeinfo);
}

char* set_up_results(char material[], char geometry[], int tot_sub_vol, double d){

	char parent_folder[] = "results";

	create_folder(parent_folder);	

	char geometry_folder[256];
	sprintf(geometry_folder, "%s/%s_%d_subvolumes", parent_folder, geometry, tot_sub_vol);

	create_folder(geometry_folder);



	char material_folder[256];
	sprintf(material_folder,"%s/%s_d_%.2e", geometry_folder, material, d);

	create_folder(material_folder);

	char time_date[80];
	get_time_date_string(time_date);
	

	char *results_folder = malloc(256*sizeof(char));
	sprintf(results_folder, "%s/%s", material_folder, time_date);

	create_folder(results_folder);
	

	return results_folder;
}
