// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Read user inputs and create folder name in DSGF framework
// Developed by RETL Lab, Department of Mechanical Engineering, The University of Utah, USA.
// LAST UPDATE: JUNE 01, 2023  
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include <stdlib.h>
#include <sys/stat.h>
#include <stdio.h>
#include <time.h>
#include <file_utils.h>

void read_user_control(char *geometry,char *material, char *solution, char *single_spectrum_analysis, int *number_bulk_objects, int *number_omega, int *number_subvolumes_per_object, char* wave_type, char *multithread, double *epsilon_ref, char *uniform_spectra, char *save_spectral_conductance, char *save_total_conductance, char *save_power_dissipated_spectral_subvolumes,char *save_power_dissipated_total_subvolumes, char *save_power_dissipated_spectral_bulk, char *save_power_dissipated_total_bulk, char *save_power_density_total_subvolumes, char *save_spectral_transmissivity){

	char dirPathUserControl[260] = "user_inputs/control.txt";
	FILE *import_control_inputs = fopen(dirPathUserControl, "r"); 
	char buffer[256]; 

	fscanf(import_control_inputs,"%s = %s\n",buffer, geometry);	
	fscanf(import_control_inputs,"%s = %d\n", buffer, number_subvolumes_per_object);
	fscanf(import_control_inputs,"%s = %d\n", buffer, number_bulk_objects);
	fscanf(import_control_inputs,"%s = %s\n",buffer, material);
	fscanf(import_control_inputs,"%s = %lf\n", buffer, epsilon_ref);
	fscanf(import_control_inputs,"%s = %d\n", buffer, number_omega);
	fscanf(import_control_inputs,"%s = %c\n",buffer, uniform_spectra);
	fscanf(import_control_inputs,"%s = %c\n",buffer, solution);
	fscanf(import_control_inputs,"%s = %c\n",buffer, multithread);
	fscanf(import_control_inputs,"%s = %c\n",buffer, single_spectrum_analysis);
	fscanf(import_control_inputs,"%s = %c\n",buffer, wave_type);
	fscanf(import_control_inputs,"%s = %c\n",buffer, save_spectral_conductance);
	fscanf(import_control_inputs,"%s = %c\n",buffer, save_total_conductance);
	fscanf(import_control_inputs,"%s = %c\n",buffer, save_power_dissipated_spectral_subvolumes);
	fscanf(import_control_inputs,"%s = %c\n",buffer, save_power_dissipated_total_subvolumes);
	fscanf(import_control_inputs,"%s = %c\n",buffer, save_power_dissipated_spectral_bulk);
	fscanf(import_control_inputs,"%s = %c\n",buffer, save_power_dissipated_total_bulk);
	fscanf(import_control_inputs,"%s = %c\n",buffer, save_power_density_total_subvolumes);
	fscanf(import_control_inputs,"%s = %c\n",buffer, save_spectral_transmissivity); 
	fclose(import_control_inputs);
}

void read_geometry_sphere(double *d, double *radius, double *T1, double *T2){

	char dirPathGeometrySphereInputs[260] = "user_inputs/Geometry/sphere.txt";
	FILE *import_sphere_inputs = fopen(dirPathGeometrySphereInputs, "r"); 
	char buffer[256]; 
	fscanf(import_sphere_inputs,"%s = %lf\n",buffer, d);
	fscanf(import_sphere_inputs,"%s = %lf\n",buffer, radius);
	fscanf(import_sphere_inputs,"%s = %lf\n",buffer, T1);
	fscanf(import_sphere_inputs,"%s = %lf\n",buffer, T2);
	fclose(import_sphere_inputs);
}

void read_geometry_thin_films(double *d, double *Lx, double *Ly, double *Lz, double *T1, double *T2){

	char dirPathGeometryThinFilmsInputs[260] = "user_inputs/Geometry/thin_films.txt";
	FILE *import_thin_films_inputs = fopen(dirPathGeometryThinFilmsInputs, "r"); 
	char buffer[256]; 
	fscanf(import_thin_films_inputs,"%s = %lf\n",buffer, d);
	fscanf(import_thin_films_inputs,"%s = %lf\n",buffer, Lx);
	fscanf(import_thin_films_inputs,"%s = %lf\n",buffer, Ly);
	fscanf(import_thin_films_inputs,"%s = %lf\n",buffer, Lz);
	fscanf(import_thin_films_inputs,"%s = %lf\n",buffer, T1);
	fscanf(import_thin_films_inputs,"%s = %lf\n",buffer, T2);	
	fclose(import_thin_films_inputs);
}

void read_geometry_user_defined(double *d,char file_name_ud[], double *T1, double *T2){

	char dirPathGeometryUserDefinedInputs[260] = "user_inputs/Geometry/user_defined.txt";
	FILE *import_user_defined_inputs = fopen(dirPathGeometryUserDefinedInputs, "r"); 
	char buffer[256]; 
	fscanf(import_user_defined_inputs,"%s = %lf\n",buffer, d);
	fscanf(import_user_defined_inputs,"%s = %s\n",buffer, file_name_ud);
	fscanf(import_user_defined_inputs,"%s = %lf\n",buffer, T1);
	fscanf(import_user_defined_inputs,"%s = %lf\n",buffer, T2);	
	fclose(import_user_defined_inputs);
	
}

void read_calculation_temperatures(int N_Tcalc, double Tcalc_vector[]){

	char dirPathT_calc[] = "user_inputs/T_calc.txt";
	FILE *import_T_calc = fopen(dirPathT_calc, "r"); 
	for (int i_T_calc=0; i_T_calc<N_Tcalc;i_T_calc++) //tot_sub_vol
	{
		fscanf(import_T_calc,"%lf\n", &Tcalc_vector[i_T_calc]);
	}	
	fclose(import_T_calc);
}

void read_calculation_split(char *frequency_set){

	char dirPathCalculation_Split[260] = "user_inputs/calculation_split.txt";
	FILE *import_calculation_split = fopen(dirPathCalculation_Split, "r"); 
	char buffer[256]; 
	fscanf(import_calculation_split,"%s = %s",buffer, frequency_set);	
	fclose(import_calculation_split);
}

void create_folder(char folder_name[]){

	struct stat st = {0};

	if (stat(folder_name, &st) == -1){
		if (0 != mkdir(folder_name, 0700)){
			printf("when executing mkdir(\"%s\")\n", folder_name);
			perror("mkdir");
			//exit(1);
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

void write_to_csv_double_imag(char file_name[], int rows, int cols, double complex matrix[rows][cols]){
	FILE * CSV_file = fopen(file_name, "w");

	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			fprintf(CSV_file, "%e+%ej\t", creal(matrix[i][j]), cimag(matrix[i][j]));
		}
		fprintf(CSV_file, "\n");
	}

	fclose(CSV_file);

}

void write_to_csv_double_matrix(char file_name[], int rows, int cols, double matrix[rows][cols]){
	FILE * CSV_file = fopen(file_name, "w");
	
	for (int i = 0; i < rows; i++){
		for (int j = 0; j < cols; j++){
			// write the number
			fprintf(CSV_file, "%e", matrix[i][j]);

			// write a comma if its not the last number
			if (j < cols-1) fprintf(CSV_file, ",");
		}
		// start a new line
		fprintf(CSV_file, "\n");
	}

	fclose(CSV_file);
}

void write_to_csv_double_array(char file_name[], int length, double array[]){
	FILE * CSV_file = fopen(file_name, "w");

	for (int i = 0; i < length; i++){
		fprintf(CSV_file, "%e\n", array[i]);
	}

	fclose(CSV_file);

}

void populate_subvol_struct(char file_name[], int array_length, subvol shape[array_length]){
	
	int i_import = 0;
	FILE *import_discretization = fopen(file_name, "r");
	while (3 == fscanf(import_discretization, "%e %e %e", &shape[i_import].x, &shape[i_import].y, &shape[i_import].z))
	{   
		i_import++;
	}
	fclose(import_discretization);
}

void populate_subvol_delta_v(char file_name[], int array_length, double delta_V_vector[array_length]){
	FILE *import_delta_v = fopen(file_name, "r");
	for (int i = 0; i < array_length; i++) {
   		fscanf(import_delta_v,"%le", &delta_V_vector[i]);
 	}
	fclose(import_delta_v);
}

void create_pos_processing(char file_name[], char material[], double initial, double end, double time_spent, double Tcalc_vector[], double Total_conductance[], int N_Tcalc){
	
	FILE * pos_processing_summary; // call the main outputs' file,  

	pos_processing_summary =fopen(file_name,"w");

	fprintf(pos_processing_summary,"Material: %s\nSpectrum range (in wavelength) = %e--%e m \n",material, initial, end); 

	for (int i = 0; i < N_Tcalc; i++) 
	{
		fprintf(pos_processing_summary,"Total conductance at %eK = %e\n",Tcalc_vector[i], Total_conductance[i]); 
	}

	fprintf(pos_processing_summary,"Time counter: %f s\n",time_spent);
	
	fclose(pos_processing_summary);
}
