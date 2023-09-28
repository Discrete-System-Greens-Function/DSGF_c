// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Discrete System Green's Function
// Near-field radiative heat transfer framework between thermal objects
// Developed by RETL group at the University of Utah, USA

// LAST UPDATE: September 27, 2023
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// General c libraries
#include <stdio.h>
#include <math.h>

#include <complex.h> // complex numbers library, code must be compiled as c99 standard https://stackoverflow.com/questions/6418807/how-to-work-with-complex-numbers-in-c
#include <stdlib.h>	 // export/import data

// Libraries for creating directories and files using a loop in C. Sources: https://stackoverflow.com/questions/46612504/creating-directories-and-files-using-a-loop-in-c and https://stackoverflow.com/questions/7430248/creating-a-new-directory-in-c
#include <string.h> // library used to concatenate 2 strings https://stackoverflow.com/questions/46612504/creating-directories-and-files-using-a-loop-in-c

#include <errno.h>
// Library with the inputs and functions for DSGF
#include "user_inputs.h"	// User inputs definitions header. No values are defined in this file.
#include "functions_DSGF.h" // Definitions of functions used in DSGF
#include "file_utils.h"		// header with definitions of read_user_inputs and read_calculation_temperatures functions
#include "computational/solvers/iterative_solver.h"
#include "array_functions.h"
#include "memory.h"

#include <malloc.h>

#include "time.h"

#include "geometry/sample.h"
//#include "geometry/sphere.h"
//#include "geometry/thin_film.h"
#include "geometry/shared.h"
#include "geometry/user_defined.h"

#include "material/SiO2.h"
#include "material/SiC.h"
#include "material/u_SiC.h"
#include "material/SiN.h"
#include "computational/GreensFunction.h"
#include "computational/solvers/direct_solver.h"
#include "computational/ThermalPower.h"

// LAPACKE libraries: https://www.netlib.org/lapack/lapacke.html ; https://extras.csc.fi/math/nag/mark21/pdf/F08/f08anf.pdf
#include <lapacke.h>
#include "lapack_header.h" //header with Lapack definitions

// #################################################################
// ################### START OF THE CODE ###########################
// #################################################################

int main()
{

	// ----------- CONSTANTS ------------------

	const double pi = 3.14159265359;			 // pi number
	const double q = 1.602176634e-19;			 // number of Joules per eV
	const double h_bar = 1.054571817e-34;		 // Planck's constant [J*s]
	const double k_b = 1.38064852e-23;			 // Boltzmann constant [J/K]
	const double epsilon_0 = 8.8542e-12;		 // Permittivity of free space [F/m]
	const double c_0 = 299792458;				 // Speed of light in vacuum [m/s]
	const double mu_0 = (4. * pi) * pow(10, -7); // Permeability of free space [H/m]
	
	baseline = get_mem_usage(); // measure baseline memory usage
	
	time_t simulation_start;
	time(&simulation_start);

	int N_subvolumes_per_object, N_bulk_objects, N_omega;

	char wave_type, multithread;

	read_user_control(geometry, material, &solution, &single_spectrum_analysis, &N_subvolumes_per_object_1, &N_subvolumes_per_object_2, &N_omega, &multithread, &epsilon_ref, &uniform_spectra, &save_spectral_conductance, &save_total_conductance, &save_power_dissipated_spectral_subvolumes, &save_power_dissipated_total_subvolumes, &save_power_dissipated_spectral_bulk, &save_power_dissipated_total_bulk, &save_power_density_total_subvolumes, &save_spectral_transmissivity); //&wave_type,
	read_calculation_temperatures(N_Tcalc, Tcalc_vector);
	
	char frequency_set[260]; // definition for the file with the frequencies 
	read_calculation_split(frequency_set); 
	
	int const const_N_bulk_objects = N_bulk_objects; // Number of objects
	int const const_N_omega = N_omega; // Number of frequencies to be computed 

	int const const_N_subvolumes_per_object = N_subvolumes_per_object_1; // Number of subvolumes per object
	int const const_N_subvolumes_per_object_2 = N_subvolumes_per_object_2;
	int tot_sub_vol = const_N_subvolumes_per_object + const_N_subvolumes_per_object_2; // Assign tot_sub_vol: Computes the total number of subvolumes in the system. tot_sub_vol is defined this way because it must be a non-variable parameter due to the computations perfomed in the code. Previously, it was defined as #define tot_sub_vol const_N_subvolumes_per_object*const_N_bulk_objects
	
	/*
	double(*G_12_omega_SGF)[N_Tcalc] = calloc(const_N_omega, sizeof(*G_12_omega_SGF)); // spectral conductance
	if (G_12_omega_SGF == NULL){
			printf("Failure with memory before spectral analysis when conductance defined. ");
			return 1;
	}
	*/
	
	// ######### Properties for thermal objects ###########
	printf("%s simulation for a total of %d subvolumes \n", geometry, tot_sub_vol);
	
	double(*omega) = malloc(sizeof *omega * const_N_omega); // radial frequency [rad/s]
	if (omega == NULL){
			printf("Failure with memory before spectral analysis when frequencies are defined. ");
			return 1;
	}
	
	double delta_V_1, delta_V_2;
	double(*delta_V_vector) = malloc(sizeof *delta_V_vector * tot_sub_vol); // Vector of all subvolume size. Combines delta_V_1 and delta_V_2 in the same array
	if (delta_V_vector == NULL){
			printf("Failure with memory before spectral analysis when the size of subvolumes are defined. ");
			return 1;
	}

	double(*R)[3] = calloc(tot_sub_vol, sizeof(*R)); // center of subvolumes for thermal objects: info imported from a .txt file
	if (R == NULL){
			printf("Failure with memory before spectral analysis when the positions of subvolumes are defined. ");
			return 1;
	}	
	
	if (strcmp(geometry, "sample") == 0)
	{ 
		// R is populated here!
		set_up_sample_geometry(pi, tot_sub_vol, const_N_subvolumes_per_object, const_N_subvolumes_per_object_2, &T1, &T2, &d, &delta_V_1, &delta_V_2, R, geometry_1, geometry_2);
	}
	/*
	if (strcmp(geometry, "sphere") == 0)
	{
		set_up_sphere_geometry(pi, tot_sub_vol, const_N_subvolumes_per_object, &T1, &T2, &d, &delta_V_1, &delta_V_2, R);
	}
	if (strcmp(geometry, "thin-films") == 0)
	{
		set_up_thin_film_geometry(tot_sub_vol, const_N_subvolumes_per_object, const_N_bulk_objects, &T1, &T2, &d, &delta_V_1, &delta_V_2, R);
	}
	*/
	if(strcmp(geometry,"user_defined")==0)
	{
		//char file_name;
		set_up_user_defined_geometry(tot_sub_vol,const_N_subvolumes_per_object,const_N_bulk_objects, &d, &T1, &T2, R, delta_V_vector);//T_vector,delta_V_vector, file_name_ud
	}
	else{ set_delta_V_vector(delta_V_1, delta_V_2, tot_sub_vol, const_N_subvolumes_per_object, delta_V_vector);}
	
	char *results_folder = set_up_results(material, geometry, tot_sub_vol, d); // Folders for results 

	char copy_control[260];
	sprintf(copy_control, "cp ./user_inputs/control.txt  ./%s\n", results_folder);
	system(copy_control);
		
	char copy_geometry[260];
	sprintf(copy_geometry, "cp ./user_inputs/Geometry/%s.txt  ./%s\n",geometry,results_folder);
	system(copy_geometry);


	if ( save_power_dissipated_spectral_subvolumes == 'Y' ||
 		 save_power_dissipated_total_subvolumes == 'Y' ||
 		 save_power_dissipated_spectral_bulk == 'Y' ||
 		 save_power_dissipated_total_bulk == 'Y' ||
 		 save_power_density_total_subvolumes == 'Y' )
	{
		// EXPORT R
		char dirPathVector_subvolumes_lattice_FileName[260];
		sprintf(dirPathVector_subvolumes_lattice_FileName, "%s/vector_subvolumes_lattice.csv", results_folder); // path where the file is stored
		write_to_csv_double_matrix(dirPathVector_subvolumes_lattice_FileName, tot_sub_vol, 3, R);

		// EXPORT delta_V_vector
		char dirPathVector_subvolumes_volume_FileName[260];
		sprintf(dirPathVector_subvolumes_volume_FileName, "%s/vector_subvolumes_volume.csv", results_folder); // path where the file is stored
		write_to_csv_double_array(dirPathVector_subvolumes_volume_FileName, tot_sub_vol, delta_V_vector);
	} // end if save_power_dissipated

	
	
	
	if(strcmp(material,"SiO2")==0 || strcmp(material,"u-SiC")==0)  // removed strcmp(material,"SiC")==0 || from uniform
	{
		if (uniform_spectra == 'Y')
		{
			double initial,final;
			initial = 5.e-6;
			final = 25.e-6;
			double lambda[const_N_omega];
			double_linspace(initial, final, const_N_omega, lambda);
			for (int i_lambda = 0; i_lambda < const_N_omega; i_lambda++)
			{
				omega[i_lambda] = 2. * pi * c_0 / lambda[i_lambda]; // Radial frequency [rad/s]
			}
		}
		if (uniform_spectra == 'N')
		{
			FILE *non_uniform_spectra; // Import non-uniform spectra
			char dirPathFileNameSpectra[260];

			sprintf(dirPathFileNameSpectra, "library/Non_uniform_spectra/SiO2_non_uniform_spectra_%d.csv", const_N_omega);
			non_uniform_spectra = fopen(dirPathFileNameSpectra, "r");
			for (int i = 0; i < const_N_omega; i++)
			{
				fscanf(non_uniform_spectra, "%lf", &omega[i]); //
			}
			fclose(non_uniform_spectra);
		}	
		
		if (uniform_spectra == 'S')
		{
			FILE *spectra_split; // Import non-uniform spectra
			char dirPathFileNameSpectraSplit[260];

			sprintf(dirPathFileNameSpectraSplit, "library/Non_uniform_spectra/%s.csv", frequency_set);
			spectra_split = fopen(dirPathFileNameSpectraSplit, "r");
			for (int i = 0; i < const_N_omega; i++)
			{
				fscanf(spectra_split, "%lf", &omega[i]); //
			}
			fclose(spectra_split);
		}
		
	}
	else if (strcmp(material, "SiC") == 0)
	{
		if (uniform_spectra == 'Y')
		{
			double initial,final;
			initial = 1.4e14;
			final = 1.9e14;
			double_linspace(initial, final, const_N_omega, omega);
		}
		if (uniform_spectra == 'N')
		{
			FILE *non_uniform_spectra; 	// Import non-uniform spectra
			char dirPathFileNameSpectra[260];

			sprintf(dirPathFileNameSpectra, "library/Non_uniform_spectra/SiC_non_uniform_spectra_%d.csv", const_N_omega);
			non_uniform_spectra = fopen(dirPathFileNameSpectra, "r");
			for (int i = 0; i < const_N_omega; i++)
			{
				fscanf(non_uniform_spectra, "%lf", &omega[i]); //
			}
			fclose(non_uniform_spectra);
		}	
	}
	else if (strcmp(material, "SiN") == 0) // cannot compare strings in C with ==; source: https://ideone.com/BrFA00
	{
		if (uniform_spectra == 'Y')
		{
		double initial,final;
		initial = 2.e13;
		final = 3.e14;
		double_linspace(initial, final, const_N_omega, omega);
		}
		if (uniform_spectra == 'N')
		{
			FILE *non_uniform_spectra; // Import non-uniform spectra
			char dirPathFileNameSpectra[260];

			sprintf(dirPathFileNameSpectra, "library/Non_uniform_spectra/SiN_non_uniform_spectra_%d.csv", const_N_omega);
			non_uniform_spectra = fopen(dirPathFileNameSpectra, "r");
			for (int i = 0; i < const_N_omega; i++)
			{
				fscanf(non_uniform_spectra, "%lf", &omega[i]); //
			}
			fclose(non_uniform_spectra);
		}
	}

	// ################## FREE-SPACE GREEN'S FUNCTION AND INTERACTION A MATRIX #####################
	//  Fill terms for G^0:		

	double(*modulo_r_i_j)[tot_sub_vol] = malloc(sizeof *modulo_r_i_j * tot_sub_vol);
	if (modulo_r_i_j == NULL){
			printf("Failure with memory before spectral analysis when the distances between subvolumes are defined. ");
			return 1;
	}
	double complex(*r_i_j_outer_r_i_j)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*r_i_j_outer_r_i_j));
	if (r_i_j_outer_r_i_j == NULL){
			printf("Failure with memory before spectral analysis when the distances between subvolumes are defined. ");
			return 1;
	} 
	setup_G_0_matrices(tot_sub_vol, modulo_r_i_j, r_i_j_outer_r_i_j, R);
	free(R);
	// #################################################################
	// ################## FREQUENCY RANGE ANALYSIS #####################
	// #################################################################
	// Loop to analyze a range of desired frequencies
	//printf("----- Spectrum range calculation -----\n");

	if(solution =='D'){printf("spectrum range calculation using direct inversion: \n");}
	if(solution =='I'){printf("spectrum range calculation using iterative solver: \n");}

	int omega_range;
	if (single_spectrum_analysis == 'Y')
		omega_range = 1;
	if (single_spectrum_analysis == 'N')
		omega_range = const_N_omega;

	#pragma omp parallel for if (multithread == 'Y')

	for (int i_omega = 0; i_omega < omega_range; i_omega++) // Frequency loop
	{

		if (i_omega == 1)
		{
			single_analysis = 'n';
		}

		double omega_value = omega[i_omega]; // omega definition is on line 212
		//printf("%d) omega = %e. ", i_omega + 1, omega_value);
		printf("%d)\n", i_omega + 1);
		double complex epsilon;
		if (strcmp(material, "SiO2") == 0) // cannot compare strings in C with ==; source: https://ideone.com/BrFA00
		{
			epsilon = calculate_epsilon_SiO2(q, omega_value, h_bar);
		}

		else if (strcmp(material, "SiC") == 0)
		{
			epsilon = calculate_epsilon_SiC(omega_value);
		}
		else if (strcmp(material, "u-SiC") == 0)
		{
			epsilon = calculate_epsilon_u_SiC(omega_value, pi, c_0);
		}

		else if (strcmp(material, "SiN") == 0)
		{
			epsilon = calculate_epsilon_SiN(omega_value, pi);
		}

		double k_0 = k_0_function(omega_value, epsilon_0, mu_0); // wave vector in free space

		double k = k_function(omega_value, epsilon_ref, epsilon_0, mu_0); // wave vector in reference medium

		//double complex alpha_0[tot_sub_vol];
		double complex (*alpha_0) = malloc(sizeof *alpha_0 * tot_sub_vol); 
		if (alpha_0 == NULL){
			printf("Failure with memory when alpha_0 is defined. ");
			exit(1);
		}
		for (int i_alpha = 0; i_alpha < tot_sub_vol; i_alpha++)
		{
			alpha_0[i_alpha] = delta_V_vector[i_alpha]*(epsilon - epsilon_ref); //Bare polarizability [m^3]
			//printf("%e +i%e\n ",creal(alpha_0[i_alpha]),cimag(alpha_0[i_alpha]));
			//printf("%e\n",delta_V_vector[i_alpha]);

		}
				
		// ######### Calculation of SGF and post-processing ###########

		double complex trans_coeff;
		double inner_sum;
		double complex G_element;
		double Q_omega_subvol;
		double sum_trans_coeff = 0; // Transmissivity to calculate spectral conductance and spectral power in subvolumes

    	pre_solver_memory = get_mem_usage()-baseline;
		if(solution =='D') // Solves the linear system AG=G^0, where G^0 and A are 3N X 3N matrices. 
		{	
			double complex (*G_0)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*G_0)); 	
			if (G_0 == NULL){
				printf("Failure with memory. Use iterative solver");
				exit(1);
			} 
			//get_G0_matrix(tot_sub_vol, G_0, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, wave_type);
			get_G0_matrix_memory(tot_sub_vol, G_0, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, wave_type);

			double complex (*A)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*A));
			if (A == NULL){
				printf("Failure with memory. Use iterative solver");
				exit(1);
			}
			get_A_matrix(tot_sub_vol, G_0, A, k_0, alpha_0); // function applicable for uniform and non-uniform discretization

			double complex (*A_direct) = calloc(3*3*tot_sub_vol*tot_sub_vol, sizeof(*A_direct));
			if (A_direct == NULL){
				printf("Failure with memory. Use iterative solver");
				exit(1);
			}
			A_direct_populator(tot_sub_vol, A, A_direct);
			free(A);
			
			double complex (*b_direct) = calloc(3*3*tot_sub_vol*tot_sub_vol, sizeof(*b_direct));
			if (b_direct == NULL){
				printf("Failure with memory. Use iterative solver");
				exit(1);
			}
			b_direct_populator(tot_sub_vol, G_0, b_direct);
			free(G_0);

			direct_solver_memory(tot_sub_vol, A_direct, b_direct); // we noticed an memory improved from the old function
			
			calculation_memory = get_mem_usage()-pre_solver_memory-baseline;
			free(A_direct);

			double complex (*G_sys)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys));
			if (G_sys == NULL){
			printf("Failure with memory. Use iterative solver");
			exit(1);
			} 
			populate_G_sys(tot_sub_vol, b_direct, G_sys);
			free(b_direct);

			
			//direct_solver(tot_sub_vol, A_direct, b_direct, G_sys);  
			//direct_solver(tot_sub_vol, A, G_0, G_sys);
			//direct_solver(tot_sub_vol, G_sys, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, wave_type, alpha_0); // we noticed an memory improved from the old function
			//spectral_post_processing(tot_sub_vol, i_omega, const_N_omega, k_0, h_bar, k_b, epsilon, omega_value, T_vector, delta_V_vector, const_N_subvolumes_per_object, pi, G_sys, &sum_trans_coeff, Q_subvol);

			for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
			{
				Q_omega_subvol = 0;
				for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
				{
					G_element = 0;
					for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
					{
						int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
						for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
						{
							int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
							// Extract one Green's function component: G_sys[ig_0][jg_0][i_subG_0]
							double complex transpose_G_sys = G_sys[ig_0_2d][jg_0_2d];
							double complex G_sys_cross = conj(transpose_G_sys);
							G_element+=G_sys[ig_0_2d][jg_0_2d]*G_sys_cross; 
						}    
					}
					// Transmissivity coefficient matrix tau(omega) for comparison with Czapla Mathematica output [dimensionless]
					trans_coeff = 4.*pow(k_0,4)*delta_V_vector[ig_0]*delta_V_vector[jg_0]*cimag(epsilon)*cimag(epsilon)*G_element; 
				
					//Trans_bulk: Transmission coefficient between two bulk objects
					// This function calculates the transmission coefficient between bulk objects given the transmission coefficient between every pair of dipoles for a given frequency.
					if(ig_0 < const_N_subvolumes_per_object && jg_0 >= const_N_subvolumes_per_object)// bulk 1 -> 2
					{
						sum_trans_coeff += trans_coeff;
					} 
					// Thermal power dissipated calculation, based on the matlab code (Using Tervo's Eq. 26)
					if (ig_0 != jg_0) 
					{
						//inner_sum = (theta_function(omega_value, T_vector[jg_0], h_bar, k_b) - theta_function(omega_value, T_vector[ig_0], h_bar, k_b)) * trans_coeff;
						double T_i,T_j;
						if (ig_0<const_N_subvolumes_per_object){ T_i=T1;}
						else{ T_i=T2;}
						if (jg_0<const_N_subvolumes_per_object){ T_j=T1;}
						else{ T_j=T2;}
						inner_sum = (theta_function(omega_value, T_j, h_bar, k_b) - theta_function(omega_value, T_i, h_bar, k_b)) * trans_coeff;
					}
					else {inner_sum = 0; }
					Q_omega_subvol += (1 / (2 * pi)) * inner_sum; // calculates the thermal power dissipated per subvolume
				}
				if ( save_power_dissipated_spectral_subvolumes == 'Y' ||
 		 		save_power_dissipated_total_subvolumes == 'Y' ||
 				save_power_dissipated_spectral_bulk == 'Y' ||
 		 		save_power_dissipated_total_bulk == 'Y' ||
 		 		save_power_density_total_subvolumes == 'Y' )
				{
					FILE *power_dissipated_spectral; // append
					char dirPathPower_dissipated_spectral_subvolumes_FileName[260];
					sprintf(dirPathPower_dissipated_spectral_subvolumes_FileName, "%s/Q_omega_subvol.csv", results_folder); // path where the file is stored
					if (ig_0 == 0 && i_omega ==0){ power_dissipated_spectral = fopen(dirPathPower_dissipated_spectral_subvolumes_FileName, "w");} // write
					else if (ig_0 == 0 && i_omega!=0)
					{
						power_dissipated_spectral = fopen(dirPathPower_dissipated_spectral_subvolumes_FileName, "a"); // append
						fprintf(power_dissipated_spectral, "\n");
					}
					else { power_dissipated_spectral = fopen(dirPathPower_dissipated_spectral_subvolumes_FileName, "a"); }// append
					fprintf(power_dissipated_spectral, "%e,", Q_omega_subvol);					
					fclose(power_dissipated_spectral);
				
				}
				/*
				if ( save_power_dissipated_spectral_subvolumes == 'Y' ||
 		 		save_power_dissipated_total_subvolumes == 'Y' ||
 				save_power_dissipated_spectral_bulk == 'Y' ||
 		 		save_power_dissipated_total_bulk == 'Y' ||
 		 		save_power_density_total_subvolumes == 'Y' )
				{ Q_subvol[ig_0][i_omega] = Q_omega_subvol; }
				*/
			}	
			free(G_sys);
		} // end if (solution =='D')

		if(solution =='I')
		{ 
			//printf("Iterative solver in progress. \n");

			//full matrix definition
			double complex (*G_sys)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys));
			    	
			/*
			//triangular matrix definition
			int size = 9*tot_sub_vol*(tot_sub_vol+1)/2; // Calculate the size of the 1D array to store the upper triangular matrix, via chatgpt
			double complex* G_sys_TriangularMatrix = (double complex*)malloc(size * sizeof(double complex)); // Allocate memory for the 1D array
			*/

			//old versions
			//iterative_solver(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_sys_old, G_sys); //  we noticed an memory improved from the old function iterative_solver(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_0, G_sys);
			//iterative_solver(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_sys,k_0, pi,modulo_r_i_j, r_i_j_outer_r_i_j,wave_type);
			
			//full matrix solution
			iterative_solver(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_sys,k_0, pi,modulo_r_i_j, r_i_j_outer_r_i_j,wave_type);
			
			/*
			char G_old_file_name[260];
			char G_sys_file_name[260];
			if (multithread =='Y')
			{
				sprintf(G_old_file_name,"%s/G_old_%d.bin",results_folder,i_omega); //
				//sprintf(G_old_file_name,"%s/G_old_%d.csv",results_folder,i_omega);
				sprintf(G_sys_file_name,"%s/G_sys_%d.bin",results_folder,i_omega); //
				//sprintf(G_sys_file_name,"%s/G_sys_%d.csv",results_folder,i_omega);
			}
			else
			{
				sprintf(G_old_file_name,"%s/G_old.bin",results_folder); //
				//sprintf(G_old_file_name,"%s/G_old.csv",results_folder);
				sprintf(G_sys_file_name,"%s/G_sys.bin",results_folder); //
				//sprintf(G_sys_file_name,"%s/G_sys.csv",results_folder);
			}

			//iterative_solver_store(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_sys,k_0, pi,modulo_r_i_j, r_i_j_outer_r_i_j,wave_type,G_sys_file_name);
			iterative_solver_with_data(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0,k_0, pi,modulo_r_i_j, r_i_j_outer_r_i_j,wave_type,G_old_file_name, G_sys_file_name);
						
			double complex (*G_sys)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys));	
			read_bin(tot_sub_vol, G_sys, G_sys_file_name);
			*/
			
			//triangular matrix solution
			//iterative_solver_memory(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, size, G_sys_TriangularMatrix,k_0, pi,modulo_r_i_j, r_i_j_outer_r_i_j,wave_type);
			
			//spectral_post_processing(tot_sub_vol, i_omega, const_N_omega, k_0, h_bar, k_b, epsilon, omega_value, T_vector, delta_V_vector, const_N_subvolumes_per_object, pi, G_sys, &sum_trans_coeff, Q_subvol);
			
			for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
			{
			Q_omega_subvol = 0;
				for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) // original
				//for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++) // mod: Sept.6,2023
				{
					G_element = 0;
					for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
					{
						int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
						for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
						{
							int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
							int index = 9 * (ig_0 * tot_sub_vol + jg_0) + 3 * i_subG_0 + j_subG_0;
							// Extract one Green's function component: G_sys[ig_0][jg_0][i_subG_0]
							double complex transpose_G_sys = G_sys[ig_0_2d][jg_0_2d]; // full matrix
							//double complex transpose_G_sys = G_sys_TriangularMatrix[index]; // triangular matrix
							double complex G_sys_cross = conj(transpose_G_sys);
							G_element+=G_sys[ig_0_2d][jg_0_2d]*G_sys_cross; // full matrix
							//G_element+=G_sys_TriangularMatrix[index]*G_sys_cross; // triangular matrix
						}    
					}
					// Transmissivity coefficient matrix tau(omega) for comparison with Czapla Mathematica output [dimensionless]
					trans_coeff = 4.*pow(k_0,4)*delta_V_vector[ig_0]*delta_V_vector[jg_0]*cimag(epsilon)*cimag(epsilon)*G_element; 
			
					//Trans_bulk: Transmission coefficient between two bulk objects
					// This function calculates the transmission coefficient between bulk objects given the transmission coefficient between every pair of dipoles for a given frequency.
					if(ig_0 < const_N_subvolumes_per_object && jg_0 >= const_N_subvolumes_per_object)// bulk 1 -> 2
					{
						sum_trans_coeff += trans_coeff;
					} 
					// Thermal power dissipated calculation, based on the matlab code (Using Tervo's Eq. 26)
					if (ig_0 != jg_0){ 
						//inner_sum = (theta_function(omega_value, T_vector[jg_0], h_bar, k_b) - theta_function(omega_value, T_vector[ig_0], h_bar, k_b)) * trans_coeff; 
						double T_i,T_j;
						if (ig_0<const_N_subvolumes_per_object){ T_i=T1;}
						else{ T_i=T2;}
						if (jg_0<const_N_subvolumes_per_object){ T_j=T1;}
						else{ T_j=T2;}
						//read file to use T_i and T_j;  !!!!!!!
						inner_sum = (theta_function(omega_value, T_j, h_bar, k_b) - theta_function(omega_value, T_i, h_bar, k_b)) * trans_coeff;
					}
					else { inner_sum = 0; }
					Q_omega_subvol += (1 / (2 * pi)) * inner_sum; // calculates the thermal power dissipated per subvolume
				}
				if ( save_power_dissipated_spectral_subvolumes == 'Y' ||
 		 		save_power_dissipated_total_subvolumes == 'Y' ||
 				save_power_dissipated_spectral_bulk == 'Y' ||
 		 		save_power_dissipated_total_bulk == 'Y' ||
 		 		save_power_density_total_subvolumes == 'Y' )
				{
				
					FILE *power_dissipated_spectral; // append
					char dirPathPower_dissipated_spectral_subvolumes_FileName[260];
					sprintf(dirPathPower_dissipated_spectral_subvolumes_FileName, "%s/Q_omega_subvol.csv", results_folder); // path where the file is stored
					if (ig_0 == 0 && i_omega ==0){ power_dissipated_spectral = fopen(dirPathPower_dissipated_spectral_subvolumes_FileName, "w");} // write
					else if (ig_0 == 0 && i_omega!=0)
					{
						power_dissipated_spectral = fopen(dirPathPower_dissipated_spectral_subvolumes_FileName, "a"); // append
						fprintf(power_dissipated_spectral, "\n");
					}
					else { power_dissipated_spectral = fopen(dirPathPower_dissipated_spectral_subvolumes_FileName, "a"); }// append
					fprintf(power_dissipated_spectral, "%e,", Q_omega_subvol);					
					fclose(power_dissipated_spectral);
				
				}
				/*
				if ( save_power_dissipated_spectral_subvolumes == 'Y' ||
 		 		save_power_dissipated_total_subvolumes == 'Y' ||
 				save_power_dissipated_spectral_bulk == 'Y' ||
 		 		save_power_dissipated_total_bulk == 'Y' ||
 		 		save_power_density_total_subvolumes == 'Y' )
				{
					Q_subvol[ig_0][i_omega] = Q_omega_subvol;
					
					//if (ig_0 < const_N_subvolumes_per_object)
					//{
					//	Q_subvol[ig_0][i_omega] = Q_omega_subvol;
					//	Q_subvol[tot_sub_vol-ig_0-1][i_omega] = -Q_omega_subvol;
					//}
					
				}
				*/
			}
			//free(G_sys_TriangularMatrix);
			free(G_sys);
		}	//end if (solution =='I')
		free(alpha_0);	
	       	
		// save spectral transmissivity
		if (save_spectral_transmissivity == 'Y')
		{
			sprintf(spectral_transmissivity_folder, "%s/spectral_transmissivity", results_folder);
			create_folder(spectral_transmissivity_folder);
			FILE *spectral_transmissivity;
			char dirPathSpectral_trans_FileName[260];
			sprintf(dirPathSpectral_trans_FileName, "%s/%d.csv", spectral_transmissivity_folder, i_omega + 1); // path where the file is stored
			spectral_transmissivity = fopen(dirPathSpectral_trans_FileName, "w");
			fprintf(spectral_transmissivity, "%e", sum_trans_coeff);
			fclose(spectral_transmissivity);
		}

		// save spectral conductance
		for (int iTcalc = 0; iTcalc < N_Tcalc; iTcalc++) // EDIT VALUE :: change 1 to N_Tcalc for the temperature loop
		{
			double dtheta_dT = dtheta_dT_function(omega_value, Tcalc_vector[iTcalc], h_bar, k_b);
			//G_12_omega_SGF[i_omega][iTcalc] = dtheta_dT * sum_trans_coeff;
			double G_12_omega_SGF = dtheta_dT * sum_trans_coeff;
			if (save_spectral_conductance == 'Y' || save_total_conductance == 'Y')
			{
				FILE *spectral_conductance; // append
				char dirPathSpectral_cond_FileName[260];
				sprintf(dirPathSpectral_cond_FileName, "%s/G_omega_bulk_12_%eK.csv", results_folder, Tcalc_vector[iTcalc]); // path where the file is stored
				if (i_omega == 0)
					spectral_conductance = fopen(dirPathSpectral_cond_FileName, "w"); // write
				else
					spectral_conductance = fopen(dirPathSpectral_cond_FileName, "a"); // append
				//fprintf(spectral_conductance, "%e , %e\n", omega_value, G_12_omega_SGF[i_omega][iTcalc]);
				fprintf(spectral_conductance, "%e , %e\n", omega_value, G_12_omega_SGF);
				fclose(spectral_conductance);
			} // end if save_spectral_conductance

		} // END FOR T_calc LOOP

		// #################################################################
		// #############  Clear values for next frequency ##################
		// #################################################################

	} // END OMEGA VALUE LOOP FOR FREQUENCY RANGE ANALYSIS

	free(modulo_r_i_j);
	free(r_i_j_outer_r_i_j);

	// #################################################################
	// ################### Total-Post processing #######################
	// #################################################################
	
	if (const_N_omega > 1)
	{
		printf("end. \n");

		// implementation of trapezoidal numerical integration in C // https://stackoverflow.com/questions/25124569/implementation-of-trapezoidal-numerical-integration-in-c
		double step = 0.;
		double trapz = 0.;
		if ( save_power_dissipated_total_subvolumes == 'Y' ||
 		 save_power_dissipated_spectral_bulk == 'Y' ||
 		 save_power_dissipated_total_bulk == 'Y' ||
 		 save_power_density_total_subvolumes == 'Y' )
		{
			double(*Q_subvol)[const_N_omega] = calloc(tot_sub_vol, sizeof(*Q_subvol));
			if (Q_subvol == NULL){
				printf("Failure with memory after spectral analysis when spectral power dissipated is defined. ");
				return 1;
			}
			FILE *power_dissipated_spectral; // append
			char dirPathPower_dissipated_spectral_subvolumes_FileName[260];
			sprintf(dirPathPower_dissipated_spectral_subvolumes_FileName, "%s/Q_omega_subvol.csv", results_folder); // path where the file is stored
			power_dissipated_spectral = fopen(dirPathPower_dissipated_spectral_subvolumes_FileName, "r"); // read
			
			for (int i_omega = 0; i_omega < const_N_omega; i_omega++)
			{
				for (int i_subvol = 0; i_subvol < tot_sub_vol; i_subvol++)
				{
					fscanf(power_dissipated_spectral,"%lf ,", &Q_subvol[i_subvol][i_omega]);
				}
				fscanf(power_dissipated_spectral,"\n");	
			}
			fclose(power_dissipated_spectral);
					
			double(*trapz_Q) = malloc(sizeof *trapz_Q * tot_sub_vol); // Definition for trapezoidal integration. Used in total power dissipated
			if (trapz_Q == NULL){
				printf("Failure with memory after spectral analysis when spectral power dissipated is defined. ");
				return 1;
			}
			double sum_Q = 0.;
			double Q_bulk_1=0;
			double Q_bulk_2=0;
		
			for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) // tot_sub_vol
			{
				trapz_Q[ig_0] = 0.;
				for (int i = 1; i < const_N_omega; ++i)
				{
					step = omega[i] - omega[i - 1];
					trapz_Q[ig_0] += ((Q_subvol[ig_0][i] + Q_subvol[ig_0][i - 1]) / 2) * step;
				}
				if (ig_0<const_N_subvolumes_per_object)
				{
					Q_bulk_1 = Q_bulk_1+ trapz_Q[ig_0];
				}
				else
				{
					Q_bulk_2 = Q_bulk_2+ trapz_Q[ig_0];
				}
							
				if (save_power_dissipated_total_subvolumes == 'Y')
				{	
					FILE *power_dissipated_total; // append
					char dirPathPower_dissipated_FileName[260];
					sprintf(dirPathPower_dissipated_FileName, "%s/Q_total_subvol.csv", results_folder); // path where the file is stored
					if (ig_0 == 0)
						power_dissipated_total = fopen(dirPathPower_dissipated_FileName, "w"); // write
					else
						power_dissipated_total = fopen(dirPathPower_dissipated_FileName, "a"); // append
					fprintf(power_dissipated_total, "%e\n", trapz_Q[ig_0]);
					fclose(power_dissipated_total);
				}
						
			} // end for ig_0 = 0 to  tot_sub_vol

			if (save_power_density_total_subvolumes == 'Y')
			{		
				double(*Q_density) = malloc(sizeof *Q_density * tot_sub_vol);
				if (Q_density == NULL){
					printf("Failure with memory when power dissipated density is defined. ");
					return 1;
				}
				FILE *power_density; // append
				char dirPathPower_density_FileName[260];
				sprintf(dirPathPower_density_FileName, "%s/Q_density_subvol.csv", results_folder); // path where the file is stored
				for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) // tot_sub_vol
				{
					if (ig_0 == 0)
						power_density = fopen(dirPathPower_density_FileName, "w"); // write
					else
						power_density = fopen(dirPathPower_density_FileName, "a"); // append
					Q_density[ig_0] = trapz_Q[ig_0]/delta_V_vector[ig_0];
					fprintf(power_density, "%e\n", Q_density[ig_0]);
					fclose(power_density);
				}
				free(Q_density);
			} // end if save_power_density_total_subvolumes
		
			free(trapz_Q);
		
			if (save_power_dissipated_total_bulk == 'Y')
			{
				FILE *power_dissipated_total_bulk; // append
				char dirPathPower_dissipated_bulk_FileName[260];
				sprintf(dirPathPower_dissipated_bulk_FileName, "%s/Q_total_bulk.csv", results_folder); // path where the file is stored
				power_dissipated_total_bulk = fopen(dirPathPower_dissipated_bulk_FileName, "w"); // write
				fprintf(power_dissipated_total_bulk, "%e,%e\n", Q_bulk_1,Q_bulk_2);
				fclose(power_dissipated_total_bulk);
			}
		
			if (save_power_dissipated_spectral_bulk == 'Y')
			{
				double(*Q_omega_bulk_1) = malloc(sizeof *Q_omega_bulk_1 * const_N_omega);
				if (Q_omega_bulk_1 == NULL){
					printf("Failure with memory when bulk power dissipated is defined. ");
					return 1;
				}
				double(*Q_omega_bulk_2) = malloc(sizeof *Q_omega_bulk_2 * const_N_omega);
				if (Q_omega_bulk_2 == NULL){
					printf("Failure with memory when bulk power dissipated is defined. ");
					return 1;
				}
				FILE *power_dissipated_spectral_bulk; // append
				char dirPathPower_dissipated_spectral_bulk_FileName[260];
				sprintf(dirPathPower_dissipated_spectral_bulk_FileName, "%s/Q_omega_bulk.csv", results_folder); // path where the file is stored
			
				for (int i_omega = 0; i_omega < const_N_omega; ++i_omega)
				{
					Q_omega_bulk_1[i_omega] = 0;
					Q_omega_bulk_2[i_omega] = 0;
					if (i_omega == 0) power_dissipated_spectral_bulk = fopen(dirPathPower_dissipated_spectral_bulk_FileName, "w"); // write
					else power_dissipated_spectral_bulk = fopen(dirPathPower_dissipated_spectral_bulk_FileName, "a"); // append
					for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
					{
						if (ig_0<const_N_subvolumes_per_object)	Q_omega_bulk_1[i_omega] = Q_omega_bulk_1[i_omega] + Q_subvol[ig_0][i_omega];
						else Q_omega_bulk_2[i_omega] = Q_omega_bulk_2[i_omega] + Q_subvol[ig_0][i_omega];
					}
					fprintf(power_dissipated_spectral_bulk, "%e, %e\n", Q_omega_bulk_1[i_omega],Q_omega_bulk_2[i_omega]);
					fclose(power_dissipated_spectral_bulk);
				}
				free(Q_omega_bulk_1);
				free(Q_omega_bulk_2);
			}
		
			free(Q_subvol);
		} // if power

		double sum_conductance = 0.;
		if (save_total_conductance == 'Y') 
		{
		double(*Total_conductance) = malloc(sizeof *Total_conductance * N_Tcalc);
		if (Total_conductance == NULL){
			printf("Failure with memory when total conductance is defined. ");
			return 1;
		}
		for (int iTcalc = 0; iTcalc < N_Tcalc; iTcalc++) // EDIT VALUE :: change 1 to N_Tcalc for the temperature loop
		{
			sum_conductance = 0.;
			step = 0.;
			trapz = 0.;
			G_12_total_SGF_from_omega = 0.;
			
			double(*G_12_omega_SGF)[N_Tcalc] = calloc(const_N_omega, sizeof(*G_12_omega_SGF)); // spectral conductance
			if (G_12_omega_SGF == NULL){
				printf("Failure with memory before spectral analysis when conductance defined. ");
				return 1;
			}
			FILE *spectral_conductance; // append
			char dirPathSpectral_cond_FileName[260];
			char buffer[260];
			sprintf(dirPathSpectral_cond_FileName, "%s/G_omega_bulk_12_%eK.csv", results_folder, Tcalc_vector[iTcalc]); // path where the file is stored
			spectral_conductance = fopen(dirPathSpectral_cond_FileName, "r"); // read	
			for (int i_omega = 0; i_omega < const_N_omega; i_omega++)
			{
				//fprintf(spectral_conductance, "%e , %e\n", omega_value, G_12_omega_SGF[i_omega][iTcalc]);
				fscanf(spectral_conductance, "%s , %lf\n", buffer, &G_12_omega_SGF[i_omega][iTcalc]);
			}
			fclose(spectral_conductance);

			for (int i = 1; i < const_N_omega; ++i)
			{
				step = omega[i] - omega[i - 1];

				sum_conductance = (G_12_omega_SGF[i][iTcalc] + G_12_omega_SGF[i - 1][iTcalc]) / 2;
				trapz += sum_conductance * step;
			}
			free(G_12_omega_SGF);
			G_12_total_SGF_from_omega = (fabs(trapz) / (2. * pi)); // Total conductance [W/K]
			Total_conductance[iTcalc] = G_12_total_SGF_from_omega;
			
			FILE *Total_conductance_file; // append
			char dirPath_Total_conductance_FileName[260];
			sprintf(dirPath_Total_conductance_FileName, "%s/G_bulk_12.csv", results_folder); // path where the file is stored
			if (iTcalc == 0)
				Total_conductance_file = fopen(dirPath_Total_conductance_FileName, "w"); // write
			else
				Total_conductance_file = fopen(dirPath_Total_conductance_FileName, "a"); // append
			fprintf(Total_conductance_file, "%e, %e\n", Tcalc_vector[iTcalc], Total_conductance[iTcalc]);
			fclose(Total_conductance_file);
			
		} // END FOR T_calc LOOP

		printf("Total conductance at %e K= %e \n", Tcalc_vector[2], Total_conductance[2]);
		free(Total_conductance);
		
		}

	} // end if const_N_omega>1
		
		
	time_t simulation_time;
	time(&simulation_time);
				

	total_memory = get_mem_usage(); // measure post-processing memory usage
	{
		FILE * memory; //append
		char dirPathMemory_FileName[260];
		sprintf(dirPathMemory_FileName, "matlab_scripts/memory/memory_analysis.csv"); // path where the file is stored
		memory = fopen(dirPathMemory_FileName, "a"); //append
		if  (multithread == 'Y') {fprintf(memory,"%d,%c,Parallel,%ld,%ld,%ld,%ld,%ld s\n",tot_sub_vol,solution,baseline, pre_solver_memory, calculation_memory, total_memory,simulation_time-simulation_start); }
		else if (multithread == 'N'){ fprintf(memory,"%d,%c,Serial,%ld,%ld,%ld,%ld,%ld s\n",tot_sub_vol,solution,baseline, pre_solver_memory, calculation_memory, total_memory,simulation_time-simulation_start);}
		fclose(memory);
	}
	 
	printf("\nThe results can be accessed in the folder:\n %s\n",results_folder);
	//free(Total_conductance);
	free(omega);
	free(delta_V_vector);

	//free(G_12_omega_SGF);

	free(results_folder);

} // end of code

// #################################################################
// ##################### END OF THE CODE ###########################
// #################################################################
