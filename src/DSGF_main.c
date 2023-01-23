// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Discrete System Green's Function 
// Near-field radiative heat transfer framework between thermal objects 
// Developed by RETL group at the University of Utah, USA

// LAST UPDATE: January 09, 2023
// 
// In this version:
//	- The following .txt files remove recompile the code need for user modifications:
//	    - const_N_subvolumes_per_object.txt defines the number of subvolumes per thermal object
//	    - const_N_bulk_objects.txt defines the number of bulk objects
//	    - const_N_omega.txt defines the number of frequencies to evaluate
//	    - T_calc.txt defines the temperature where the conductance is calculated
//	    - control.txt contains other user definitions such as type of geometry, material, solution as saving options
//	    - after selecting the geometry, go to the chosen geometry .txt file to set the dimensions and temperatures of the thermal objects
//  - More info in read_me-user_inputs.txt 
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// General c libraries
#include<stdio.h>
#include<math.h> 

#include<complex.h>  // complex numbers library, code must be compiled as c99 standard https://stackoverflow.com/questions/6418807/how-to-work-with-complex-numbers-in-c
#include <stdlib.h> // export/import data

// Libraries for creating directories and files using a loop in C. Sources: https://stackoverflow.com/questions/46612504/creating-directories-and-files-using-a-loop-in-c and https://stackoverflow.com/questions/7430248/creating-a-new-directory-in-c 
#include <string.h> // library used to concatenate 2 strings https://stackoverflow.com/questions/46612504/creating-directories-and-files-using-a-loop-in-c

// Library with the inputs and functions for DSGF
#include "user_inputs.h" // User inputs definitions header. No values are defined in this file.  
#include "functions_DSGF.h" // Definitions of functions used in DSGF
#include "file_utils.h" // header with definitions of read_user_inputs and read_calculation_temperatures functions
#include "computational/solvers/iterative_solver.h"
#include "array_functions.h"

#include "time.h"

#include "geometry/sphere.h"
#include "geometry/thin_film.h"
#include "geometry/shared.h"

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

//#include <omp.h> // library for OpenMP

// #################################################################
// ################### START OF THE CODE ###########################
// #################################################################


int main()
{

	// ----------- CONSTANTS ------------------

	const double pi = 3.14159265359;         // pi number   
	const double q = 1.602176634e-19;        // number of Joules per eV
	const double h_bar = 1.054571817e-34;    // Planck's constant [J*s]
	const double k_b = 1.38064852e-23;       // Boltzmann constant [J/K]
	const double epsilon_0 = 8.8542e-12;     // Permittivity of free space [F/m]
	const double c_0 = 299792458;            // Speed of light in vacuum [m/s]
	const double mu_0 = (4.*pi)*pow(10,-7);        // Permeability of free space [H/m]
	const double epsilon_ref = 1.;             	 // dielectric function of the background reference medium


	long baseline = get_mem_usage(); // measure baseline memory usage
	clock_t begin = clock();  /* set timer here, do your time-consuming job */

	int N_subvolumes_per_object, N_bulk_objects, N_omega;

	read_user_control(geometry, material, &solution, &single_spectrum_analysis, &save_spectral_conductance, &save_spectral_transmissivity, &save_power_dissipated, &N_bulk_objects, &N_omega, &N_subvolumes_per_object);

	read_calculation_temperatures(N_Tcalc, Tcalc_vector);

	int const const_N_subvolumes_per_object = N_subvolumes_per_object;

	int const const_N_bulk_objects = N_bulk_objects;

	int const const_N_omega = N_omega;

	int tot_sub_vol = const_N_subvolumes_per_object*const_N_bulk_objects; // Assign tot_sub_vol: Computes the total number of subvolumes in the system. tot_sub_vol is defined this way because it must be a non-variable parameter due to the computations perfomed in the code. Previously, it was defined as #define tot_sub_vol const_N_subvolumes_per_object*const_N_bulk_objects


	// ####################################    
	// #### Dynamic memory allocation: ####    
	double (*R)[3] = calloc(tot_sub_vol, sizeof(*R)); // center of subvolumes for thermal objects: info imported from a .txt file

	// radial frequency [rad/s]
	double (*omega) = malloc(sizeof *omega *const_N_omega); 

	double (*delta_V_vector) = malloc(sizeof *delta_V_vector *tot_sub_vol);  //Vector of all subvolume size. Combines delta_V_1 and delta_V_2 in the same array
	double (*T_vector) = malloc(sizeof *T_vector *tot_sub_vol); // (N x 1) vector of all subvolume temperatures [K]
	double complex (*alpha_0) = malloc(sizeof *alpha_0 *tot_sub_vol); //Bare polarizability [m^3]

	double (*G_12_omega_SGF)[N_Tcalc] = calloc(const_N_omega, sizeof(*G_12_omega_SGF));

	double (*modulo_r_i_j)[tot_sub_vol] = malloc(sizeof *modulo_r_i_j * tot_sub_vol);
	double (*sum_trans_coeff) = calloc(const_N_omega, sizeof(*sum_trans_coeff)); 

	double (*Q_subvol)[const_N_omega] = calloc(tot_sub_vol, sizeof(*Q_subvol));

	// ######### Properties for thermal objects ###########
	printf("Simulation for a total of %d dipoles in %d thermal objects\n",tot_sub_vol,const_N_bulk_objects);

	double delta_V_1, delta_V_2;
	if(strcmp(geometry,"sphere")==0)
	{
		set_up_sphere_geometry(pi, tot_sub_vol, const_N_subvolumes_per_object, &T1, &T2, &d, &delta_V_1, &delta_V_2, R);
	}   

	if(strcmp(geometry,"thin-films")==0)
	{
		set_up_thin_film_geometry(tot_sub_vol, const_N_subvolumes_per_object, const_N_bulk_objects, &T1, &T2, &d, &delta_V_1, &delta_V_2, R);
	}
	set_delta_V_vector_T_vector(T1, T2, delta_V_1, delta_V_2, tot_sub_vol, const_N_subvolumes_per_object, T_vector, delta_V_vector);


	printf("d = %e m \n",d);

	char *results_folder = set_up_results(material, geometry, tot_sub_vol, d); // Folders for results 

	if(save_power_dissipated =='Y'){
		//EXPORT R
		char dirPathVector_subvolumes_lattice_FileName[260];
		sprintf(dirPathVector_subvolumes_lattice_FileName, "%s/vector_subvolumes_lattice.csv",results_folder); // path where the file is stored
		write_to_csv_double_matrix(dirPathVector_subvolumes_lattice_FileName, tot_sub_vol, 3, R);
		
		//EXPORT delta_V_vector
		char dirPathVector_subvolumes_volume_FileName[260];
		sprintf(dirPathVector_subvolumes_volume_FileName, "%s/vector_subvolumes_volume.csv",results_folder); // path where the file is stored
		write_to_csv_double_array(dirPathVector_subvolumes_volume_FileName, tot_sub_vol, delta_V_vector);
	} // end if save_power_dissipated


	double initial,final;

	if(strcmp(material,"SiO2")==0 || strcmp(material,"SiC")==0 || strcmp(material,"u-SiC")==0) 
	{
		//Uniform spectrum
		initial = 5.e-6;
		final = 25.e-6;
		double lambda[const_N_omega];
		double_linspace(initial, final, const_N_omega, lambda);
		for(int i_lambda = 0; i_lambda < const_N_omega; i_lambda++)
		{
			omega[i_lambda] = 2.*pi*c_0/lambda[i_lambda];  // Radial frequency [rad/s]
		}
	}	

	else if(strcmp(material,"SiN")==0) //cannot compare strings in C with ==; source: https://ideone.com/BrFA00
	{
		//Uniform spectrum: 
		initial = 2.e13;
		final = 3.e14;
		double_linspace(initial, final, const_N_omega, omega);
	}


	//################## FREE-SPACE GREEN'S FUNCTION AND INTERACTION A MATRIX ##################### 
	// Fill terms for G^0: 

	//Definitions for G^0: 
	double complex (*r_i_j_outer_r_i_j)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*r_i_j_outer_r_i_j));  
	setup_G_0_matrices(tot_sub_vol, modulo_r_i_j, r_i_j_outer_r_i_j, R);

	// #################################################################
	// ################## FREQUENCY RANGE ANALYSIS #####################
	// #################################################################
	//Loop to analyze a range of desired frequencies
	printf("----- Spectrum range calculation -----\n");

	double omega_range;
	if(single_spectrum_analysis =='Y') omega_range=1;
	if(single_spectrum_analysis =='N') omega_range=const_N_omega;

	//	#pragma omp parallel for
	for (int i_omega = 0; i_omega < omega_range; i_omega++) // Frequency loop
	{

		if (i_omega==1)
		{
			single_analysis = 'n';
		}


		double omega_value = omega[i_omega]; // omega definition is on line 212
		printf("%d) omega = %e. ", i_omega+1,omega_value);

		double complex epsilon;
		if(strcmp(material,"SiO2")==0) //cannot compare strings in C with ==; source: https://ideone.com/BrFA00
		{
			epsilon = calculate_epsilon_SiO2(q, omega_value, h_bar);
		}

		else if(strcmp(material,"SiC")==0) 
		{  
			epsilon = calculate_epsilon_SiC(omega_value);
		}
		else if(strcmp(material,"u-SiC")==0) 
		{  
			epsilon = calculate_epsilon_u_SiC(omega_value, pi, c_0);
		}

		else if(strcmp(material,"SiN")==0) 
		{  
			epsilon = calculate_epsilon_SiN(omega_value, pi);
		}


		double k_0=k_0_function(omega_value, epsilon_0, mu_0) ; //wave vector in free space

		double k=k_function(omega_value, epsilon_ref, epsilon_0, mu_0); //wave vector in reference medium

		for (int i_alpha = 0; i_alpha < tot_sub_vol; i_alpha++)
		{
			alpha_0[i_alpha] = delta_V_vector[i_alpha]*(epsilon - epsilon_ref); //Bare polarizability [m^3]
		}

		// ################### MATRICES STRUCTURE LOOPS ###########################
		// 3N X 3N Matrices structure loops for G^0 and A:

		// Linear system AG=G^0 
		double complex (*G_0)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*G_0)); 
		double complex (*A)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*A));

		get_G0_A_matrices(tot_sub_vol, G_0, A, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, alpha_0, delta_V_vector);

		//printf("##################### \n SOLVE LINEAR SYSTEM AG=G^0 \n##################### \n");
		//printf("##################### \n LAPACK/LAPACKE ZGELS ROUTINE \n##################### \n");
		//Description of ZGELS: https://extras.csc.fi/math/nag/mark21/pdf/F08/f08anf.pdf

		double complex (*G_sys)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys));

		if(solution =='D')
		{
			printf("Direct inversion status: ");
			direct_solver(tot_sub_vol, A, G_0, G_sys);
			printf("concluded\n");
		}
		free(A);


		if(solution =='I')
		{ 
			printf("Iterative status:\n m= ");
			iterative_solver(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_0, G_sys);
			printf("Final solution:\n");
		}

		free(G_0);


		// #################################################################
		//  Spectral and total conductance between bulk objects at temp. T
		// ###################  Thermal power dissipated ###################
		// #################################################################

		spectral_total_conductance_thermal_power(tot_sub_vol, i_omega, const_N_omega, k_0, h_bar, k_b, epsilon, omega_value, T_vector, delta_V_vector, const_N_subvolumes_per_object, pi, G_sys, sum_trans_coeff, Q_subvol);
		free(G_sys);


		if(save_spectral_transmissivity =='Y'){

			sprintf(spectral_transmissivity_folder, "%s/spectral_transmissivity",results_folder); 
			create_folder(spectral_transmissivity_folder);

			{
				FILE * spectral_transmissivity; 
				char dirPathSpectral_trans_FileName[260];
				sprintf(dirPathSpectral_trans_FileName, "%s/%d.csv",spectral_transmissivity_folder,i_omega+1); // path where the file is stored
				spectral_transmissivity =fopen(dirPathSpectral_trans_FileName,"w"); 
				fprintf(spectral_transmissivity,"%e",sum_trans_coeff[i_omega]);      
				fclose(spectral_transmissivity);
			}
		}

		for (int iTcalc = 0; iTcalc < N_Tcalc; iTcalc++) // EDIT VALUE :: change 1 to N_Tcalc for the temperature loop
		{
			double dtheta_dT = dtheta_dT_function(omega_value,Tcalc_vector[iTcalc], h_bar, k_b); 
			G_12_omega_SGF[i_omega][iTcalc] = dtheta_dT*sum_trans_coeff[i_omega];
			if(save_spectral_conductance =='Y'){
				{
					FILE * spectral_conductance; //append
					char dirPathSpectral_cond_FileName[260];
					sprintf(dirPathSpectral_cond_FileName, "%s/spectral_conductance_%eK.csv", results_folder, Tcalc_vector[iTcalc]); // path where the file is stored
					if(i_omega == 0) spectral_conductance =fopen(dirPathSpectral_cond_FileName,"w"); //write
					else spectral_conductance = fopen(dirPathSpectral_cond_FileName, "a"); //append
					fprintf(spectral_conductance,"%e ; %e\n",omega_value,G_12_omega_SGF[i_omega][iTcalc]); 
					fclose(spectral_conductance);
				}
			} // end if save_spectral_conductance

		} // END FOR T_calc LOOP


		// #################################################################
		// #############  Clear values for next frequency ##################
		// #################################################################

		memset(sum_trans_coeff, 0, sizeof sum_trans_coeff);

	} // END OMEGA VALUE LOOP FOR FREQUENCY RANGE ANALYSIS

	free(R);
	free(sum_trans_coeff);
	free(modulo_r_i_j);
	free(r_i_j_outer_r_i_j);

	free(alpha_0);

	double (*Total_conductance) = malloc(sizeof *Total_conductance *N_Tcalc); 
	double (*trapz_Q) = malloc(sizeof *trapz_Q *tot_sub_vol); // Definition for trapezoidal integration. Used in total power dissipated

	double trapz=0.;
	if( const_N_omega > 1)
	{
		printf("\nEnd of frequency loop\n");
		//printf("----- Total conductance -----\n");


		// implementation of trapezoidal numerical integration in C // https://stackoverflow.com/questions/25124569/implementation-of-trapezoidal-numerical-integration-in-c
		double sum_conductance = 0.;
		double step = 0.;
		double sum_Q = 0.;

		for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
		{
			trapz_Q[ig_0] = 0.;
			for (int i=1; i < const_N_omega; ++i)
			{
				step = omega[i] - omega[i-1];
				trapz_Q[ig_0] +=  ((Q_subvol[ig_0][i]+Q_subvol[ig_0][i-1])/2)*step;
			}

			if(save_power_dissipated =='Y'){
				{
					FILE * power_dissipated; //append
					char dirPathPower_dissipated_FileName[260];
					sprintf(dirPathPower_dissipated_FileName, "%s/power_dissipated.csv",results_folder); // path where the file is stored
					if(ig_0 == 0) power_dissipated =fopen(dirPathPower_dissipated_FileName,"w"); //write
					else power_dissipated = fopen(dirPathPower_dissipated_FileName, "a"); //append
					fprintf(power_dissipated,"%e\n",trapz_Q[ig_0]); 
					fclose(power_dissipated);
				}
			} // end if save_power_dissipated

		}



		for (int iTcalc = 0; iTcalc < N_Tcalc; iTcalc++) // EDIT VALUE :: change 1 to N_Tcalc for the temperature loop
		{      
			sum_conductance=0.;
			step=0.;
			trapz=0.;
			G_12_total_SGF_from_omega=0.;

			for (int i=1; i < const_N_omega; ++i)
			{

				step = omega[i]-omega[i-1];

				sum_conductance = (G_12_omega_SGF[i][iTcalc]+G_12_omega_SGF[i-1][iTcalc])/2;
				trapz += sum_conductance*step;

			}
			G_12_total_SGF_from_omega = (fabs(trapz)/(2.*pi)); // Total conductance [W/K]
			Total_conductance[iTcalc]=G_12_total_SGF_from_omega;

		} // END FOR T_calc LOOP

		printf("Total conductance at %e K= %e \n",Tcalc_vector[2],Total_conductance[2]);    

	} //end if const_N_omega>1

	printf("\n");

	sprintf(dirPathpos_processing_summary_FileName, "%s/pos_processing_summary.txt",results_folder); // path where the file is stored

	clock_t end = clock(); //end timer
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC; //calculate time for the code

	create_pos_processing(dirPathpos_processing_summary_FileName, material, initial, final, time_spent, Tcalc_vector, Total_conductance, N_Tcalc);


	free(Total_conductance); 

	printf("Usage: %ld + %ld = %ld kb\n", baseline, get_mem_usage()-baseline,get_mem_usage());

	free(omega);
	free(delta_V_vector);
	free(T_vector);

	free(G_12_omega_SGF);

	free(trapz_Q);
	free(Q_subvol);
	free(results_folder);

} //end of code 

// #################################################################
// ##################### END OF THE CODE ###########################
// #################################################################
