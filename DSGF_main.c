// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Discrete System Green's Function 
// Near-field radiative heat transfer framework between thermal objects 
// Developed by RETL group at the University of Utah, USA

// LAST UPDATE: SEPTEMBER 19, 2022
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
#include "iterative_solver.h"
#include "indexing_util.h"

// debugging utils
#include "debugging/debugging_utils.h"

// LAPACKE libraries: https://www.netlib.org/lapack/lapacke.html ; https://extras.csc.fi/math/nag/mark21/pdf/F08/f08anf.pdf
#include <lapacke.h> 
#include "lapack_header.h" //header with Lapack definitions

//#include <omp.h> // library for OpenMP

FILE * pos_processing_summary; // call the main outputs' file,  

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
	
	//read_user_inputs(material, geometry, discretization_thin_film, &d, &radius, &Lx, &Ly, &Lz, &T1, &T2, &solution, &single_spectrum_analysis, &save_A_matrix, &save_G0_matrix, &save_SGF_matrix, &save_spectral_conductance, &save_spectral_transmissivity, &save_power_dissipated);
	
	read_user_control(geometry, material, &solution, &single_spectrum_analysis, &save_spectral_conductance, &save_spectral_transmissivity, &save_power_dissipated);

	read_calculation_temperatures(N_Tcalc, Tcalc_vector);

	int const const_N_subvolumes_per_object = read_int_from_file(N_subvolumes_per_object_file);

	int const const_N_bulk_objects = read_int_from_file(N_bulk_objects_file);

	int const const_N_omega = read_int_from_file(N_omega_file);

	size_t tot_sub_vol = const_N_subvolumes_per_object*const_N_bulk_objects; // Assign tot_sub_vol: Computes the total number of subvolumes in the system. tot_sub_vol is defined this way because it must be a non-variable parameter due to the computations perfomed in the code. Previously, it was defined as #define tot_sub_vol const_N_subvolumes_per_object*const_N_bulk_objects

	subvol shape_file[const_N_subvolumes_per_object]; //typedef struct node 
	subvol shape_filetf[tot_sub_vol];

	// ####################################    
	// #### Dynamic memory allocation: ####    
	// Links: https://stackoverflow.com/questions/13534966/how-to-dynamically-allocate-a-contiguous-block-of-memory-for-a-2d-array https://stackoverflow.com/questions/39108092/allocating-contiguous-memory-for-a-3d-array-in-c

	//int (*a)[sz[1]][sz[2]] = calloc(sz[0], sizeof(*a));

	double (*R)[3] = calloc(tot_sub_vol, sizeof(*R));//double R[tot_sub_vol][3]; // center of subvolumes for thermal objects: info imported from a .txt file

	// radial frequency [rad/s]
	double (*x_omega) = malloc(sizeof *x_omega *const_N_omega); 
	double (*lambda) = malloc(sizeof *lambda *const_N_omega); 
	double (*omega) = malloc(sizeof *omega *const_N_omega); 
							  //double E_joules[const_N_omega];
							  //double E_ev[const_N_omega];

	double (*delta_V_vector) = malloc(sizeof *delta_V_vector *tot_sub_vol);  //Vector of all subvolume size. Combines delta_V_1 and delta_V_2 in the same array
	double (*T_vector) = malloc(sizeof *T_vector *tot_sub_vol); // (N x 1) vector of all subvolume temperatures [K]
	double complex (*alpha_0) = malloc(sizeof *alpha_0 *tot_sub_vol); //Bare polarizability [m^3]

	//double (*G_12_omega_SGF) = malloc(sizeof *G_12_omega_SGF *const_N_omega); //double G_12_omega_SGF[const_N_omega]; //  //double G_12_omega_SGF[tot_sub_vol][tot_sub_vol];
	double (*G_12_omega_SGF)[N_Tcalc] = calloc(const_N_omega, sizeof(*G_12_omega_SGF));

	double (*modulo_r_i_j)[tot_sub_vol] = malloc(sizeof *modulo_r_i_j * tot_sub_vol);
	double complex (*G_element)[tot_sub_vol] = malloc(sizeof *G_element * tot_sub_vol); 
	double (*sum_trans_coeff) = malloc(sizeof *sum_trans_coeff *const_N_omega); 

	double complex (*Q_omega_subvol) = malloc(sizeof * Q_omega_subvol * tot_sub_vol); // power dissipated per subvolume
	double (*Q_omega_thermal_object)[const_N_omega] = malloc(sizeof * Q_omega_thermal_object * const_N_bulk_objects);  
	double (*Q_subvol)[const_N_omega] = calloc(tot_sub_vol, sizeof(*Q_subvol));
	/* 
	// The SGF calculation considering the wavelength is not considered, so it is commented. 
	double (*G_12_lambda_SGF)[tot_sub_vol] = calloc(tot_sub_vol, sizeof(*G_12_lambda_SGF)); //double G_12_lambda_SGF[tot_sub_vol][tot_sub_vol];
	double (*trans_coeff_lambda)[tot_sub_vol] = calloc(tot_sub_vol, sizeof(*trans_coeff_lambda)); // double trans_coeff_lambda[tot_sub_vol][tot_sub_vol];
	free(G_12_lambda_SGF);
	free(trans_coeff_lambda);
	*/

	// ######### Properties for thermal objects ###########
	printf("Simulation for a total of %d dipoles in %d thermal objects\n",tot_sub_vol,const_N_bulk_objects);

	
	// Import data from txt file into C program: https://stackoverflow.com/questions/22745457/import-data-from-txt-file-into-c-program ; https://stackoverflow.com/questions/49563003/importing-data-from-txt-file-into-c-program ; https://www.cs.utah.edu/~germain/PPS/Topics/C_Language/file_IO.html
	FILE *import_discretization;  // https://stackoverflow.com/questions/22745457/import-data-from-txt-file-into-c-program
	int i_import = 0;

	char dirPathFileNameDISCRETIZATION[260]; // https://stackoverflow.com/questions/46612504/creating-directories-and-files-using-a-loop-in-c 

	if(strcmp(geometry,"sphere")==0) //cannot compare strings in C with ==; source: https://ideone.com/BrFA00
	{
		read_geometry_sphere(&d, &radius, &T1, &T2);
		radius1 = radius; // perfect same-sized spheres
		radius2 = radius; // perfect same-sized spheres
		vol1 = vol_sphere(radius1, pi); // calls function that calculates the volume for the sphere 1
		vol2 = vol_sphere(radius2, pi); // calls function that calculates the volume for the sphere 2
		delta_V_1 = vol1/const_N_subvolumes_per_object; // defines the subvolumes' volume for sphere 1
		delta_V_2 = vol2/const_N_subvolumes_per_object; // defines the subvolumes' volume for sphere 2
	
		sprintf(dirPathFileNameDISCRETIZATION, "discretizations/sphere_subvol_%d.txt",const_N_subvolumes_per_object); // path where the file is stored
		import_discretization = fopen(dirPathFileNameDISCRETIZATION, "r"); // "r" = read, "w" = write
		while (3 == fscanf(import_discretization, "%e %e %e", &shape_file[i_import].x, &shape_file[i_import].y, &shape_file[i_import].z))
		{   
			i_import++;
		}
		double origin1[3] = {radius1,radius1,radius1}; // first index is 0
		double origin2[3]= {origin1[0]+radius1+d+radius2,origin1[1]+radius2-radius1,origin1[2]+radius2-radius1}; // first index is 0
		for (int i_subvol=0; i_subvol<tot_sub_vol;i_subvol++) //tot_sub_vol
		{
			if(i_subvol<const_N_subvolumes_per_object)
			{
				R[i_subvol][0] = shape_file[i_subvol].x *pow(delta_V_1,1./3) + origin1[0]; // Lindsay's code: R1*pow(delta_V_1,1/3) + origin1[i]
				R[i_subvol][1] = shape_file[i_subvol].y *pow(delta_V_1,1./3) + origin1[1]; // Lindsay's code: R1*pow(delta_V_1,1/3) + origin1[i]
				R[i_subvol][2] = shape_file[i_subvol].z *pow(delta_V_1,1./3) + origin1[2]; // Lindsay's code: R1*pow(delta_V_1,1/3) + origin1[i]
													   //printf("delta_V_1 = %e\n",pow(delta_V_1,1./3));
													   //printf("R[%d] = %e, %e, %e \n",i_subvol,R[i_subvol][0],R[i_subvol][1],R[i_subvol][2]);
			}
			else
			{
				R[i_subvol][0] = shape_file[i_subvol-const_N_subvolumes_per_object].x* pow(delta_V_2,1./3) + origin2[0]; // Lindsay's code: R2*pow(delta_V_2,1/3) + origin2[i]
				R[i_subvol][1] = shape_file[i_subvol-const_N_subvolumes_per_object].y*pow(delta_V_2,1./3) + origin2[1]; // Lindsay's code: R2*pow(delta_V_2,1/3) + origin2[i]
				R[i_subvol][2] = shape_file[i_subvol-const_N_subvolumes_per_object].z*pow(delta_V_2,1./3) + origin2[2]; // Lindsay's code: R2*pow(delta_V_2,1/3) + origin2[i]    
																  //printf("R = %e, %e, %e \n",R[i_subvol][0],R[i_subvol][1],R[i_subvol][2]);
			}

		}   
	}   


	if(strcmp(geometry,"thin-films")==0) //cannot compare strings in C with ==; source: https://ideone.com/BrFA00
	{
		read_geometry_thin_films(&d, &Lx, &Ly, &Lz, &T1, &T2);
		
		vol1 = Lx*Ly*Lz; // calculates the volume for membrane 1
		vol2 = vol1;     // defines the volume of membrane 2 = membrane 1
		delta_V_1 = vol1/const_N_subvolumes_per_object; // defines the subvolumes' volume for membrane 1
		delta_V_2 = vol2/const_N_subvolumes_per_object;  // defines the subvolumes' volume for membrane 1
		
		int Lx_int, Ly_int, Lz_int, d_int;
		Lx_int = Lx*pow(10,9); 
		Ly_int = Ly*pow(10,9); 
		Lz_int = Lz*pow(10,9); 
		d_int = d*pow(10,9); 
		//sprintf(dirPathFileNameDISCRETIZATION, "discretizations/%s.txt", discretization_thin_film);
		sprintf(dirPathFileNameDISCRETIZATION, "discretizations/%d_thin_films_Lx%dnm_Ly%dnm_Lz%dnm_d%dnm_N%d_discretization.txt",const_N_bulk_objects,  Lx_int, Ly_int, Lz_int, d_int, tot_sub_vol); //under modification
		import_discretization = fopen(dirPathFileNameDISCRETIZATION, "r"); // "r" = read, "w" = write
		while (3 == fscanf(import_discretization, "%e %e %e", &shape_filetf[i_import].x, &shape_filetf[i_import].y, &shape_filetf[i_import].z))
		{   
			i_import++;
		}
		for (int i_subvol=0; i_subvol<tot_sub_vol;i_subvol++) //tot_sub_vol
		{
			R[i_subvol][0] = shape_filetf[i_subvol].x ; // Lindsay's code: R1*pow(delta_V_1,1/3) + origin1[i]
			R[i_subvol][1] = shape_filetf[i_subvol].y ; // Lindsay's code: R1*pow(delta_V_1,1/3) + origin1[i]
			R[i_subvol][2] = shape_filetf[i_subvol].z ; // Lindsay's code: R1*pow(delta_V_1,1/3) + origin1[i]
								    //printf("R[%d] = %e, %e, %e \n",i_subvol,R[i_subvol][0],R[i_subvol][1],R[i_subvol][2]);
		} 
	}



	fclose(import_discretization);

	printf("d = %e m \n",d);
	
	char *results_folder = set_up_results(material, geometry, tot_sub_vol, d); // Folders for results 

	
	for (int i_vec=0; i_vec<tot_sub_vol; i_vec++)
	{
		if (i_vec < const_N_subvolumes_per_object) // 2-body case
		{
			delta_V_vector[i_vec] = delta_V_1;
			T_vector[i_vec] = T1;
		}
		else
		{
			delta_V_vector[i_vec] = delta_V_2;
			T_vector[i_vec] = T2;
		}
	}  

	double initial,final;
	
	if(strcmp(material,"SiO2")==0) 
	{
	//Uniform spectrum
	double initial_lambda = 5.e-6;//5.e-6;
	double final_lambda = 25.e-6; //25.e-6;
	initial = initial_lambda;
	final = final_lambda;
	double step_lambda = (final-initial)/(const_N_omega-1); // linspace in C: https://stackoverflow.com/questions/60695284/linearly-spaced-array-in-c
	for(int i_lambda = 0; i_lambda < const_N_omega; i_lambda++)
	{
		lambda[i_lambda]= initial + i_lambda*step_lambda; // Wavelength [m]
		omega[i_lambda] = 2.*pi*c_0/lambda[i_lambda];  // Radial frequency [rad/s]
							       //printf("lambda = %e; omega = %e\n ", lambda[i_lambda],omega[i_lambda]);
							       //E_joules[i_lambda] = h_bar/omega[i_lambda];   // Wave energy [J] h_bar;//
							       //E_ev[i_lambda] = E_joules[i_lambda]/q;        // Wave energy [eV]
							       //printf("E = %e J; E = %e eV \n",,E_joules[i_lambda],E_ev[i_lambda]);
	}
	}
	
	if(strcmp(material,"SiC")==0) 
	{
	//Uniform spectrum: SAME SPECTRUM AS SiO2
	double initial_lambda = 5.e-6;//5.e-6;
	double final_lambda = 25.e-6; //25.e-6;
	initial = initial_lambda;
	final = final_lambda;
	double step_lambda = (final-initial)/(const_N_omega-1); // linspace in C: https://stackoverflow.com/questions/60695284/linearly-spaced-array-in-c
	for(int i_lambda = 0; i_lambda < const_N_omega; i_lambda++)
	{
		lambda[i_lambda]= initial + i_lambda*step_lambda; // Wavelength [m]
		omega[i_lambda] = 2.*pi*c_0/lambda[i_lambda];  // Radial frequency [rad/s]
							       //printf("lambda = %e; omega = %e\n ", lambda[i_lambda],omega[i_lambda]);
							       //E_joules[i_lambda] = h_bar/omega[i_lambda];   // Wave energy [J] h_bar;//
							       //E_ev[i_lambda] = E_joules[i_lambda]/q;        // Wave energy [eV]
							       //printf("E = %e J; E = %e eV \n",,E_joules[i_lambda],E_ev[i_lambda]);
	}
	}
	
	if(strcmp(material,"u-SiC")==0) 
	{
	//Uniform spectrum:  SAME SPECTRUM AS SiO2
	double initial_lambda = 5.e-6;//5.e-6;
	double final_lambda = 25.e-6; //25.e-6;
	initial = initial_lambda;
	final = final_lambda;
	double step_lambda = (final-initial)/(const_N_omega-1); // linspace in C: https://stackoverflow.com/questions/60695284/linearly-spaced-array-in-c
	for(int i_lambda = 0; i_lambda < const_N_omega; i_lambda++)
	{
		lambda[i_lambda]= initial + i_lambda*step_lambda; // Wavelength [m]
		omega[i_lambda] = 2.*pi*c_0/lambda[i_lambda];  // Radial frequency [rad/s]
							       //printf("lambda = %e; omega = %e\n ", lambda[i_lambda],omega[i_lambda]);
							       //E_joules[i_lambda] = h_bar/omega[i_lambda];   // Wave energy [J] h_bar;//
							       //E_ev[i_lambda] = E_joules[i_lambda]/q;        // Wave energy [eV]
							       //printf("E = %e J; E = %e eV \n",,E_joules[i_lambda],E_ev[i_lambda]);
	}
	}
	
	
	else if(strcmp(material,"SiN")==0) //cannot compare strings in C with ==; source: https://ideone.com/BrFA00
	{
	//Uniform spectrum: 
	double initial_omega = 2.e13;//5.e-6;
	double final_omega = 3.e14; //25.e-6;
	initial = initial_omega;
	final = final_omega;
	//even though we refer to step_lambda, for SiN we use radial frequency omega
	double step_omega = (final-initial)/(const_N_omega-1); // linspace in C: https://stackoverflow.com/questions/60695284/linearly-spaced-array-in-c
	for(int i_omega = 0; i_omega < const_N_omega; i_omega++)
	{
		//lambda[i_lambda]= initial + i_lambda*step_lambda; // Wavelength [m]
		omega[i_omega] = initial + i_omega*step_omega;  // Radial frequency [rad/s]
							       //printf("lambda = %e; omega = %e\n ", lambda[i_lambda],omega[i_lambda]);
							       //E_joules[i_lambda] = h_bar/omega[i_lambda];   // Wave energy [J] h_bar;//
							       //E_ev[i_lambda] = E_joules[i_lambda]/q;        // Wave energy [eV]
							       //printf("E = %e J; E = %e eV \n",,E_joules[i_lambda],E_ev[i_lambda]);
	}
	}
	

	
	// #################################################################
	// ################## FREQUENCY RANGE ANALYSIS #####################
	// #################################################################
	//Loop to analyze a range of desired frequencies
	printf("----- Spectrum range calculation -----\n");

	if(single_spectrum_analysis =='Y') omega_range=1;
	if(single_spectrum_analysis =='N') omega_range=const_N_omega;
	for (int i_omega = 0; i_omega < omega_range; i_omega++) // Frequency loop
	{

		if (i_omega==1)
		{
			single_analysis = 'n';
		}


		omega_value = omega[i_omega]; // omega definition is on line 212
		printf("%d) omega = %e. ", i_omega+1,omega_value);

		if(strcmp(material,"SiO2")==0) //cannot compare strings in C with ==; source: https://ideone.com/BrFA00
		{ 
			// Precision of this code is defined on the dielectric function, with 5 digits of precision.
			// Dielectric function for SiO2
			// Defined in Czapla and Narayanaswamy, JQSRT, 2019. DOI: https://doi.org/10.1016/j.jqsrt.2019.01.020 
			double epsilon_inf = 2.03843;
			double h_bar_omega_0[] = {0.05624, 0.09952, 0.13355};    // [eV]
			double S[] = {0.93752, 0.05050, 0.60642};                // [-]
			double Gamma[] = {0.09906, 0.05511,0.05246};
			double omega_0[3]; 
			double epsilon_term1[3];
			double epsilon_term2[3];
			double complex epsilon_num;
			double epsilon_denom; 
			epsilon = epsilon_inf + 0.*I ;
			for(int i_epsilon = 0; i_epsilon < 3; i_epsilon++)
			{
				omega_0[i_epsilon] = h_bar_omega_0[i_epsilon]*q/h_bar;
				epsilon_term1[i_epsilon] = pow(omega_value/omega_0[i_epsilon],2);
				epsilon_term2[i_epsilon] = Gamma[i_epsilon]*omega_value/omega_0[i_epsilon];
				epsilon_num = 0. + 0.*I;
				epsilon_num = S[i_epsilon]*(1-epsilon_term1[i_epsilon]) + S[i_epsilon]*epsilon_term2[i_epsilon]*I;
				epsilon_denom=0.;
				epsilon_denom = pow(1. - epsilon_term1[i_epsilon],2) + pow(epsilon_term2[i_epsilon],2); 
				epsilon += epsilon_num/epsilon_denom ;  
			}
		}

		else if(strcmp(material,"SiC")==0) 
		{  
			// Dielectric function for SiC -  Lorentz oscillator model:
			// Defined in Francoeur et al., PRB, 2011. DOI: https://doi.org/10.1103/PhysRevB.84.075436           
			double epsilon_inf = 6.7;              // [-]
			double Gamma = 8.966e11; 		//[1/s]  Drude relaxation time
			double omega_TO = 1.494e14;		//[rad/s]		
			double omega_LO = 1.825e14;		//[rad/s]
			epsilon = epsilon_inf*(pow(omega_value,2)-pow(omega_LO,2)+I*Gamma*omega_value)/(pow(omega_value,2)-pow(omega_TO,2)+I*Gamma*omega_value) ; // omega_TO omega_LO Gamma epsilon_inf
		}
		else if(strcmp(material,"u-SiC")==0) 
		{  
			// Dielectric function for u-SiC -  Lorentz oscillator model:
			// Defined in St-Gelais et al., Nature nanotechnology, 2016. DOI: https://doi.org/10.1038/nnano.2016.20    
			double to_rad_s = 100*2*pi*c_0;       
			double epsilon_inf = 8;              // [-]
			double Gamma = 20*to_rad_s; 		//[1/s]  Drude relaxation time
			double omega_TO = 789*to_rad_s;		//[rad/s]		
			double omega_LO = 956*to_rad_s;		//[rad/s]
			epsilon = epsilon_inf*(pow(omega_value,2)-pow(omega_LO,2)+I*Gamma*omega_value)/(pow(omega_value,2)-pow(omega_TO,2)+I*Gamma*omega_value) ; // omega_TO omega_LO Gamma epsilon_inf
		}
		
		else if(strcmp(material,"SiN")==0) 
		{  
			// Dielectric function for SiN -  Lorentz oscillator model:
			// Defined in Cataldo et al, Optics letters 2012. DOI: https://doi.org/10.1364/OL.37.004200
			int M=5;
			double complex epsilon_j[] = {7.582+ 0*I, 6.754+0.3759*I, 6.601+ 0.0041*I, 5.430+0.1179*I, 4.601+0.2073*I, 4.562 + 0.0124*I};
			double complex epsilon_inf = epsilon_j[5];
			double omega_T[] = {13.913, 15.053, 24.521, 26.440, 31.724};
			double Gamma[] = {5.810, 6.436, 2.751, 3.482, 5.948};
			double alpha[] = {0.0001, 0.3427, 0.0006, 0.0002, 0.0080};
			double omega_T_converted, Gamma_converted;
			double complex delta_epsilon, num_1, denom_1,Gamma_prime,num_2, denom_2;
		
			double complex summation = 0. + 0.*I;
			for(int i_epsilon = 0; i_epsilon < M; i_epsilon++)
			{				
				omega_T_converted = omega_T[i_epsilon]*(2*pi)*(1.e12);
				Gamma_converted =Gamma[i_epsilon]*(2*pi)*(1.e12);
				num_1 = pow(omega_T_converted,2)- pow(omega_value,2);
				denom_1 = omega_value*Gamma_converted;
    				
    				Gamma_prime = Gamma_converted*cexp(-alpha[i_epsilon]*pow(num_1/denom_1,2));
    								
				delta_epsilon = epsilon_j[i_epsilon] - epsilon_j[i_epsilon+1]; 
				
				num_2 = delta_epsilon*pow(omega_T_converted,2);
				denom_2 = pow(omega_T_converted,2) - pow(omega_value,2) - omega_value*Gamma_prime*I;
				summation += num_2/denom_2;
				
			}
			epsilon = epsilon_inf + summation;

		}
		
		
		//printf("Epsilon  = %e + i %e \n",creal(epsilon),cimag(epsilon));

		//printf("##################### \n Build A = I - k_0^2*G^0*alpha_0 \n##################### \n");

		double k_0=k_0_function(omega_value, epsilon_0, mu_0) ; //wave vector in free space
						       //printf("k_0 = %e; ",k_0);

		double k=k_function(omega_value, epsilon_ref, epsilon_0, mu_0); //wave vector in reference medium
						  //printf("k = %e; ",k);

						  //printf("alpha_0 = [ ");
		for (int i_alpha = 0; i_alpha < tot_sub_vol; i_alpha++)
		{
			alpha_0[i_alpha] = delta_V_vector[i_alpha]*(epsilon - epsilon_ref); //Bare polarizability [m^3]
											    //  printf("%e + i %e ;", creal(alpha_0[i_alpha]), cimag(alpha_0[i_alpha]));
		}
		//printf("] \n");

		//################## FREE-SPACE GREEN'S FUNCTION AND INTERACTION A MATRIX ##################### 
		// Fill terms for G^0: 
		//printf("Each G_ij^0 is computed considering its analytical solution \n");

		//Definitions for G^0: 

		double (*r)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*r));  // Complete discretized lattice [m]
		double (*abs_r_ij)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*abs_r_ij)); 
		double (*unit_r_ij)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*unit_r_ij));  
		double (*transpose)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*transpose));
		double complex (*unit_conj_r_ij)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*unit_conj_r_ij)); 
		double complex (*r_i_j_outer_r_i_j)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*r_i_j_outer_r_i_j));  

		for (int i=0; i<tot_sub_vol; i++)
		{   
			for (int j=0; j<tot_sub_vol; j++)
			{
				for (int i_alpha=0; i_alpha<3; i_alpha++)
				{
					//r[i][j][i_alpha]=R2[i_alpha] - R1[i_alpha]; // for 1 dipole p/ sphere case // G^0 function matlab code, line 34 
					r[i][j][i_alpha]= R[i][i_alpha] - R[j][i_alpha] ; //general case 
											  //printf("%e ; ",r[i][j][i_alpha]);
					abs_r_ij[i][j][i_alpha] = fabs(r[i][j][i_alpha]);
					//abs_r_ij[i][j][i_alpha] = fabs(r[i][j][i_alpha]);   // absolute value for double data type
					//printf("r%d%d = %e and abs(r%d%d) = %e ; ",i,j,r[i][j][i_alpha],i,j,abs_r_ij[i][j][i_alpha]);
					//printf("\n");
					//modulo_r_i_j[i][j] += pow(abs_r_ij[i][j][i_alpha],2);
					modulo_r_i_j[i][j] += pow(abs_r_ij[i][j][i_alpha],2);
				}
				modulo_r_i_j[i][j] = sqrt(modulo_r_i_j[i][j]);
			}
		}

		//printf("]\n");
		//printf("Abs r_i_j = %e ; \n",modulo_r_i_j[0][1]); // value used on free space green function solution

		for (int i_i = 0; i_i < tot_sub_vol; i_i++)
		{ 
			for (int i_j = 0; i_j < tot_sub_vol; i_j++)
			{ 
				for (int i_alpha = 0; i_alpha < 3; i_alpha++)
				{  
					unit_r_ij[i_i][i_j][i_alpha] = r[i_i][i_j][i_alpha]/modulo_r_i_j[i_i][i_j]; // ˆr -- unit distance 
					transpose[i_j][i_i][i_alpha] = r[i_i][i_j][i_alpha]/modulo_r_i_j[i_i][i_j]; // https://www.programiz.com/c-programming/examples/matrix-transpose        
														    //printf("%d%d: ",i_alpha,i_alpha);//,j_alpha);
					unit_conj_r_ij[i_i][i_j][i_alpha] = conj(transpose[i_j][i_i][i_alpha]);// r^+ Conjugate transpose unit distance 
													       //printf("unit= %e ; ", unit_r_ij[i_i][i_j][i_alpha]);
													       //printf("unit= %e + i%e; conj= %e + i%e; ", creal(unit_r_ij[i_i][i_j][i_alpha]),cimag(unit_r_ij[i_i][i_j][i_alpha]),creal(unit_conj_r_ij[i_i][i_j][i_alpha]),cimag(unit_conj_r_ij[i_i][i_j][i_alpha]));
													       //printf("\n");  
				}
			}
		}

		for (int i_i = 0; i_i < tot_sub_vol; i_i++)
		{ 
			for (int i_j = 0; i_j < tot_sub_vol; i_j++)
			{ 
				//printf("%d%d:\n",i_i,i_j);
				//if (i_i!=i_j)
				//{
				for (int i_alpha = 0; i_alpha < 3; i_alpha++)
				{  
					for (int j_alpha = 0; j_alpha < 3; j_alpha++)
					{ 
						//printf("%d%d: ",i_alpha,j_alpha);
						r_i_j_outer_r_i_j[i_i][i_j][i_alpha][j_alpha] = unit_r_ij[i_i][i_j][i_alpha]*unit_conj_r_ij[i_i][i_j][j_alpha];//*unit_conj_r_ij[i_i][i_j][j_alpha]; // [j_alpha] outer product
																			       //printf("outer= %e +i%e; ", creal(r_i_j_outer_r_i_j[i_i][i_j][i_alpha][j_alpha]),cimag(r_i_j_outer_r_i_j[i_i][i_j][i_alpha][j_alpha]));    //
					}
					//printf("\n");  
				}    
				//} 
			}
		}
		
		// ################### MATRICES STRUCTURE LOOPS ###########################
		// 3N X 3N Matrices structure loops for G^0 and A:
		//printf("  ----- Build G^0 and A matrices -----\n");

		//G^0_ij when i=j:
		double (*a_j) = malloc(sizeof *a_j *tot_sub_vol);  //double a_j[tot_sub_vol];
		double (*part1ii) = malloc(sizeof *part1ii *tot_sub_vol); //double part1ii[tot_sub_vol];
		double complex (*part2ii) = malloc(sizeof *part2ii *tot_sub_vol); //double complex part2ii[tot_sub_vol]; 
		double complex (*part2iiexp) = malloc(sizeof *part2iiexp *tot_sub_vol); //double complex part2iiexp[tot_sub_vol]; 
		double complex (*part3ii) = malloc(sizeof *part3ii *tot_sub_vol); //double complex part3ii[tot_sub_vol];

		//G^0_ij when i!=j:
		double complex (*part1ij) = malloc(sizeof *part1ij *tot_sub_vol); //double complex part1ij[tot_sub_vol];
		double complex (*part1aij) = malloc(sizeof *part1aij *tot_sub_vol); //double complex part1aij[tot_sub_vol];
		double complex (*part1aijexp) = malloc(sizeof *part1aijexp *tot_sub_vol); //double complex part1aijexp[tot_sub_vol];
		double complex (*part2ij) = malloc(sizeof *part2ij *tot_sub_vol); //double complex part2ij[tot_sub_vol];
		double complex (*part3ij) = malloc(sizeof *part3ij *tot_sub_vol); //double complex part3ij[tot_sub_vol];

		double (*eyeG_0)[3] = calloc(3, sizeof(*eyeG_0)); //double eyeG_0[3][3];
		double (*eyeA)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*eyeA)); //double eyeA[tot_sub_vol][tot_sub_vol][3][3];

		// Linear system AG=G^0 
		double complex (*G_0)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*G_0)); //double complex G_0[tot_sub_vol][tot_sub_vol][3][3]; //[tot_sub_vol][tot_sub_vol][3][3]; [size_mat][size_mat]
		double complex (*A)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*A)); //double complex A[tot_sub_vol][tot_sub_vol][3][3];

		// eq. 25 from Lindsay's paper 
		for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
		{
			for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
			{
				if (ig_0!=jg_0) // if i!=j:
				{
					part1aij[ig_0] = 0.+ k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0]*I; // com i term 
					part1aijexp[ig_0]= cexp(k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0]*I);
					//printf("exp=%e + i %e \n",creal(part1aijexp[ig_0]),cimag(part1aijexp[ig_0])); // 0.990461 + i*0.13779 
					part1ij[ig_0] = part1aijexp[ig_0]/(4.*pi*modulo_r_i_j[ig_0][jg_0]);
					//printf("Part 1: %e + i %e \n",creal(part1ij[ig_0]),cimag(part1ij[ig_0])); // bate com Mathematica e Lindsay
					denom1 = epsilon_ref*pow(k_0*modulo_r_i_j[ig_0][jg_0],2);
					denom2 = k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0];
					part2ij[ig_0] = (1. - 1./denom1 + 1.*I/denom2 ) ;
					part3ij[ig_0] = (1. - 3./denom1 + 3.*I/denom2) ;
					for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
					{
						for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
						{
							if (i_subG_0 == j_subG_0)
							{
								eyeG_0[i_subG_0][j_subG_0] = 1.;    // 3x3 Identity matrix for G^0:
								eyeA[ig_0][jg_0][i_subG_0][j_subG_0] = 0.; // 3Nx3N identity matrix for A:
													   //printf("%e ", eyeA[ig_0][jg_0][i_subG_0][j_subG_0]);
													   //G_0[ig_0][jg_0][i_subG_0][j_subG_0] = part1ij[ig_0]*((part2ij[ig_0]*eyeG_0[i_subG_0][j_subG_0])-(part3ij[ig_0]*r_i_j_outer_r_i_j[ig_0][jg_0][j_subG_0]));  //*r_i_j_outer_r_i_j[ig_0][jg_0]
							}
							else
							{
								eyeG_0[i_subG_0][j_subG_0] = 0.;     // 3x3 Identity matrix for G^0:
								eyeA[ig_0][jg_0][i_subG_0][j_subG_0] = 0.; // 3Nx3N identity matrix for A: 
													   //printf("%e ", eyeA[ig_0][jg_0][i_subG_0][j_subG_0]);
													   //G_0[ig_0][jg_0][i_subG_0][j_subG_0] = part1ij[ig_0]*((part2ij[ig_0]*eyeG_0[i_subG_0][j_subG_0])-(part3ij[ig_0]*r_i_j_outer_r_i_j[ig_0][jg_0][j_subG_0])); //0. ;
							}
							G_0[ig_0][jg_0][i_subG_0][j_subG_0] = part1ij[ig_0]*((part2ij[ig_0]*eyeG_0[i_subG_0][j_subG_0])-(part3ij[ig_0]*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  //*r_i_j_outer_r_i_j[ig_0][jg_0]
																											     //printf("%e + i %e ; ", creal(G_0[ig_0][jg_0][i_subG_0][j_subG_0]),cimag(G_0[ig_0][jg_0][i_subG_0][j_subG_0]));//[i_subG_0][j_subG_0]
							A[ig_0][jg_0][i_subG_0][j_subG_0] = eyeA[ig_0][jg_0][i_subG_0][j_subG_0] - pow(k_0,2)*alpha_0[ig_0]*G_0[ig_0][jg_0][i_subG_0][j_subG_0]; 
							//printf("%e + i %e ; ", creal(A[ig_0][jg_0][i_subG_0][j_subG_0]),cimag(A[ig_0][jg_0][i_subG_0][j_subG_0]));   
						}    
					}
				}    
			}    
		} //end ig_0 

		// eq. 26 from Lindsay's paper: 


		for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
		{
			a_j[ig_0] = a_j_function(delta_V_vector[ig_0], pi);
			//printf("a_j = %e; ",a_j[ig_0]);
			part1ii[ig_0] = 1./(3.*delta_V_vector[ig_0]*epsilon_ref*pow(k_0,2));
			//printf("Cte 1: %e \n",part1ii[ig_0]);
			part2ii[ig_0] = a_j[ig_0]*k_0*sqrt(epsilon_ref)*I; // com i term 
									   //printf("Cte 2: %e \n",cimag(part2ii[ig_0]));
			part2iiexp[ig_0] = cexp(0. + a_j[ig_0]*k_0*sqrt(epsilon_ref)*I); // complex exp: cexp()
											 //printf("exp=%e + i %e \n",creal(part2iiexp[ig_0]),cimag(part2iiexp[ig_0])); // 0.998027 + i 0.0627906
											 // part3ii is inside brackets
			part3ii[ig_0] = part2iiexp[ig_0]*(1-part2ii[ig_0]) - 1. ;
			//printf("G_%d%d^0 = \n", ig_0+1,jg_0+1);
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
				{ 
					if (ig_0==jg_0) // if i=j:
					{
						for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
						{
							if (i_subG_0 == j_subG_0)
							{
								eyeG_0[i_subG_0][j_subG_0] = 1.;     // 3x3 Identity matrix for G^0:
								eyeA[ig_0][jg_0][i_subG_0][j_subG_0] = 1.; // 3Nx3N identity matrix for A:
													   //printf("%e ", eyeA[ig_0][jg_0][i_subG_0][j_subG_0]);
							}
							else
							{
								eyeG_0[i_subG_0][j_subG_0] = 0.;     // 3x3 Identity matrix for G^0:
								eyeA[ig_0][jg_0][i_subG_0][j_subG_0] = 0.; // 3Nx3N identity matrix for A:
													   //printf("%e ", eyeA[ig_0][jg_0][i_subG_0][j_subG_0]);
							}
							G_0[ig_0][jg_0][i_subG_0][j_subG_0] = eyeG_0[i_subG_0][j_subG_0]*part1ii[ig_0]*(2.*part3ii[ig_0]-1.); 
							//printf("%e + i %e ; ", creal(G_0[ig_0][jg_0][i_subG_0][j_subG_0]),cimag(G_0[ig_0][jg_0][i_subG_0][j_subG_0]));
							A[ig_0][jg_0][i_subG_0][j_subG_0] = eyeA[ig_0][jg_0][i_subG_0][j_subG_0] - pow(k_0,2)*alpha_0[ig_0]*G_0[ig_0][jg_0][i_subG_0][j_subG_0]; 
							//printf("%e + i %e ; ", creal(A[ig_0][jg_0][i_subG_0][j_subG_0]),cimag(A[ig_0][jg_0][i_subG_0][j_subG_0]));
						}
					} 
				}   //end jg_0   
				    //  printf("\n"); 
			} //end i_subG_0     
		} //end ig_0 


		free(eyeG_0);
		//free(eyeA);
		free(r);
		free(abs_r_ij);
		free(unit_r_ij);
		free(transpose);
		free(unit_conj_r_ij);
		free(r_i_j_outer_r_i_j);

		free(a_j);
		free(part1ii);
		free(part2ii);
		free(part2iiexp);
		free(part3ii);
		free(part1ij);
		free(part1aij);
		free(part1aijexp);
		free(part2ij);
		free(part3ij);

		memset(modulo_r_i_j, 0, sizeof modulo_r_i_j);
		//free(modulo_r_i_j);

/*		
		//printf("##################### \n EXPORT DATA \n##################### \n");
		//char text[1024];
		//mkdir(array[i_omega], 0700);  // Create directory: https://stackoverflow.com/questions/7430248/creating-a-new-directory-in-c

		if(save_A_matrix =='Y'|save_G0_matrix =='Y'|save_SGF_matrix =='Y'){
			sprintf(matrices_folder, "%s/matrices_folder", results_folder); // How to store words in an array in C? https://www.geeksforgeeks.org/how-to-store-words-in-an-array-in-c/
			create_folder(matrices_folder);

			// Folder of frequencies
			sprintf(frequency_folder, "%s/omega_%d", matrices_folder,i_omega+1); // How to store words in an array in C? https://www.geeksforgeeks.org/how-to-store-words-in-an-array-in-c/
			create_folder(frequency_folder);


			//printf("----- Export data -----\n");
			if(save_A_matrix =='Y'){
				FILE * Amatrix; // A matrix .csv file
				char dirPathAFileName[260]; // https://stackoverflow.com/questions/46612504/creating-directories-and-files-using-a-loop-in-c
				sprintf(dirPathAFileName, "%s/%s", frequency_folder, "A_matrix.csv"); // path where the file is stored
				Amatrix = fopen(dirPathAFileName,"w");//array[i_omega]"/A_matrix.txt" // EN: https://www.tutorialspoint.com/c_standard_library/c_function_fopen.htm
				for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
				{
					for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
					{
						for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
						{
							for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
							{
								fprintf(Amatrix,"%e + i %e \t ", creal(A[ig_0][jg_0][i_subG_0][j_subG_0]),cimag(A[ig_0][jg_0][i_subG_0][j_subG_0])); // PS: note that the function is fprintf not printf    
							}    
						}
						fprintf(Amatrix,"\n"); 
					}        
				}
				fclose(Amatrix);
			}
			if(save_G0_matrix =='Y'){
				FILE * G0matrix; // G^0 matrix .csv file
				char dirPathG0FileName[260];
				sprintf(dirPathG0FileName, "%s/%s", frequency_folder,"G0_matrix.csv"); // path where the file is stored
				G0matrix =fopen(dirPathG0FileName,"w"); // PT: https://terminalroot.com.br/2014/12/linguagem-c-utilizando-as-funcoes-fopen.html
				for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
				{
					for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
					{
						for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
						{
							for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
							{
								fprintf(G0matrix,"%e + i %e \t ", creal(G_0[ig_0][jg_0][i_subG_0][j_subG_0]),cimag(G_0[ig_0][jg_0][i_subG_0][j_subG_0]));
							}    
						}
						fprintf(G0matrix,"\n");
					}        
				} 
				fclose(G0matrix); 
			}  

		} //end if export data 
*/
		
		//printf("##################### \n SOLVE LINEAR SYSTEM AG=G^0 \n##################### \n");
		//printf("##################### \n LAPACK/LAPACKE ZGELS ROUTINE \n##################### \n");
		//Description of ZGELS: https://extras.csc.fi/math/nag/mark21/pdf/F08/f08anf.pdf
		//printf("  ----- Solve linear system AG=G^0 -----\n");

		double complex (*G_sys)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys));
		//double complex (*G_sys)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*G_sys)); //double complex G_sys[tot_sub_vol][tot_sub_vol][3][3];

		if(solution =='D')
		{
			printf("Direct inversion status: ");
			//double complex Alapack[lda*n], blapack[ldb*nrhs], work[lwork]; // based on example from https://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=506&p=1692&hilit=zgels#p1692
			double complex (*A_direct) = malloc(sizeof *A_direct *lda*n); //double complex Alapack[lda*n];
			double complex (*b_direct) = malloc(sizeof *b_direct *ldb*nrhs); //double complex blapack[ldb*nrhs];
			//double complex (*work) = malloc(sizeof *work *lwork); //double complex work[lwork];

			ipack=0; 
			for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
			{
				for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
				{
					for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
					{
						for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
						{
							A_direct[ipack] = A[ig_0][jg_0][i_subG_0][j_subG_0];
							b_direct[ipack]= G_0[ig_0][jg_0][i_subG_0][j_subG_0];
							ipack = ipack + 1;
						}    
					}
				}        
			}   

			//printf("  --------- ZGELS -----\n");
			//double complex G_sys[tot_sub_vol][tot_sub_vol][3][3];
			gpack=0;
			// F08ANF (ZGELS) solves linear least-squares problems using a QR or LQ factorization of A
			info = LAPACKE_zgels(LAPACK_ROW_MAJOR,'N',m,n,nrhs,A_direct,lda,b_direct,ldb); 
			free(A_direct); 
			
			for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
			{
				for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
				{
					ig_0_2d = (3*ig_0 + i_subG_0);
					for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
					{
						for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
						{
							jg_0_2d = (3*jg_0 + j_subG_0);
							G_sys[ig_0_2d][jg_0_2d] = b_direct[gpack];
							//G_sys[ig_0][jg_0][i_subG_0][j_subG_0] = blapack[gpack];
							//printf("%e + i %e ; ", creal(G[ig_0][jg_0][i_subG_0][j_subG_0]),cimag(G[ig_0][jg_0][i_subG_0][j_subG_0]));
							gpack=gpack + 1;
						}    
					}
					//printf("\n");
				}        
			} 

			free(b_direct); 
			//free(work);
			printf("concluded\n");
		}
		free(A); //A[3Nx3N] matrix is not used in the iterative solver, Amm[3x3] is used instead.

		
		if(solution =='I')
		{ 
			printf("Iterative status:\n m= ");

			double complex (*G_sys_old)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys_old)); 
			double complex (*G_sys_new)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys_old)); 
			//double (*eyeA_2d)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*eyeA_2d)); // DO WE NEED THIS???
			double complex (*A_2d)[3] = calloc(3, sizeof(*A_2d)); // Amm
			double complex (*A_iterative) = malloc(sizeof *A_iterative *3*3); // direct inversion when i=m;
			double complex (*b_iterative) = malloc(sizeof *b_iterative *3*3); //direct inversion when i=m;
   					
			double eye_iter[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}}; // 3-by-3 unit matrix used in iterative solver
						
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// Calculate background medium Green's function 
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
			
			// this reshapes 2 4D matrices to 2 2D matrices
			// G0 and eyeA are 4D
			// G_sys_old and eyeA_2d are the respective 2D matrices
			//matrix_reshape(3, tot_sub_vol, G_sys_old, eyeA_2d, G_0, eyeA);
			matrix_reshape(3, tot_sub_vol, G_sys_old, G_0);

    	//	print_matrix(3*tot_sub_vol, 3*tot_sub_vol, G_sys_old);
		
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// Calculate system Green's function 
			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			// First, solve ii = mm system of equations.
			for (int mm = 0; mm < tot_sub_vol; mm++) //tot_sub_vol
			{
				printf("%d - ",mm+1);
				mm_2d =0;
        
				epsilon_s = (epsilon - epsilon_ref); // Scattering dielectric function
				
				A2d_solver(epsilon_s, mm, tot_sub_vol, eye_iter, delta_V_vector[mm], G_sys_old, A_2d, k);	
				/*
				// %%%%%%%%%%%%%%%%%% Manual inversion of Amm %%%%%%%%%%%%%%%%%%%% 
					det =  (Amm[0][0] * (Amm[1][1] * Amm[2][2] - Amm[2][1] * Amm[1][2])) - (Amm[0][1] * (Amm[1][0] * Amm[2][2] - Amm[1][2] * Amm[2][0])) + (Amm[0][2] * (Amm[1][0] * Amm[2][1] - Amm[1][1] * Amm[2][0]));

 					// computes the inverse of Amm
					Amm_inv[0][0] = (Amm[1][1] * Amm[2][2] - Amm[2][1] * Amm[1][2]) / det;
					Amm_inv[0][1] = (Amm[0][2] * Amm[2][1] - Amm[0][1] * Amm[2][2]) / det;
					Amm_inv[0][2] = (Amm[0][1] * Amm[1][2] - Amm[0][2] * Amm[1][1]) / det;
					Amm_inv[1][0] = (Amm[1][2] * Amm[2][0] - Amm[1][0] * Amm[2][2]) / det;
					Amm_inv[1][1] = (Amm[0][0] * Amm[2][2] - Amm[0][2] * Amm[2][0]) / det;
					Amm_inv[1][2] = (Amm[1][0] * Amm[0][2] - Amm[0][0] * Amm[1][2]) / det;
					Amm_inv[2][0] = (Amm[1][0] * Amm[2][1] - Amm[2][0] * Amm[1][1]) / det;
					Amm_inv[2][1] = (Amm[2][0] * Amm[0][1] - Amm[0][0] * Amm[2][1]) / det;
					Amm_inv[2][2] = (Amm[0][0] * Amm[1][1] - Amm[1][0] * Amm[0][1]) / det;
				// %%%%%%%%%%%%%%%%%% end Manual inversion of Amm %%%%%%%%%%%%%%%%%%	
				*/	
								
				for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) // Only loop through remaining perturbations
				{

					 A_b_iterative_populator(tot_sub_vol, A_iterative, b_iterative, A_2d, G_sys_old, mm, jg_0);

					// %%%%%%%%%%% G_new using Linear inversion using LAPACK %%%%%%%%%%%%%%%%%
					info = LAPACKE_zgels(LAPACK_ROW_MAJOR,'N',3,3,3,A_iterative,3,b_iterative,3); 
					gpack=0;     	

					for (int mm_sub = 0; mm_sub < 3; mm_sub++) // 3D coordinate positions
					{
						mm_2d = (3*mm + mm_sub);
						for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
						{
							jg_0_2d = (3*jg_0 + j_subG_0); // Set indices

							G_sys_new[mm_2d][jg_0_2d] = b_iterative[gpack]; // stores G^1_11
							//G_sys_new[mm][jg_0][i_subG_0][j_subG_0] = b1lapack[gpack]; // stores G^1_11

							gpack=gpack + 1;
						}  
					}
					// %%%%%%%%%%% End G_new using linear inversion using LAPACK %%%%%%%%%%%%%%%%%
					

					/*				
					// %%%%%%%%%%%%%%%%%% G_new using manual inversion of Amm %%%%%%%%%%%%%%%%%%%% 
										
					for (int mm_sub = 0; mm_sub < 3; mm_sub++) // 3D coordinate positions
					{
						mm_2d = (3*mm + mm_sub);
						for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
						{
							jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
							G_sys_new[mm_2d][jg_0_2d] = 0;
							for(int k_sub = 0; k_sub < 3; k_sub++)
							{
								k_sub_2d = (3*mm + k_sub); 
								G_sys_new[mm_2d][jg_0_2d] += Amm_inv[mm_sub][k_sub]*G_sys_old[k_sub_2d][jg_0_2d]; // stores G^1_11	
							}	
						}  
					}
					// %%%%%%%%%%%%%%%%%% end G_new using manual inversion of Amm %%%%%%%%%%%%%%%%%%%% 
					*/


				} // end jg_0					
				
				memset(A_2d, 0, sizeof *A_2d * 3);
				memset(A_iterative, 0, sizeof *A_iterative *3*3); //double complex Alapack[lda*n];
				memset(b_iterative, 0, sizeof *b_iterative *3*3);
								
				// Next, solve all systems of equations for ii not equal to mm		
				for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
				{ 
					//mm_sub_n = 0;  
					for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
					{
						ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
						for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol  ig_0 //lower triangular matrix
						{
							if(ig_0 != mm)
							{
								for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
								{
									jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
									G_sys_prod = 0.;
									//G_sys_prod[ig_0_2d][jg_0_2d] = 0.;
									for(int m_sub = 0;  m_sub < 3;  m_sub++)//loop for matricial multiplication
									{
										mm_2d = (3*mm + m_sub);
										G_sys_prod += G_sys_old[ig_0_2d][mm_2d]*G_sys_new[mm_2d][jg_0_2d];// 
										//G_sys_prod[ig_0_2d][jg_0_2d] += G_sys_old[ig_0_2d][mm_2d]*G_sys_new[mm_2d][jg_0_2d];//
									}

									G_sys_new[ig_0_2d][jg_0_2d] = G_sys_old[ig_0_2d][jg_0_2d] + pow(k,2)*alpha_0[mm]*G_sys_prod;
									//G_sys_new[ig_0_2d][jg_0_2d] = G_sys_old[ig_0_2d][jg_0_2d] + pow(k,2)*alpha_0[mm]*G_sys_prod[ig_0_2d][jg_0_2d];		

								} // j_subG_0            				
							} // if(ig_0 != mm) 
						} // jg_0   	
					}// i_subG_0   	                	
				} // ig_0    
								
				memcpy(G_sys_old,G_sys_new,3*tot_sub_vol*3*tot_sub_vol*sizeof(double complex)); // Update G_old = G_new for next iteration.
				

			}//end mm loop 
			
			printf("Final solution:\n");

			free(A_iterative); 
			free(b_iterative);
			free(A_2d);
			free(G_sys_old);
			
			memcpy(G_sys,G_sys_new,3*tot_sub_vol*3*tot_sub_vol*sizeof(double complex)); //Populate G_sys with G_new 
			free(G_sys_new);	   

		}

		free(G_0);
		
		/*
		if(save_SGF_matrix =='Y'){   
			printf("  --------- Export G_sys -----\n");

			FILE * G_sys_matrix; // G matrix .csv file
			char dirPathG_sys_FileName[260];
			sprintf(dirPathG_sys_FileName, "%s/%s", frequency_folder,"G_sys_matrix.csv"); // path where the file is stored
			G_sys_matrix =fopen(dirPathG_sys_FileName,"w"); // PT: https://terminalroot.com.br/2014/12/linguagem-c-utilizando-as-funcoes-fopen.html
			for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
			{
				for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
				{
					for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
					{
						for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
						{
							fprintf(G_sys_matrix,"%e + i %e \t ", creal(G_sys[ig_0][jg_0][i_subG_0][j_subG_0]),cimag(G_sys[ig_0][jg_0][i_subG_0][j_subG_0]));    
						}    
					}
					fprintf(G_sys_matrix,"\n");
				}        
			} 
			fclose(G_sys_matrix); 

		}
		*/
		
		// #################################################################
		//  Spectral and total conductance between bulk objects at temp. T
		// ###################  Thermal power dissipated ###################
		// #################################################################

		double complex (*trans_coeff)[tot_sub_vol] = calloc(tot_sub_vol, sizeof(*trans_coeff)); //double complex trans_coeff[tot_sub_vol][tot_sub_vol];//[3][3]
		
		double complex (*G_sys_cross)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*G_sys_cross)); 
		double complex (*transpose_G_sys)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*transpose_G_sys)); 
		double(*inner_sum) = malloc(sizeof * inner_sum * tot_sub_vol); //double inner_sum[jg_0];

		//printf("  ----- Spectral transmissivity, spectral conductance, and thermal power dissipated -----\n");
		int counter = 0;
		sum_trans_coeff[i_omega] = 0.;    

		for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
		{
			//printf("Theta(%e,%f) = %e\n",omega_value,T_vector[ig_0],theta[ig_0]);
			for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
			{
				for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
				{
					ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
					{
						jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
						// Extract one Green's function component: G_sys[ig_0][jg_0][i_subG_0]
						//transpose_G_sys[jg_0][ig_0][i_subG_0][j_subG_0] = G_sys[ig_0][jg_0][i_subG_0][j_subG_0];
						transpose_G_sys[jg_0][ig_0][i_subG_0][j_subG_0] = G_sys[ig_0_2d][jg_0_2d];
						G_sys_cross[ig_0][jg_0][i_subG_0][j_subG_0]= conj(transpose_G_sys[jg_0][ig_0][i_subG_0][j_subG_0]); // G_element'						
																		    //if (i_subG_0==i_subG_0) //
																		    //{
						//G_element[ig_0][jg_0]+=G_sys[ig_0][jg_0][i_subG_0][j_subG_0]*G_sys_cross[ig_0][jg_0][i_subG_0][j_subG_0]; //sum(diag((G_element*(G_element'))
						G_element[ig_0][jg_0]+=G_sys[ig_0_2d][jg_0_2d]*G_sys_cross[ig_0][jg_0][i_subG_0][j_subG_0]; //sum(diag((G_element*(G_element'))
																			  //}

					}    
				}
				// Transmissivity coefficient matrix tau(omega) for comparison with Czapla Mathematica output [dimensionless]
				//Matlab code: trans_coeff= 4*(k_0^4)*delta_V(i)*delta_V(j)*imag(epsilon(i))*imag(epsilon(j))*sum(diag((G_element*(G_element')))) 
				trans_coeff[ig_0][jg_0] = 4.*pow(k_0,4)*delta_V_vector[ig_0]*delta_V_vector[jg_0]*cimag(epsilon)*cimag(epsilon)*G_element[ig_0][jg_0] ;  
				//printf("%e ; ",trans_coeff[ig_0][jg_0]); // values match with Lindsay's results

				// Thermal power dissipated calculation, based on the matlab code (Using Tervo's Eq. 26)
				if (ig_0 != jg_0) 
				{
					inner_sum[jg_0] = (theta_function(omega_value, T_vector[jg_0], h_bar, k_b) - theta_function(omega_value, T_vector[ig_0], h_bar, k_b)) * trans_coeff[ig_0][jg_0]; //Q_omega_subvol_function(theta_1,theta_2, trans_coeff[ig_0][jg_0]);
				}
				else {
					inner_sum[jg_0] = 0;
				}

				//Trans_bulk: Transmission coefficient between two bulk objects
				// This function calculates the transmission coefficient between bulk objects given the transmission coefficient between every pair of dipoles for a given frequency.
				if(ig_0 < const_N_subvolumes_per_object && jg_0 >= const_N_subvolumes_per_object)// bulk 1 -> 2
				{
					sum_trans_coeff[i_omega] += trans_coeff[ig_0][jg_0];
					counter+=1;
				} 

				Q_omega_subvol[ig_0] += (1 / (2 * pi)) * inner_sum[jg_0]; // calculates the thermal power dissipated per subvolume

			} 
			Q_subvol[ig_0][i_omega] = Q_omega_subvol[ig_0];
			//printf("Q_subvol_%d_%d= %e \n",ig_0+1,i_omega+1, Q_subvol[ig_0][i_omega]); //matches matlab code

			if(ig_0 < const_N_subvolumes_per_object)// Thermal power dissipated was not verified yet!!!!
			{
				bulk=0;
				Q_omega_thermal_object[bulk][i_omega] += Q_omega_subvol[ig_0];
			}
			else if(const_N_subvolumes_per_object<=ig_0<tot_sub_vol)
			{
				bulk=1;
				Q_omega_thermal_object[bulk][i_omega] += Q_omega_subvol[ig_0];
			} 
			//printf("Q_bulk_%d_%d= %e \n",bulk+1, i_omega+1, Q_omega_thermal_object[bulk][i_omega]); //matches matlab code

		}       
		free(G_sys);

		free(trans_coeff);
		free(G_sys_cross);
		free(transpose_G_sys);

		free(inner_sum);


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

		double (*dtheta_dT) = malloc(sizeof *dtheta_dT *N_Tcalc); //double dtheta_dT[N_Tcalc]; // function used to calculate conductance, modified for several temperatures

		for (int iTcalc = 0; iTcalc < N_Tcalc; iTcalc++) // EDIT VALUE :: change 1 to N_Tcalc for the temperature loop
		{
			// printf("T_calc = %e; \n",Tcalc_vector[iTcalc]);
			//dtheta_dT[0] = dtheta_dT_function(omega_value,T_calc); 
			dtheta_dT[iTcalc] = dtheta_dT_function(omega_value,Tcalc_vector[iTcalc], h_bar, k_b); 
			//printf("dT = %e; \n",dtheta_dT);            
			//G_12_omega_SGF[i_omega][0] = dtheta_dT[0]*sum_trans_coeff[i_omega];
			G_12_omega_SGF[i_omega][iTcalc] = dtheta_dT[iTcalc]*sum_trans_coeff[i_omega];

			//G_12_lambda_SGF[i_omega] = dtheta_dT*sum_trans_coeff_lambda[i_omega]; // Spectral conductance (function of frequency) [W/(K*m)]


			if(save_spectral_conductance =='Y'){
				{
					FILE * spectral_conductance; //append
					char dirPathSpectral_cond_FileName[260];
					sprintf(dirPathSpectral_cond_FileName, "%s/spectral_conductance_%eK.csv", results_folder, Tcalc_vector[iTcalc]); // path where the file is stored
					if(i_omega == 0) spectral_conductance =fopen(dirPathSpectral_cond_FileName,"w"); //write
					else spectral_conductance = fopen(dirPathSpectral_cond_FileName, "a"); //append
													       //fprintf(spectral_conductance,"%e ; %e\n",omega_value,G_12_omega_SGF[i_omega]); 
					fprintf(spectral_conductance,"%e ; %e\n",omega_value,G_12_omega_SGF[i_omega][iTcalc]); 
					fclose(spectral_conductance);
				}
			} // end if save_spectral_conductance

		} // END FOR T_calc LOOP
		free(dtheta_dT);

		/*
		//double omega_SGF = omega_value;                                         // [rad/s]
		//double lambda_SGF = lambda[i_omega];                                       // [m]
		//double omega_eV = E_eV;                                         // [eV]
		//double Trans_12_omega_SGF;// = (Trans_12_omega);               // Transmissivity coefficient tau(omega) for comparison with Czapla Mathematica output [dimensionless]
		//double G_12_omega_SGF = dtheta_dT*Trans_12_omega_SGF;                  // Spectral conductance (function of frequency) [W/(K*(rad/s))]
		//double Trans_12_lambda_SGF = Trans_12_omega_SGF*c_0/(pow(lambda[i_omega],2));  // Transmissivity coefficient, tau(lambda) [1/(m*s)]
		//double G_12_lambda_SGF = dtheta_dT*Trans_12_lambda_SGF;                // Spectral conductance (function of frequency) [W/(K*m)]
		//printf("G_12: %e\n",G_12_omega_SGF);
		*/

		// #################################################################
		// #############  Clear values for next frequency ##################
		// #################################################################

		memset(modulo_r_i_j, 0, sizeof *modulo_r_i_j * tot_sub_vol); //modulo_r_i_j
		memset(G_element, 0, sizeof *G_element * tot_sub_vol); //G_element
		memset(sum_trans_coeff, 0, sizeof sum_trans_coeff);
		memset(Q_omega_subvol, 0, sizeof* Q_omega_subvol* tot_sub_vol); 
		//memset(Q_omega_thermal_object,0, sizeof Q_omega_thermal_object * const_N_bulk_objects * const_N_omega); // //double Q_omega_thermal_object[N_bulk_objects][const_N_omega]

	} // END OMEGA VALUE LOOP FOR FREQUENCY RANGE ANALYSIS, loop started in line 481

	free(R);
	free(sum_trans_coeff);
	free(modulo_r_i_j);
	free(G_element);
	//free(modulo_r_i_j);
	//free(G_element);

	free(alpha_0);

	double (*Total_conductance) = malloc(sizeof *Total_conductance *N_Tcalc); //double Total_conductance[N_Tcalc];
	double(*Q_tot_thermal_object) = malloc(sizeof * Q_tot_thermal_object * const_N_bulk_objects); //double Q_tot_thermal_object[const_N_bulk_objects];           
	double (*trapz_Q) = malloc(sizeof *trapz_Q *tot_sub_vol); // Definition for trapezoidal integration. Used in total power dissipated //double Total_conductance[N_Tcalc];
								  // int i_frequency = 0;

	trapz=0.;
	if( const_N_omega > 1)
	{
		printf("\nEnd of frequency loop\n");
		//printf("----- Total conductance -----\n");


		// implementation of trapezoidal numerical integration in C // https://stackoverflow.com/questions/25124569/implementation-of-trapezoidal-numerical-integration-in-c
		double sum_conductance = 0.;
		double step = 0.;
		//double sum_Q_subvol =0.;
		//double sum_Q_bulk =0.;
		double sum_Q = 0.;

		for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
		{
			trapz_Q[ig_0] = 0.;
			for (int i=1; i < const_N_omega; ++i)
			{
				step = omega[i] - omega[i-1];
				trapz_Q[ig_0] +=  ((Q_subvol[ig_0][i]+Q_subvol[ig_0][i-1])/2)*step;
			}

			//printf("Power dissipated per subvolume %d = %e\n", ig_0+1,trapz_Q[ig_0]);

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
		//Q_tot_thermal_object[0] = (trapz_Q)/(2.*pi);
		//printf("Q_tot_thermal_object[1] = %e; Q_tot_thermal_object[2] = %e\n", Q_tot_thermal_object[0],Q_tot_thermal_object[1]);



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
									   //G_12_total_SGF_from_lambda = fabs(trapz(lambda[i_omega], G_12_lambda_SGF));          // Total conductance [W/K]

									   //printf("Total conductance at temperature %e = %e \n",T_calc,G_12_total_SGF_from_omega);
									   //Total_conductance[0]=G_12_total_SGF_from_omega;
			Total_conductance[iTcalc]=G_12_total_SGF_from_omega;

		} // END FOR T_calc LOOP

		printf("Total conductance at %e K= %e \n",Tcalc_vector[2],Total_conductance[2]);    

	} //end if const_N_omega>1

	printf("\n");


	sprintf(dirPathpos_processing_summary_FileName, "%s/pos_processing_summary.txt",results_folder); // path where the file is stored

	pos_processing_summary =fopen(dirPathpos_processing_summary_FileName,"w"); 
	fprintf(pos_processing_summary,"Material: %s\n Spectrum range(lambda) = %e--%e m \n Volume of each subvolume = %e \n",material, initial,final,delta_V_1);  // \nThermal power dissipated for thermal object 1= 
																				   //fopen(dirPathpos_processing_summary_FileName,"a"); // Opens file for appending. (File remains unchanged, file pointer gets moved to end) https://stackoverflow.com/questions/16427664/trying-not-to-overwrite-a-file-in-c/16427698
																				   //fprintf(pos_processing_summary,"Time counter: %f s\n",time_spent);
	fclose(pos_processing_summary);

	for (int iTcalc = 0; iTcalc < N_Tcalc; iTcalc++) 
	{
		pos_processing_summary =fopen(dirPathpos_processing_summary_FileName,"a"); 
		//fprintf(pos_processing_summary,"Temperature %eK: Total conductance = %e\n",T_calc, Total_conductance[0]); 
		fprintf(pos_processing_summary,"Temperature %eK: Total conductance = %e\n",Tcalc_vector[iTcalc], Total_conductance[iTcalc]); 
		fclose(pos_processing_summary);
	}

	free(Q_tot_thermal_object);
	free(Total_conductance); 


	clock_t end = clock(); //end timer
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC; //calculate time for the code
								    //printf("Time counter: %f s\n",time_spent);

	fopen(dirPathpos_processing_summary_FileName,"a"); // Opens file for appending. (File remains unchanged, file pointer gets moved to end) https://stackoverflow.com/questions/16427664/trying-not-to-overwrite-a-file-in-c/16427698
	fprintf(pos_processing_summary,"Time counter: %f s\n",time_spent);
	fclose(pos_processing_summary);


	// Mac shows the usage in bytes and Linux in kilobytes
	//   if (__linux__)
	//   {
	printf("Usage: %ld + %ld = %ld kb\n", baseline, get_mem_usage()-baseline,get_mem_usage());
	//   }
	//   else //MAC 
	//   {
	//       printf("Usage: %ld + %ld = %ld bytes\n", baseline, get_mem_usage()-baseline,get_mem_usage());
	//   }

	/*   // Condition used to compute memory_usage for one frequency.
	     if ( single_analysis =='y')
	     {
	     char memory_folder[100];
	     sprintf(memory_folder, "memory_usage"); // How to store words in an array in C? https://www.geeksforgeeks.org/how-to-store-words-in-an-array-in-c/

	     if (stat(memory_folder, &st) == -1) // This condition does not allow that new folders with the same name are created 
	     {
	     mkdir(memory_folder, 0700); // Create directory for the frequency that is analyzed
	     }

	     FILE * memory_usage; // 
	     char dirPathMem_FileName[260];
	     sprintf(dirPathMem_FileName, "%s/memory_usage_%s.csv",memory_folder,folder_subvol); // path where the file is stored
	     memory_usage =fopen(dirPathMem_FileName,"w"); 
	     fprintf(memory_usage,"%ld",get_mem_usage());      
	     fclose(memory_usage);
	     }
	     */ 

	free(x_omega);
	free(lambda);
	free(omega);
	free(delta_V_vector);
	free(T_vector);

	free(G_12_omega_SGF);

	free(trapz_Q);
	free(Q_subvol);
	free(Q_omega_subvol);
	free(Q_omega_thermal_object);
	free(results_folder);

} //end of code 

// #################################################################
// ##################### END OF THE CODE ###########################
// #################################################################
