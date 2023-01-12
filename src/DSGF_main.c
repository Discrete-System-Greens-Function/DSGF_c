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
#include "iterative_solver.h"
#include "array_functions.h"

#include "time.h"

#include "geometry/sphere.h"
#include "geometry/thin_film.h"
#include "material/SiO2.h"
#include "material/SiC.h"
#include "material/u_SiC.h"
#include "material/SiN.h"

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
	double (*lambda) = malloc(sizeof *lambda *const_N_omega); 
	double (*omega) = malloc(sizeof *omega *const_N_omega); 

	double (*delta_V_vector) = malloc(sizeof *delta_V_vector *tot_sub_vol);  //Vector of all subvolume size. Combines delta_V_1 and delta_V_2 in the same array
	double (*T_vector) = malloc(sizeof *T_vector *tot_sub_vol); // (N x 1) vector of all subvolume temperatures [K]
	double complex (*alpha_0) = malloc(sizeof *alpha_0 *tot_sub_vol); //Bare polarizability [m^3]

	double (*G_12_omega_SGF)[N_Tcalc] = calloc(const_N_omega, sizeof(*G_12_omega_SGF));

	double (*modulo_r_i_j)[tot_sub_vol] = malloc(sizeof *modulo_r_i_j * tot_sub_vol);
	double complex (*G_element)[tot_sub_vol] = malloc(sizeof *G_element * tot_sub_vol); 
	double (*sum_trans_coeff) = malloc(sizeof *sum_trans_coeff *const_N_omega); 

	double complex (*Q_omega_subvol) = malloc(sizeof * Q_omega_subvol * tot_sub_vol); // power dissipated per subvolume
	double (*Q_omega_thermal_object)[const_N_omega] = malloc(sizeof * Q_omega_thermal_object * const_N_bulk_objects);  
	double (*Q_subvol)[const_N_omega] = calloc(tot_sub_vol, sizeof(*Q_subvol));

	// ######### Properties for thermal objects ###########
	printf("Simulation for a total of %d dipoles in %d thermal objects\n",tot_sub_vol,const_N_bulk_objects);

	char dirPathFileNameDISCRETIZATION[260]; // https://stackoverflow.com/questions/46612504/creating-directories-and-files-using-a-loop-in-c 

	if(strcmp(geometry,"sphere")==0)
	{
		read_geometry_sphere(&d, &radius, &T1, &T2);
		radius1 = radius; // perfect same-sized spheres
		radius2 = radius; // perfect same-sized spheres
		vol1 = vol_sphere(radius1, pi); // calls function that calculates the volume for the sphere 1
		vol2 = vol_sphere(radius2, pi); // calls function that calculates the volume for the sphere 2
		delta_V_1 = vol1/const_N_subvolumes_per_object; // defines the subvolumes' volume for sphere 1
		delta_V_2 = vol2/const_N_subvolumes_per_object; // defines the subvolumes' volume for sphere 2

		subvol shape_file[const_N_subvolumes_per_object];

		sprintf(dirPathFileNameDISCRETIZATION, "discretizations/sphere/sphere_subvol_%d.txt",const_N_subvolumes_per_object); // path where the file is stored
		populate_subvol_struct(dirPathFileNameDISCRETIZATION, const_N_subvolumes_per_object, shape_file);	

		double origin1[3] = {radius1,radius1,radius1};
		double origin2[3]= {origin1[0]+radius1+d+radius2,origin1[1]+radius2-radius1,origin1[2]+radius2-radius1};

		populate_R_sphere(tot_sub_vol, const_N_subvolumes_per_object, origin1, origin2, delta_V_1, delta_V_2, shape_file, R);
	}   

	if(strcmp(geometry,"thin-films")==0)
	{
		read_geometry_thin_films(&d, &Lx, &Ly, &Lz, &T1, &T2);

		vol1 = Lx*Ly*Lz; // calculates the volume for membrane 1
		vol2 = vol1;     // defines the volume of membrane 2 = membrane 1
		delta_V_1 = vol1/const_N_subvolumes_per_object; // defines the subvolumes' volume for membrane 1
		delta_V_2 = vol2/const_N_subvolumes_per_object;  // defines the subvolumes' volume for membrane 1

		subvol shape_filetf[tot_sub_vol];
	
		int Lx_int, Ly_int, Lz_int, d_int;
		Lx_int = Lx*pow(10,9); 
		Ly_int = Ly*pow(10,9); 
		Lz_int = Lz*pow(10,9); 
		d_int = d*pow(10,9); 
		sprintf(dirPathFileNameDISCRETIZATION, "discretizations/thin-film/%d_thin_films_Lx%dnm_Ly%dnm_Lz%dnm_d%dnm_N%d_discretization.txt",const_N_bulk_objects,  Lx_int, Ly_int, Lz_int, d_int, tot_sub_vol);	
		populate_subvol_struct(dirPathFileNameDISCRETIZATION, tot_sub_vol, shape_filetf);

		populate_R_thin_film(tot_sub_vol, shape_filetf, R);
	}

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

	if(save_power_dissipated =='Y'){
		//EXPORT R
		for (int i_subvol=0; i_subvol<tot_sub_vol;i_subvol++) //tot_sub_vol
		{
			{
				FILE * vector_subvolumes_lattice; //append
				char dirPathVector_subvolumes_lattice_FileName[260];
				sprintf(dirPathVector_subvolumes_lattice_FileName, "%s/vector_subvolumes_lattice.csv",results_folder); // path where the file is stored
				if(i_subvol == 0) vector_subvolumes_lattice =fopen(dirPathVector_subvolumes_lattice_FileName,"w"); //write
				else vector_subvolumes_lattice = fopen(dirPathVector_subvolumes_lattice_FileName, "a"); //append
				fprintf(vector_subvolumes_lattice,"%e,%e,%e \n",R[i_subvol][0],R[i_subvol][1],R[i_subvol][2]); 
				fclose(vector_subvolumes_lattice);
			}
		} 
		//EXPORT delta_V_vector
		for (int i_subvol=0; i_subvol<tot_sub_vol; i_subvol++)
		{
			{
				FILE * vector_subvolumes_volume; //append
				char dirPathVector_subvolumes_volume_FileName[260];
				sprintf(dirPathVector_subvolumes_volume_FileName, "%s/vector_subvolumes_volume.csv",results_folder); // path where the file is stored
				if(i_subvol == 0) vector_subvolumes_volume =fopen(dirPathVector_subvolumes_volume_FileName,"w"); //write
				else vector_subvolumes_volume = fopen(dirPathVector_subvolumes_volume_FileName, "a"); //append
				fprintf(vector_subvolumes_volume,"%e\n",delta_V_vector[i_subvol]); 
				fclose(vector_subvolumes_volume);
			}
		}		
	} // end if save_power_dissipated


	double initial,final;

	if(strcmp(material,"SiO2")==0 || strcmp(material,"SiC")==0 || strcmp(material,"u-SiC")==0) 
	{
		//Uniform spectrum
		initial = 5.e-6;
		final = 25.e-6;
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

		//################## FREE-SPACE GREEN'S FUNCTION AND INTERACTION A MATRIX ##################### 
		// Fill terms for G^0: 

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
					r[i][j][i_alpha]= R[i][i_alpha] - R[j][i_alpha] ; //general case 
					abs_r_ij[i][j][i_alpha] = fabs(r[i][j][i_alpha]);
					modulo_r_i_j[i][j] += pow(abs_r_ij[i][j][i_alpha],2);
				}
				modulo_r_i_j[i][j] = sqrt(modulo_r_i_j[i][j]);
			}
		}

		for (int i_i = 0; i_i < tot_sub_vol; i_i++)
		{ 
			for (int i_j = 0; i_j < tot_sub_vol; i_j++)
			{ 
				for (int i_alpha = 0; i_alpha < 3; i_alpha++)
				{  
					unit_r_ij[i_i][i_j][i_alpha] = r[i_i][i_j][i_alpha]/modulo_r_i_j[i_i][i_j]; // ˆr -- unit distance 
					transpose[i_j][i_i][i_alpha] = r[i_i][i_j][i_alpha]/modulo_r_i_j[i_i][i_j]; // https://www.programiz.com/c-programming/examples/matrix-transpose        
					unit_conj_r_ij[i_i][i_j][i_alpha] = conj(transpose[i_j][i_i][i_alpha]);// r^+ Conjugate transpose unit distance 
				}
			}
		}

		for (int i_i = 0; i_i < tot_sub_vol; i_i++)
		{ 
			for (int i_j = 0; i_j < tot_sub_vol; i_j++)
			{ 
				for (int i_alpha = 0; i_alpha < 3; i_alpha++)
				{  
					for (int j_alpha = 0; j_alpha < 3; j_alpha++)
					{ 
						r_i_j_outer_r_i_j[i_i][i_j][i_alpha][j_alpha] = unit_r_ij[i_i][i_j][i_alpha]*unit_conj_r_ij[i_i][i_j][j_alpha];

					}
				}    
			}
		}

		// ################### MATRICES STRUCTURE LOOPS ###########################
		// 3N X 3N Matrices structure loops for G^0 and A:

		//G^0_ij when i=j:
		double (*a_j) = malloc(sizeof *a_j *tot_sub_vol); 
		double (*part1ii) = malloc(sizeof *part1ii *tot_sub_vol); 
		double complex (*part2ii) = malloc(sizeof *part2ii *tot_sub_vol);
		double complex (*part2iiexp) = malloc(sizeof *part2iiexp *tot_sub_vol); 
		double complex (*part3ii) = malloc(sizeof *part3ii *tot_sub_vol); 

		//G^0_ij when i!=j:
		double complex (*part1ij) = malloc(sizeof *part1ij *tot_sub_vol); 
		double complex (*part1aij) = malloc(sizeof *part1aij *tot_sub_vol);
		double complex (*part1aijexp) = malloc(sizeof *part1aijexp *tot_sub_vol);
		double complex (*part2ij) = malloc(sizeof *part2ij *tot_sub_vol); 
		double complex (*part3ij) = malloc(sizeof *part3ij *tot_sub_vol); 

		double (*eyeG_0)[3] = calloc(3, sizeof(*eyeG_0)); 
		double (*eyeA)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*eyeA)); 

		// Linear system AG=G^0 
		double complex (*G_0)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*G_0)); 
		double complex (*A)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*A)); 

		// eq. 25 from Lindsay's paper 
		for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
		{
			for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
			{
				if (ig_0!=jg_0)
				{
					part1aij[ig_0] = 0.+ k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0]*I; // com i term 
					part1aijexp[ig_0]= cexp(k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0]*I);
					part1ij[ig_0] = part1aijexp[ig_0]/(4.*pi*modulo_r_i_j[ig_0][jg_0]);
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

							}
							else
							{
								eyeG_0[i_subG_0][j_subG_0] = 0.;     // 3x3 Identity matrix for G^0:
								eyeA[ig_0][jg_0][i_subG_0][j_subG_0] = 0.; // 3Nx3N identity matrix for A: 
							}
							G_0[ig_0][jg_0][i_subG_0][j_subG_0] = part1ij[ig_0]*((part2ij[ig_0]*eyeG_0[i_subG_0][j_subG_0])-(part3ij[ig_0]*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
							A[ig_0][jg_0][i_subG_0][j_subG_0] = eyeA[ig_0][jg_0][i_subG_0][j_subG_0] - pow(k_0,2)*alpha_0[ig_0]*G_0[ig_0][jg_0][i_subG_0][j_subG_0]; 
						}    
					}
				}    
			}    
		} //end ig_0 

		// eq. 26 from Lindsay's paper: 


		for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
		{
			a_j[ig_0] = a_j_function(delta_V_vector[ig_0], pi);
			part1ii[ig_0] = 1./(3.*delta_V_vector[ig_0]*epsilon_ref*pow(k_0,2));
			part2ii[ig_0] = a_j[ig_0]*k_0*sqrt(epsilon_ref)*I; // com i term 
			part2iiexp[ig_0] = cexp(0. + a_j[ig_0]*k_0*sqrt(epsilon_ref)*I); 
			// part3ii is inside brackets
			part3ii[ig_0] = part2iiexp[ig_0]*(1-part2ii[ig_0]) - 1. ;
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
							}
							else
							{
								eyeG_0[i_subG_0][j_subG_0] = 0.;     // 3x3 Identity matrix for G^0:
								eyeA[ig_0][jg_0][i_subG_0][j_subG_0] = 0.; // 3Nx3N identity matrix for A:
							}
							G_0[ig_0][jg_0][i_subG_0][j_subG_0] = eyeG_0[i_subG_0][j_subG_0]*part1ii[ig_0]*(2.*part3ii[ig_0]-1.); 
							A[ig_0][jg_0][i_subG_0][j_subG_0] = eyeA[ig_0][jg_0][i_subG_0][j_subG_0] - pow(k_0,2)*alpha_0[ig_0]*G_0[ig_0][jg_0][i_subG_0][j_subG_0]; 
						}
					} 
				}   //end jg_0   
			} //end i_subG_0     
		} //end ig_0 


		free(eyeG_0);
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

		//printf("##################### \n SOLVE LINEAR SYSTEM AG=G^0 \n##################### \n");
		//printf("##################### \n LAPACK/LAPACKE ZGELS ROUTINE \n##################### \n");
		//Description of ZGELS: https://extras.csc.fi/math/nag/mark21/pdf/F08/f08anf.pdf

		double complex (*G_sys)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys));

		if(solution =='D')
		{
			printf("Direct inversion status: ");
			double complex (*A_direct) = malloc(sizeof *A_direct *lda*n); 
			double complex (*b_direct) = malloc(sizeof *b_direct *ldb*nrhs);

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
							gpack=gpack + 1;
						}    
					}
				}        
			} 

			free(b_direct); 
			printf("concluded\n");
		}
		free(A);


		if(solution =='I')
		{ 
			printf("Iterative status:\n m= ");

			double complex (*G_sys_old)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys_old)); 
			double complex (*G_sys_new)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys_old)); 
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
			matrix_reshape(3, tot_sub_vol, G_sys_old, G_0);


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
				memset(A_iterative, 0, sizeof *A_iterative *3*3);
				memset(b_iterative, 0, sizeof *b_iterative *3*3);

				// Next, solve all systems of equations for ii not equal to mm		
				for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
				{ 
					for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
					{
						ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
						for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //lower triangular matrix
						{
							if(ig_0 != mm)
							{
								for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
								{
									jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
									G_sys_prod = 0.;
									for(int m_sub = 0;  m_sub < 3;  m_sub++)//loop for matricial multiplication
									{
										mm_2d = (3*mm + m_sub);
										G_sys_prod += G_sys_old[ig_0_2d][mm_2d]*G_sys_new[mm_2d][jg_0_2d];
									}

									G_sys_new[ig_0_2d][jg_0_2d] = G_sys_old[ig_0_2d][jg_0_2d] + pow(k,2)*alpha_0[mm]*G_sys_prod;

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


		// #################################################################
		//  Spectral and total conductance between bulk objects at temp. T
		// ###################  Thermal power dissipated ###################
		// #################################################################

		double complex (*trans_coeff)[tot_sub_vol] = calloc(tot_sub_vol, sizeof(*trans_coeff));

		double complex (*G_sys_cross)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*G_sys_cross)); 
		double complex (*transpose_G_sys)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*transpose_G_sys)); 
		double(*inner_sum) = malloc(sizeof * inner_sum * tot_sub_vol);

		int counter = 0;
		sum_trans_coeff[i_omega] = 0.;    

		for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
		{
			for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
			{
				for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
				{
					ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
					{
						jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
						// Extract one Green's function component: G_sys[ig_0][jg_0][i_subG_0]
						transpose_G_sys[jg_0][ig_0][i_subG_0][j_subG_0] = G_sys[ig_0_2d][jg_0_2d];
						G_sys_cross[ig_0][jg_0][i_subG_0][j_subG_0]= conj(transpose_G_sys[jg_0][ig_0][i_subG_0][j_subG_0]);
						G_element[ig_0][jg_0]+=G_sys[ig_0_2d][jg_0_2d]*G_sys_cross[ig_0][jg_0][i_subG_0][j_subG_0]; 

					}    
				}
				// Transmissivity coefficient matrix tau(omega) for comparison with Czapla Mathematica output [dimensionless]
				trans_coeff[ig_0][jg_0] = 4.*pow(k_0,4)*delta_V_vector[ig_0]*delta_V_vector[jg_0]*cimag(epsilon)*cimag(epsilon)*G_element[ig_0][jg_0] ;  
				// Thermal power dissipated calculation, based on the matlab code (Using Tervo's Eq. 26)
				if (ig_0 != jg_0) 
				{
					inner_sum[jg_0] = (theta_function(omega_value, T_vector[jg_0], h_bar, k_b) - theta_function(omega_value, T_vector[ig_0], h_bar, k_b)) * trans_coeff[ig_0][jg_0];
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

		double (*dtheta_dT) = malloc(sizeof *dtheta_dT *N_Tcalc);// function used to calculate conductance, modified for several temperatures

		for (int iTcalc = 0; iTcalc < N_Tcalc; iTcalc++) // EDIT VALUE :: change 1 to N_Tcalc for the temperature loop
		{
			dtheta_dT[iTcalc] = dtheta_dT_function(omega_value,Tcalc_vector[iTcalc], h_bar, k_b); 
			G_12_omega_SGF[i_omega][iTcalc] = dtheta_dT[iTcalc]*sum_trans_coeff[i_omega];
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
		free(dtheta_dT);


		// #################################################################
		// #############  Clear values for next frequency ##################
		// #################################################################

		memset(modulo_r_i_j, 0, sizeof *modulo_r_i_j * tot_sub_vol); //modulo_r_i_j
		memset(G_element, 0, sizeof *G_element * tot_sub_vol); //G_element
		memset(sum_trans_coeff, 0, sizeof sum_trans_coeff);
		memset(Q_omega_subvol, 0, sizeof* Q_omega_subvol* tot_sub_vol); 

	} // END OMEGA VALUE LOOP FOR FREQUENCY RANGE ANALYSIS, loop started in line 481

	free(R);
	free(sum_trans_coeff);
	free(modulo_r_i_j);
	free(G_element);

	free(alpha_0);

	double (*Total_conductance) = malloc(sizeof *Total_conductance *N_Tcalc); 
	double(*Q_tot_thermal_object) = malloc(sizeof * Q_tot_thermal_object * const_N_bulk_objects); 
	double (*trapz_Q) = malloc(sizeof *trapz_Q *tot_sub_vol); // Definition for trapezoidal integration. Used in total power dissipated

	trapz=0.;
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


	free(Q_tot_thermal_object);
	free(Total_conductance); 

	printf("Usage: %ld + %ld = %ld kb\n", baseline, get_mem_usage()-baseline,get_mem_usage());

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
