#include "computational/GreensFunction.h"
#include <math.h>
#include "functions_DSGF.h"
#include <stdlib.h>
#include "computational/solvers/iterative_solver.h"
#include "mkl.h"

void setup_G_0_matrices(int tot_sub_vol, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double R[][3]){
	
	double (*r)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*r));
	if (r == NULL){
		printf("Failure with memory=%ld before spectral analysis during setup_G_0_matrices. ",get_mem_usage());
		//return 1;
		exit(1);
	}
	/*
	double (*abs_r_ij)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*abs_r_ij));
	if (abs_r_ij == NULL){
		printf("Failure with memory=%ld before spectral analysis during setup_G_0_matrices. ",get_mem_usage());
		//return 1;
		exit(1);
	}
	*/

	for (int i=0; i<tot_sub_vol; i++)
	{   
		for (int j=0; j<tot_sub_vol; j++)
		{
			for (int i_alpha=0; i_alpha<3; i_alpha++)
			{
				r[i][j][i_alpha]= R[i][i_alpha] - R[j][i_alpha] ; //general case; can be a scalar 
				//abs_r_ij[i][j][i_alpha] = fabs(r[i][j][i_alpha]); // can be a scalar
				double abs_r_ij = fabs(r[i][j][i_alpha]);
				//modulo_r_i_j[i][j] += pow(abs_r_ij[i][j][i_alpha],2); // change the name in this step, can be a scalar
				modulo_r_i_j[i][j] += pow(abs_r_ij,2); // change the name in this step, can be a scalar 
			}
			modulo_r_i_j[i][j] = sqrt(modulo_r_i_j[i][j]);
		}
	}
	//free(abs_r_ij);

	
	double (*unit_r_ij)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*unit_r_ij));
	if (unit_r_ij == NULL){
		printf("Failure with memory=%ld before spectral analysis during setup_G_0_matrices. ",get_mem_usage());
		//return 1;
		exit(1);
	}
	double complex (*unit_conj_r_ij)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*unit_conj_r_ij));
	if (unit_conj_r_ij == NULL){
		printf("Failure with memory=%ld before spectral analysis during setup_G_0_matrices. ",get_mem_usage());
		//return 1;
		exit(1);
	}
	

	for (int i_i = 0; i_i < tot_sub_vol; i_i++)
	{ 
		for (int i_j = 0; i_j < tot_sub_vol; i_j++)
		{ 
			for (int i_alpha = 0; i_alpha < 3; i_alpha++)
			{  
				//double unit_r_ij = r[i_i][i_j][i_alpha]/modulo_r_i_j[i_i][i_j];
				unit_r_ij[i_i][i_j][i_alpha] = r[i_i][i_j][i_alpha]/modulo_r_i_j[i_i][i_j]; // ˆr -- unit distance 
				double transpose = r[i_i][i_j][i_alpha]/modulo_r_i_j[i_i][i_j]; // https://www.programiz.com/c-programming/examples/matrix-transpose        
				unit_conj_r_ij[i_i][i_j][i_alpha] = conj(transpose);// r^+ Conjugate transpose unit distance 
				//double complex unit_conj_r_ij = conj(transpose); 
			}
		}
	}
	free(r);


	//#pragma omp parallel for collapse(3)
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
	free(unit_r_ij);
	free(unit_conj_r_ij);
}

void get_G0_A_matrices(int tot_sub_vol, double complex G_0[tot_sub_vol][tot_sub_vol][3][3], double complex A[tot_sub_vol][tot_sub_vol][3][3], double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double complex alpha_0[tot_sub_vol], double delta_V_vector[tot_sub_vol],char wave_type){

	// ################### MATRICES STRUCTURE LOOPS ###########################
	// 3N X 3N Matrices structure loops for G^0 and A:
	double denom_1, denom_2 ; // used in G^0_ij function
	double complex const_1, const_2, const_3;

	double eyeG_0[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};


	// eq. 25 from Lindsay's paper 

	//for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++)//
	for (int jg_0 = 0; jg_0 < tot_sub_vol-1; jg_0++) //tot_sub_vol
	{
		//for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++)//
		for (int ig_0 = jg_0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
		{
			//if (ig_0!=jg_0)
			//{
			const_1 = cexp(k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0]*I)/(4.*pi*modulo_r_i_j[ig_0][jg_0]); 
			denom_1 = epsilon_ref*pow(k_0*modulo_r_i_j[ig_0][jg_0],2);
			denom_2 = k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0];
			const_2 = (1. - 1./denom_1 + 1.*I/denom_2 ) ;
			const_3 = (1. - 3./denom_1 + 3.*I/denom_2) ;
			
			//split of G_0:
			/*
			if (jg_0<tot_sub_vol/2 && ig_0>=tot_sub_vol/2) //subvolumes in different objects
			{
				const_2 = (1. - 1./denom_1 + 1.*I/denom_2 ) ; // total 
				const_3 = (1. - 3./denom_1 + 3.*I/denom_2) ;  // total 
				
				if(wave_type == 'T')  
				{ 
					const_2 = (1. - 1./denom_1 + 1.*I/denom_2 ) ; // total 
					const_3 = (1. - 3./denom_1 + 3.*I/denom_2) ;  // total 
				}	
				//Goal: compute only the propagating wave contribution of DSGF
				if(wave_type == 'P') 
				{
					const_2 = (1. ) ; //propagating only
					const_3 = (1. ) ; //propagating only
				}	
				//Goal: compute only the evanescent wave contribution of DSGF
				if(wave_type == 'E') 
				{
					const_2 = (- 1./denom_1 + 1.*I/denom_2 ) ; //evanescent only
					const_3 = (- 3./denom_1 + 3.*I/denom_2) ; //evanescent only
				}
					
			}
			else 
			{
				const_2 = (1. - 1./denom_1 + 1.*I/denom_2 ) ;
				const_3 = (1. - 3./denom_1 + 3.*I/denom_2) ; 
			}
			*/
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					G_0[ig_0][jg_0][i_subG_0][j_subG_0] = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
					G_0[jg_0][ig_0][i_subG_0][j_subG_0] = G_0[ig_0][jg_0][i_subG_0][j_subG_0];
					A[ig_0][jg_0][i_subG_0][j_subG_0] = 0 - pow(k_0,2)*alpha_0[ig_0]*G_0[ig_0][jg_0][i_subG_0][j_subG_0]; 
					A[jg_0][ig_0][i_subG_0][j_subG_0] =A[ig_0][jg_0][i_subG_0][j_subG_0];
				}    
			}
			//}
		}    
	} //end j_subG_0


	// eq. 26 from Lindsay's paper: 

	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		double a_j = a_j_function(delta_V_vector[ig_0], pi);
		double part1ii = 1./(3.*delta_V_vector[ig_0]*epsilon_ref*pow(k_0,2));
		double complex part2ii = a_j*k_0*sqrt(epsilon_ref)*I; // com i term 
		double complex part2iiexp = cexp(0. + a_j*k_0*sqrt(epsilon_ref)*I); 
		// part3ii is inside brackets
		double complex part3ii = part2iiexp*(1-part2ii) - 1. ;
		for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
		{
			for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
			{ 
				if (ig_0==jg_0) // if i=j:
				{
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
					{
						G_0[ig_0][jg_0][i_subG_0][j_subG_0] = eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.); 
						A[ig_0][jg_0][i_subG_0][j_subG_0] = eyeG_0[i_subG_0][j_subG_0] - pow(k_0,2)*alpha_0[ig_0]*G_0[ig_0][jg_0][i_subG_0][j_subG_0]; 
					}
				} 
			}   //end jg_0   
		} //end i_subG_0     
	} //end ig_0 
}

void get_G0_matrix(int tot_sub_vol, double complex G_0[tot_sub_vol][tot_sub_vol][3][3], double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char wave_type){ //double complex alpha_0[tot_sub_vol],

	// ################### MATRICES STRUCTURE LOOPS ###########################
	// 3N X 3N Matrices structure loops for G^0 and A:
	double denom_NF, denom_IF ; // used in G^0_ij function
	double complex const_1, const_2, const_3;

	double eyeG_0[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};


	// eq. 25 from Lindsay's paper 
	//#pragma omp parallel for collapse(2)
	for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++)//
	//for (int jg_0 = 0; jg_0 < tot_sub_vol-1; jg_0++) //tot_sub_vol
	{
		for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++)//
		//for (int ig_0 = jg_0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
		{
			//if (ig_0!=jg_0)
			//{
			const_1 = cexp(k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0]*I)/(4.*pi*modulo_r_i_j[ig_0][jg_0]); 
			denom_NF = epsilon_ref*pow(k_0*modulo_r_i_j[ig_0][jg_0],2);
			denom_IF = k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0];
			const_2 = (1. - 1./denom_NF + 1.*I/denom_IF ) ;
			const_3 = (1. - 3./denom_NF + 3.*I/denom_IF) ;
			
			//split of G_0:
			/*
			if (jg_0<tot_sub_vol/2 && ig_0>=tot_sub_vol/2) //subvolumes in different objects
			{
				if(wave_type == 'T')  
				{ 
					const_2 = (1. - 1./denom_NF + 1.*I/denom_IF ) ; // total 
					const_3 = (1. - 3./denom_NF + 3.*I/denom_IF) ;  // total 
				}	
				//Goal: compute only the propagating wave contribution of DSGF
				if(wave_type == 'P') 
				{
					const_2 = (1. ) ; //propagating only
					const_3 = (1. ) ; //propagating only
				}	
				//Goal: compute only the evanescent wave contribution of DSGF
				if(wave_type == 'E') 
				{
					//const_2 = (- 1./denom_NF + 1.*I/denom_IF ) ; //evanescent only
					//const_3 = (- 3./denom_NF + 3.*I/denom_IF) ; //evanescent only
					const_2 = ( - 1./denom_NF ) ; //test NF term
					const_3 = ( - 3./denom_NF ) ; //test NF term
				}
					
			}
			else 
			{
				const_2 = (1. - 1./denom_NF + 1.*I/denom_IF ) ;
				const_3 = (1. - 3./denom_NF + 3.*I/denom_IF) ; 
			}
			*/
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					G_0[ig_0][jg_0][i_subG_0][j_subG_0] = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
					G_0[jg_0][ig_0][i_subG_0][j_subG_0] = G_0[ig_0][jg_0][i_subG_0][j_subG_0];
					//A[ig_0][jg_0][i_subG_0][j_subG_0] = 0 - pow(k_0,2)*alpha_0[ig_0]*G_0[ig_0][jg_0][i_subG_0][j_subG_0]; 
					//A[jg_0][ig_0][i_subG_0][j_subG_0] =A[ig_0][jg_0][i_subG_0][j_subG_0];
				}    
			}
			//}
		}    
	} //end j_subG_0


	// eq. 26 from Lindsay's paper: 
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		double a_j = a_j_function(delta_V_vector[ig_0], pi);
		double part1ii = 1./(3.*delta_V_vector[ig_0]*epsilon_ref*pow(k_0,2));
		double complex part2ii = a_j*k_0*sqrt(epsilon_ref)*I; // com i term 
		double complex part2iiexp = cexp(0. + a_j*k_0*sqrt(epsilon_ref)*I); 
		// part3ii is inside brackets
		double complex part3ii = part2iiexp*(1-part2ii) - 1. ;
		for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
		{
			for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //tot_sub_vol
			{ 
				if (ig_0==jg_0) // if i=j:
				{
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
					{
						G_0[ig_0][jg_0][i_subG_0][j_subG_0] = eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.); 
						//A[ig_0][jg_0][i_subG_0][j_subG_0] = eyeG_0[i_subG_0][j_subG_0] - pow(k_0,2)*alpha_0[ig_0]*G_0[ig_0][jg_0][i_subG_0][j_subG_0]; 
					}
				} 
			}   //end jg_0   
		} //end i_subG_0     
	} //end ig_0 
}


void get_G0_matrix_memory(int tot_sub_vol, double complex G_0[tot_sub_vol][tot_sub_vol][3][3], double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char wave_type){

	// ################### MATRICES STRUCTURE LOOPS ###########################
	// 3N X 3N Matrices structure loops for G^0 and A:
	double denom_NF, denom_IF ; // used in G^0_ij function
	double complex const_1, const_2, const_3;

	double eyeG_0[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};

	double a_j, part1ii;
	double complex part2ii, part2iiexp,part3ii;
	
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	//for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++)//
	//for (int jg_0 = 0; jg_0 < tot_sub_vol-1; jg_0++) //tot_sub_vol
	{
		for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++)//
		//for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++)//
		//for (int ig_0 = jg_0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
		{
			if (ig_0!=jg_0) // eq. 25 from Walter et al., PRB 2021
			{
			const_1 = cexp(k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0]*I)/(4.*pi*modulo_r_i_j[ig_0][jg_0]); 
			denom_NF = epsilon_ref*pow(k_0*modulo_r_i_j[ig_0][jg_0],2);
			denom_IF = k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0];
			const_2 = (1. - 1./denom_NF + 1.*I/denom_IF ) ;
			const_3 = (1. - 3./denom_NF + 3.*I/denom_IF) ;
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					G_0[ig_0][jg_0][i_subG_0][j_subG_0] = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
					G_0[jg_0][ig_0][i_subG_0][j_subG_0] = G_0[ig_0][jg_0][i_subG_0][j_subG_0];
				}    
			}
			}//end if (ig_0!=jg_0)
			else //if (ig_0==jg_0) // eq. 26 from Walter et al., PRB 2021 
			{
				a_j = a_j_function(delta_V_vector[ig_0], pi);
				part1ii = 1./(3.*delta_V_vector[ig_0]*epsilon_ref*pow(k_0,2));
				part2ii = a_j*k_0*sqrt(epsilon_ref)*I; // com i term 
				part2iiexp = cexp(0. + a_j*k_0*sqrt(epsilon_ref)*I); 
				part3ii = part2iiexp*(1-part2ii) - 1. ; // part3ii is inside brackets
				for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
				{
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
					{
						G_0[ig_0][jg_0][i_subG_0][j_subG_0] = eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.); 
					}
				} //end i_subG_0
			}//end if (ig_0=jg_0)	
		}    
	} //end j_subG_0

}

void get_A_matrix(int tot_sub_vol, double complex G_0[tot_sub_vol][tot_sub_vol][3][3], double complex A[tot_sub_vol][tot_sub_vol][3][3], double k_0, double complex alpha_0[tot_sub_vol]){ //,char wave_type
	
	double complex (*alpha_0_matrix)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*alpha_0_matrix));
	
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++)//
	{
		for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
		{
			for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++)//
			{
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					if (ig_0==jg_0 && i_subG_0==j_subG_0)
					{
						alpha_0_matrix[ig_0][jg_0][i_subG_0][j_subG_0] = alpha_0[ig_0];
					}
					else
					{
						alpha_0_matrix[ig_0][jg_0][i_subG_0][j_subG_0] = 0.;
					}
				}
			}
		}
	}
	
	double complex (*alpha_0_matrix_2D)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*alpha_0_matrix_2D));
	matrix_reshape(3, tot_sub_vol, alpha_0_matrix_2D, alpha_0_matrix); // this reshapes 2 4D matrices to 2 2D matrices, where G0 and eyeA are 4D and G_sys_old and eyeA_2d are the respective 2D matrices
	free(alpha_0_matrix);

	double complex (*G_0_matrix_2D)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_0_matrix_2D));
	matrix_reshape(3, tot_sub_vol, G_0_matrix_2D, G_0); // this reshapes 2 4D matrices to 2 2D matrices, where G0 and eyeA are 4D and G_sys_old and eyeA_2d are the respective 2D matrices
	
	double complex (*prod)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*prod));

	double complex alpha_parameter = pow(k_0,2) + 0.0 * I;  // Scaling factor for A*B
	double complex beta_parameter = 0.0 + 0.0 * I;   // Scaling factor for C
	
	//cblas_zgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans, 3*tot_sub_vol, 3*tot_sub_vol, 3*tot_sub_vol,&alpha_parameter, G_0_matrix_2D,  3*tot_sub_vol, alpha_0_matrix_2D, 3*tot_sub_vol,&beta_parameter,prod, 3*tot_sub_vol);
	cblas_zgemm(
		CblasRowMajor,CblasNoTrans, CblasNoTrans, 
		3*tot_sub_vol, 3*tot_sub_vol, 3*tot_sub_vol,  	// Dimensions of matrices
		&alpha_parameter, 								// Scaling factor for A*B
		G_0_matrix_2D,  3*tot_sub_vol, 					// Matrix A
		alpha_0_matrix_2D, 3*tot_sub_vol, 				// Matrix B
		&beta_parameter, 								// Scaling factor for C
		prod, 3*tot_sub_vol								 // Matrix C
	);
	
	free(alpha_0_matrix_2D);
	free(G_0_matrix_2D);
		
	/*
	// slow matrix multiplication 
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{ 
		for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
		{
			for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //lower triangular matrix
			{
				for(int j_subG_0 = 0;  j_subG_0 < 3;  j_subG_0++)//loop for matricial multiplication
				{
							double complex G_0_alpha_matrix_prod = 0.;
							for (int object = 0; object < tot_sub_vol; object++) //slow part of the code!!!!!
							{
								for (int axis = 0; axis < 3; axis++)
								{
									G_0_alpha_matrix_prod += pow(k_0,2)*G_0[ig_0][object][i_subG_0][axis]*alpha_0_matrix[object][jg_0][axis][j_subG_0];
									//G_0_alpha_matrix_prod += pow(k_0,2)*G_0[ig_0][object][i_subG_0][axis]*alpha_0_matrix_transpose[jg_0][object][j_subG_0][axis];
								}
							}
							prod[ig_0][jg_0][i_subG_0][j_subG_0] = G_0_alpha_matrix_prod;	
				}// j_subG_0  
			} // jg_0   	
		}// i_subG_0   	                	
	} // ig_0    
	*/

	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++)//
	{
		for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
		{
			int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
			for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++)//
			{
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					if (ig_0==jg_0 && i_subG_0==j_subG_0) // if i=j:
					{
						//A[ig_0][jg_0][i_subG_0][j_subG_0] = 1. - pow(k_0,2)*alpha_0_matrix[ig_0][jg_0][i_subG_0][j_subG_0]*G_0[ig_0][jg_0][i_subG_0][j_subG_0]; // eq. 26 from Lindsay's paper
						//A[ig_0][jg_0][i_subG_0][j_subG_0] = 1. -  prod[ig_0][jg_0][i_subG_0][j_subG_0]; // eq. 26 from Lindsay's paper
						A[ig_0][jg_0][i_subG_0][j_subG_0] = 1. -  prod[ig_0_2d][jg_0_2d]; // eq. 26 from Lindsay's paper
					
					}	
					else
					{
						//A[ig_0][jg_0][i_subG_0][j_subG_0] = 0. - pow(k_0,2)*alpha_0_matrix[ig_0][jg_0][i_subG_0][j_subG_0]*G_0[ig_0][jg_0][i_subG_0][j_subG_0]; // eq. 26 from Lindsay's paper
						//A[ig_0][jg_0][i_subG_0][j_subG_0] = 0. -  prod[ig_0][jg_0][i_subG_0][j_subG_0]; // eq. 26 from Lindsay's paper
						A[ig_0][jg_0][i_subG_0][j_subG_0] = 0. -  prod[ig_0_2d][jg_0_2d]; // eq. 26 from Lindsay's paper
					}
				} //end j_subG_0 
			} //end jg_0  
		}  //end i_subG_0     
	} //end ig_0    

free(prod);

}


void get_G_old_matrix_memory(int tot_sub_vol, double complex G_old[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char wave_type){

	// ################### MATRICES STRUCTURE LOOPS ###########################
	// 3N X 3N Matrices structure loops for G^0 and A:
	double denom_NF, denom_IF ; // used in G^0_ij function
	double complex const_1, const_2, const_3;

	double eyeG_0[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};

	double a_j, part1ii;
	double complex part2ii, part2iiexp,part3ii;
	
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++)//
		{
			if (ig_0!=jg_0) // eq. 25 from Walter et al., PRB 2021
			{
			const_1 = cexp(k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0]*I)/(4.*pi*modulo_r_i_j[ig_0][jg_0]); 
			denom_NF = epsilon_ref*pow(k_0*modulo_r_i_j[ig_0][jg_0],2);
			denom_IF = k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0];
			const_2 = (1. - 1./denom_NF + 1.*I/denom_IF ) ;
			const_3 = (1. - 3./denom_NF + 3.*I/denom_IF) ;
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					G_old[ig_0_2d][jg_0_2d] = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
					G_old[jg_0_2d][ig_0_2d] = G_old[ig_0_2d][jg_0_2d];
				}    
			}
			}//end if (ig_0!=jg_0)
			else //if (ig_0==jg_0) // eq. 26 from Walter et al., PRB 2021 
			{
				a_j = a_j_function(delta_V_vector[ig_0], pi);
				part1ii = 1./(3.*delta_V_vector[ig_0]*epsilon_ref*pow(k_0,2));
				part2ii = a_j*k_0*sqrt(epsilon_ref)*I; // com i term 
				part2iiexp = cexp(0. + a_j*k_0*sqrt(epsilon_ref)*I); 
				part3ii = part2iiexp*(1-part2ii) - 1. ; // part3ii is inside brackets
				for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
				{
					int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
					{
						int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
						G_old[ig_0_2d][jg_0_2d] = eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.); 
					}
				} //end i_subG_0
			}//end if (ig_0=jg_0)	
		}    
	} //end j_subG_0

}


void get_G0_triangular_matrix(int tot_sub_vol, int size, double complex upperTriangularMatrix[size], double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char wave_type){
	//printf("G_0: ");
	// ################### MATRICES STRUCTURE LOOPS ###########################
	// 3N X 3N Matrices structure loops for G^0 and A:
	double denom_NF, denom_IF ; // used in G^0_ij function
	double complex const_1, const_2, const_3;

	double eyeG_0[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};

	double a_j, part1ii;
	double complex part2ii, part2iiexp,part3ii;
	int index = 0;
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
		{
			for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++)//
			{
				if (ig_0!=jg_0) // eq. 25 from Walter et al., PRB 2021
				{
					const_1 = cexp(k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0]*I)/(4.*pi*modulo_r_i_j[ig_0][jg_0]); 
					denom_NF = epsilon_ref*pow(k_0*modulo_r_i_j[ig_0][jg_0],2);
					denom_IF = k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0];
					const_2 = (1. - 1./denom_NF + 1.*I/denom_IF ) ;
					const_3 = (1. - 3./denom_NF + 3.*I/denom_IF) ;

					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
					{
						//printf("%d -",index );
						//G_0[ig_0][jg_0][i_subG_0][j_subG_0] = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
						//G_0[jg_0][ig_0][i_subG_0][j_subG_0] = G_0[ig_0][jg_0][i_subG_0][j_subG_0];
						upperTriangularMatrix[index] = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
						//printf(" G_0[%d]= (%e + i%e) ",index, creal(upperTriangularMatrix[index]),cimag(upperTriangularMatrix[index]));
						index++;
					}    
				}
				else //if (ig_0==jg_0) // eq. 26 from Walter et al., PRB 2021 
				{
					a_j = a_j_function(delta_V_vector[ig_0], pi);
					part1ii = 1./(3.*delta_V_vector[ig_0]*epsilon_ref*pow(k_0,2));
					part2ii = a_j*k_0*sqrt(epsilon_ref)*I; // com i term 
					part2iiexp = cexp(0. + a_j*k_0*sqrt(epsilon_ref)*I); 
					part3ii = part2iiexp*(1-part2ii) - 1. ; // part3ii is inside brackets
					//for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
					for(int j_subG_0 = i_subG_0; j_subG_0 < 3; j_subG_0++) // modified loop
					{
						//printf("%d -",index );
						upperTriangularMatrix[index] = eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.); 
						//printf(" G_0[%d]= (%e + i%e) ",index, creal(upperTriangularMatrix[index]),cimag(upperTriangularMatrix[index]));
						index++;
					}
				}	
			}//end jg_0
			//printf("\n"); 
		}   //end i_subG_0
	} //end ig_0

/*
	index =0;
	printf("\n");
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
		for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++)//
		{
				for(int j_subG_0 = i_subG_0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					printf(" G_0[%d]= (%e + i%e) ",index, creal(upperTriangularMatrix[index]),cimag(upperTriangularMatrix[index]));
					index++;
				}
			}
			printf("\n");
		}
	}
*/

}




void get_G_old_struct_matrix_memory(int tot_sub_vol, double complex G_old[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char wave_type){

	// ################### MATRICES STRUCTURE LOOPS ###########################
	// 3N X 3N Matrices structure loops for G^0 and A:
	//double denom_NF, denom_IF ; // used in G^0_ij function
	//double complex const_1, const_2, const_3;
	//double a_j, part1ii;
	// double complex part2ii, part2iiexp,part3ii;
	double eyeG_0[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};
	// #pragma omp parallel for collapse(2)
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++)//
		{
			if (ig_0!=jg_0) // eq. 25 from Walter et al., PRB 2021
			{
			double complex const_1 = cexp(k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0]*I)/(4.*pi*modulo_r_i_j[ig_0][jg_0]); 
			double denom_NF = epsilon_ref*pow(k_0*modulo_r_i_j[ig_0][jg_0],2);
			double denom_IF = k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0];
			double complex const_2 = (1. - 1./denom_NF + 1.*I/denom_IF ) ;
			double complex const_3 = (1. - 3./denom_NF + 3.*I/denom_IF) ;
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					//G_old[ig_0][jg_0][i_subG_0][j_subG_0]=const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
					//G_old[jg_0][ig_0][j_subG_0][i_subG_0]=G_old[ig_0][jg_0][i_subG_0][j_subG_0];
					G_old[ig_0_2d][jg_0_2d] = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
					G_old[jg_0_2d][ig_0_2d] = G_old[ig_0_2d][jg_0_2d];
				}    
			}
			}//end if (ig_0!=jg_0)
			else //if (ig_0==jg_0) // eq. 26 from Walter et al., PRB 2021 
			{
				double a_j = a_j_function(delta_V_vector[ig_0], pi);
				double part1ii = 1./(3.*delta_V_vector[ig_0]*epsilon_ref*pow(k_0,2));
				double complex part2ii = a_j*k_0*sqrt(epsilon_ref)*I; // com i term 
				double complex part2iiexp = cexp(0. + a_j*k_0*sqrt(epsilon_ref)*I); 
				double complex part3ii = part2iiexp*(1-part2ii) - 1. ; // part3ii is inside brackets
				for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
				{
					int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
					{
						int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
						//G_old[ig_0][jg_0][i_subG_0][j_subG_0]=eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.);
						G_old[ig_0_2d][jg_0_2d] = eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.); 
					}
				} //end i_subG_0
			}//end if (ig_0=jg_0)	
		}    
	} //end j_subG_0

}

void get_G_old_struct_matrix_memory_file(int tot_sub_vol, double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char wave_type, char* G_old_file_name){

	// ################### MATRICES STRUCTURE LOOPS ###########################
	// 3N X 3N Matrices structure loops for G^0 and A:
	double denom_NF, denom_IF ; // used in G^0_ij function
	double complex const_1, const_2, const_3;

	double eyeG_0[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};
	double complex G_oldValue;
	double a_j, part1ii;
	double complex part2ii, part2iiexp,part3ii;
	FILE * G_old_export = fopen(G_old_file_name, "wb"); 
	if (G_old_export == NULL) {
    	perror("Error opening binary file");
    	exit(1); // Exit with an error code
	}
	
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++)//
		{
			if (ig_0!=jg_0) // eq. 25 from Walter et al., PRB 2021
			{
			const_1 = cexp(k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0]*I)/(4.*pi*modulo_r_i_j[ig_0][jg_0]); 
			denom_NF = epsilon_ref*pow(k_0*modulo_r_i_j[ig_0][jg_0],2);
			denom_IF = k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0];
			const_2 = (1. - 1./denom_NF + 1.*I/denom_IF ) ;
			const_3 = (1. - 3./denom_NF + 3.*I/denom_IF) ;
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					//G_old[ig_0_2d][jg_0_2d] = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
					G_oldValue = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
					int position_old_ij = 9*tot_sub_vol*ig_0+3*tot_sub_vol*i_subG_0+3*jg_0+j_subG_0; //seems to be correct
					int position_old_ji = 9*tot_sub_vol*jg_0+3*tot_sub_vol*j_subG_0+3*ig_0+i_subG_0; //seems to be correct
					//printf("position =%d , ",position); 
					
					fseek(G_old_export, position_old_ij * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
        			size_t elements_ij_read = fwrite(&G_oldValue, sizeof(double complex), 1, G_old_export); // Read the matrix data from the binary file into the struct
					if (elements_ij_read != 1) { // Handle error or add debugging information
            			printf("Error reading data at position %d, when i!=m while writing G_new_ij\n", position_old_ij);
						exit(1); // Exit with an error code
        			}
					//G_old[jg_0_2d][ig_0_2d] = G_old[ig_0_2d][jg_0_2d];
					fseek(G_old_export, position_old_ji * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
        			size_t element_ji_read = fwrite(&G_oldValue, sizeof(double complex), 1, G_old_export); // Read the matrix data from the binary file into the struct
					if (element_ji_read != 1) { // Handle error or add debugging information
            			printf("Error reading data at position %d, when i!=m while writing G_new_ij\n", position_old_ji);
						exit(1); // Exit with an error code
        			}

				}    
			}
			}//end if (ig_0!=jg_0)
			else //if (ig_0==jg_0) // eq. 26 from Walter et al., PRB 2021 
			{
				a_j = a_j_function(delta_V_vector[ig_0], pi);
				part1ii = 1./(3.*delta_V_vector[ig_0]*epsilon_ref*pow(k_0,2));
				part2ii = a_j*k_0*sqrt(epsilon_ref)*I; // com i term 
				part2iiexp = cexp(0. + a_j*k_0*sqrt(epsilon_ref)*I); 
				part3ii = part2iiexp*(1-part2ii) - 1. ; // part3ii is inside brackets
				for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
				{
					int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
					{
						int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
						//G_old[ig_0_2d][jg_0_2d] = eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.); 
						int position_old_ij = 9*tot_sub_vol*ig_0+3*tot_sub_vol*i_subG_0+3*jg_0+j_subG_0; //seems to be correct
						G_oldValue = eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.); 
						
						fseek(G_old_export, position_old_ij * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
        				size_t elements_ij_read = fwrite(&G_oldValue, sizeof(double complex), 1, G_old_export); // Read the matrix data from the binary file into the struct
						if (elements_ij_read != 1) { // Handle error or add debugging information
            				printf("Error reading data at position %d, when i!=m while writing G_new_ij\n", position_old_ij);
							exit(1); // Exit with an error code
        				}

					} //end j_subG_0
				} //end i_subG_0
			} // end if (ig_0=jg_0)	
		}    // end jg_0
	} // end ig_0

	fclose(G_old_export);

}

void get_G0_matrix_memory_2D(int tot_sub_vol, double complex G_0[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char wave_type){

	// ################### MATRICES STRUCTURE LOOPS ###########################
	// 3N X 3N Matrices structure loops for G^0 and A:
	double denom_NF, denom_IF ; // used in G^0_ij function
	double complex const_1, const_2, const_3;

	double eyeG_0[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};

	double a_j, part1ii;
	double complex part2ii, part2iiexp,part3ii;
	
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	//for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++)//
	//for (int jg_0 = 0; jg_0 < tot_sub_vol-1; jg_0++) //tot_sub_vol
	{
		for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++)//
		//for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++)//
		//for (int ig_0 = jg_0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
		{
			if (ig_0!=jg_0) // eq. 25 from Walter et al., PRB 2021
			{
			const_1 = cexp(k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0]*I)/(4.*pi*modulo_r_i_j[ig_0][jg_0]); 
			denom_NF = epsilon_ref*pow(k_0*modulo_r_i_j[ig_0][jg_0],2);
			denom_IF = k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0];
			const_2 = (1. - 1./denom_NF + 1.*I/denom_IF ) ;
			const_3 = (1. - 3./denom_NF + 3.*I/denom_IF) ;
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					G_0[ig_0_2d][jg_0_2d] = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
					G_0[jg_0_2d][ig_0_2d] = G_0[ig_0_2d][jg_0_2d];
				}    
			}
			}//end if (ig_0!=jg_0)
			else //if (ig_0==jg_0) // eq. 26 from Walter et al., PRB 2021 
			{
				a_j = a_j_function(delta_V_vector[ig_0], pi);
				part1ii = 1./(3.*delta_V_vector[ig_0]*epsilon_ref*pow(k_0,2));
				part2ii = a_j*k_0*sqrt(epsilon_ref)*I; // com i term 
				part2iiexp = cexp(0. + a_j*k_0*sqrt(epsilon_ref)*I); 
				part3ii = part2iiexp*(1-part2ii) - 1. ; // part3ii is inside brackets
				for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
				{					
					int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
					{
						int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
						G_0[ig_0_2d][jg_0_2d] = eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.); 
					}
				} //end i_subG_0
			}//end if (ig_0=jg_0)	
		}    
	} //end j_subG_0

}

void get_A_matrix_2D(int tot_sub_vol, double complex G_0[3*tot_sub_vol][3*tot_sub_vol], double complex A[3*tot_sub_vol][3*tot_sub_vol], double k_0, double complex alpha_0[tot_sub_vol]){ //,char wave_type
	
	double complex (*alpha_0_matrix)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*alpha_0_matrix));
	if (alpha_0_matrix == NULL){
		printf("Failure with memory when generating A. Use iterative solver");
		exit(1);
	} 
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++)//
	{
		for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++)//
		{
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					if (ig_0==jg_0 && i_subG_0==j_subG_0)
					{
						alpha_0_matrix[ig_0_2d][jg_0_2d] = alpha_0[ig_0];
					}
					else
					{
						alpha_0_matrix[ig_0_2d][jg_0_2d] = 0.;
					}
				}
			}
		}
	}
	
	double complex (*prod)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*prod));
	if (prod == NULL){
		printf("Failure with memory when generating A. Use iterative solver");
		exit(1);
	} 

	double complex alpha_parameter = pow(k_0,2) + 0.0 * I;  // Scaling factor for A*B
	double complex beta_parameter = 0.0 + 0.0 * I;   // Scaling factor for C
	
	//cblas_zgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans, 3*tot_sub_vol, 3*tot_sub_vol, 3*tot_sub_vol,&alpha_parameter, G_0_matrix_2D,  3*tot_sub_vol, alpha_0_matrix_2D, 3*tot_sub_vol,&beta_parameter,prod, 3*tot_sub_vol);
	cblas_zgemm(
		CblasRowMajor,CblasNoTrans, CblasNoTrans, 
		3*tot_sub_vol, 3*tot_sub_vol, 3*tot_sub_vol,  	// Dimensions of matrices
		&alpha_parameter, 								// Scaling factor for A*B
		G_0,  3*tot_sub_vol, 					// Matrix A
		alpha_0_matrix, 3*tot_sub_vol, 				// Matrix B
		&beta_parameter, 								// Scaling factor for C
		prod, 3*tot_sub_vol								 // Matrix C
	);
	
	free(alpha_0_matrix);

	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++)//
	{
		for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++)//
		{
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					if (ig_0==jg_0 && i_subG_0==j_subG_0) // if i=j:
					{
						A[ig_0_2d][jg_0_2d] = 1. -  prod[ig_0_2d][jg_0_2d]; // eq. 26 from Lindsay's paper
					
					}	
					else
					{
						A[ig_0_2d][jg_0_2d] = 0. -  prod[ig_0_2d][jg_0_2d]; // eq. 26 from Lindsay's paper
					}
				} //end j_subG_0 
			} //end jg_0  
		}  //end i_subG_0     
	} //end ig_0    

free(prod);

}

void get_alpha_prod_for_A_2D(int tot_sub_vol, double complex G_0[3*tot_sub_vol][3*tot_sub_vol], double k_0, double complex alpha_0[tot_sub_vol], double complex prod[3*tot_sub_vol][3*tot_sub_vol]){
	
	double complex (*alpha_0_matrix)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*alpha_0_matrix));
	if (alpha_0_matrix == NULL){
		printf("Failure with memory when generating A. Use iterative solver");
		exit(1);
	} 
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++)//
	{
		for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++)//
		{
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					if (ig_0==jg_0 && i_subG_0==j_subG_0)
					{
						alpha_0_matrix[ig_0_2d][jg_0_2d] = alpha_0[ig_0];
					}
					else
					{
						alpha_0_matrix[ig_0_2d][jg_0_2d] = 0.;
					}
				}
			}
		}
	}
	
	double complex alpha_parameter = pow(k_0,2) + 0.0 * I;  // Scaling factor for A*B
	double complex beta_parameter = 0.0 + 0.0 * I;   // Scaling factor for C
	
	//cblas_zgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans, 3*tot_sub_vol, 3*tot_sub_vol, 3*tot_sub_vol,&alpha_parameter, G_0_matrix_2D,  3*tot_sub_vol, alpha_0_matrix_2D, 3*tot_sub_vol,&beta_parameter,prod, 3*tot_sub_vol);
	cblas_zgemm(
		CblasRowMajor,CblasNoTrans, CblasNoTrans, 
		3*tot_sub_vol, 3*tot_sub_vol, 3*tot_sub_vol,  	// Dimensions of matrices
		&alpha_parameter, 								// Scaling factor for A*B
		G_0,  3*tot_sub_vol, 					// Matrix A
		alpha_0_matrix, 3*tot_sub_vol, 				// Matrix B
		&beta_parameter, 								// Scaling factor for C
		prod, 3*tot_sub_vol								 // Matrix C
	);
	
	free(alpha_0_matrix);
}

void fill_A_matrix_2D(int tot_sub_vol, double complex A[3*tot_sub_vol][3*tot_sub_vol], double complex prod[3*tot_sub_vol][3*tot_sub_vol]){ 
	
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++)//
	{
		for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++)//
		{
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					if (ig_0==jg_0 && i_subG_0==j_subG_0) // if i=j:
					{
						A[ig_0_2d][jg_0_2d] = 1. -  prod[ig_0_2d][jg_0_2d]; // eq. 26 from Lindsay's paper
					
					}	
					else
					{
						A[ig_0_2d][jg_0_2d] = 0. -  prod[ig_0_2d][jg_0_2d]; // eq. 26 from Lindsay's paper
					}
				} //end j_subG_0 
			} //end jg_0  
		}  //end i_subG_0     
	} //end ig_0  
}


void setup_G_0_matrices(int tot_sub_vol, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double R[][3], char multithread){
	
	double (*r)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*r));
	if (r == NULL){
		printf("Failure with memory=%ld before spectral analysis during setup_G_0_matrices. ",get_mem_usage());
		//return 1;
		exit(1);
	}
	/*
	double (*abs_r_ij)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*abs_r_ij));
	if (abs_r_ij == NULL){
		printf("Failure with memory=%ld before spectral analysis during setup_G_0_matrices. ",get_mem_usage());
		//return 1;
		exit(1);
	}
	*/

	for (int i=0; i<tot_sub_vol; i++)
	{   
		for (int j=0; j<tot_sub_vol; j++)
		{
			for (int i_alpha=0; i_alpha<3; i_alpha++)
			{
				r[i][j][i_alpha]= R[i][i_alpha] - R[j][i_alpha] ; //general case; can be a scalar 
				//abs_r_ij[i][j][i_alpha] = fabs(r[i][j][i_alpha]); // can be a scalar
				double abs_r_ij = fabs(r[i][j][i_alpha]);
				//modulo_r_i_j[i][j] += pow(abs_r_ij[i][j][i_alpha],2); // change the name in this step, can be a scalar
				modulo_r_i_j[i][j] += pow(abs_r_ij,2); // change the name in this step, can be a scalar 
			}
			modulo_r_i_j[i][j] = sqrt(modulo_r_i_j[i][j]);
			printf("%e\n",modulo_r_i_j[i][j]);
		}
	}
	//free(abs_r_ij);

	
	double (*unit_r_ij)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*unit_r_ij));
	if (unit_r_ij == NULL){
		printf("Failure with memory=%ld before spectral analysis during setup_G_0_matrices. ",get_mem_usage());
		//return 1;
		exit(1);
	}
	double complex (*unit_conj_r_ij)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*unit_conj_r_ij));
	if (unit_conj_r_ij == NULL){
		printf("Failure with memory=%ld before spectral analysis during setup_G_0_matrices. ",get_mem_usage());
		//return 1;
		exit(1);
	}
	
	#pragma omp parallel for // if (multithread == 'Y')
	for (int i_i = 0; i_i < tot_sub_vol; i_i++)
	{ 
		for (int i_j = 0; i_j < tot_sub_vol; i_j++)
		{ 
			for (int i_alpha = 0; i_alpha < 3; i_alpha++)
			{  
				unit_r_ij[i_i][i_j][i_alpha] = r[i_i][i_j][i_alpha]/modulo_r_i_j[i_i][i_j]; // ˆr -- unit distance 
				//double transpose = r[i_i][i_j][i_alpha]/modulo_r_i_j[i_i][i_j]; // https://www.programiz.com/c-programming/examples/matrix-transpose        
				//unit_conj_r_ij[i_i][i_j][i_alpha] = conj(transpose);// r^+ Conjugate transpose unit distance 
				//double complex unit_conj_r_ij = conj(transpose); 
				unit_conj_r_ij[i_i][i_j][i_alpha] = conj(r[i_i][i_j][i_alpha]/modulo_r_i_j[i_i][i_j]);// r^+ Conjugate transpose unit distance 
			}
			for (int i_alpha = 0; i_alpha < 3; i_alpha++)
			{
				for (int j_alpha = 0; j_alpha < 3; j_alpha++)
				{ 
					r_i_j_outer_r_i_j[i_i][i_j][i_alpha][j_alpha] = unit_r_ij[i_i][i_j][i_alpha]*unit_conj_r_ij[i_i][i_j][j_alpha];
				}
			}	
		}
	}
	free(r);

	free(unit_r_ij);
	free(unit_conj_r_ij);
}


void set_up_get_G0_2D(int tot_sub_vol, double complex G_0[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi, double epsilon_ref, double delta_V_vector[tot_sub_vol], double R[][3]){
//void set_up_get_G0_2D(int tot_sub_vol, double complex G_0[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi, double epsilon_ref, double delta_V_vector[tot_sub_vol],char multithread, double R[][3]){

	//printf("set-up and calculation of free space Green's function\n");
	// ################### MATRICES STRUCTURE LOOPS ###########################
	// 3N X 3N Matrices structure loops for G^0 and A:
	//double denom_NF, denom_IF ; // used in G^0_ij function
	//double complex const_1, const_2, const_3;
	//double a_j, part1ii;
	//double complex part2ii, part2iiexp,part3ii;
	double eyeG_0[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};
	double (*r_ij)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*r_ij));
	if (r_ij == NULL){
		printf("Failure with memory=%ld during setup_G_0_matrices. ",get_mem_usage());
		//return 1;
		exit(1);
	}
	#pragma omp parallel for // if (multithread == 'Y')	// PARALLELIZE HERE
	//for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++)//
	//for (int jg_0 = 0; jg_0 < tot_sub_vol-1; jg_0++) //tot_sub_vol
	{
		//for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++)//
		//for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++)//
		for (int ig_0 = jg_0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
		{
			if (ig_0!=jg_0) // eq. 25 from Walter et al., PRB 2021
			{
				//printf("(i=%d j=%d)",ig_0,jg_0);
				double r_ij_mag;
				double r_ij_mag_pre=0;
				for (int i_alpha=0; i_alpha<3; i_alpha++)
				{
					//printf("%e, ",R[ig_0][i_alpha]);
					r_ij[ig_0][jg_0][i_alpha]= R[ig_0][i_alpha] - R[jg_0][i_alpha] ; //general case; can be a scalar 
					//abs_r_ij[i][j][i_alpha] = fabs(r[i][j][i_alpha]); // can be a scalar
					double abs_r_ij = fabs(r_ij[ig_0][jg_0][i_alpha]);
					//modulo_r_i_j[i][j] += pow(abs_r_ij[i][j][i_alpha],2); // change the name in this step, can be a scalar
					r_ij_mag_pre += pow(abs_r_ij,2); // change the name in this step, can be a scalar 
				}
				//printf("\n");
				r_ij_mag = sqrt(r_ij_mag_pre);
				//printf("%e\n",r_ij);
				double r_ij_hat[3];
				double complex conj_r_ij_hat[3];
				for (int i_alpha = 0; i_alpha < 3; i_alpha++)
				{  
					r_ij_hat[i_alpha] = r_ij[ig_0][jg_0][i_alpha]/r_ij_mag; // ˆr -- unit distance 
					conj_r_ij_hat[i_alpha] = conj(r_ij[ig_0][jg_0][i_alpha]/r_ij_mag);// r^+ Conjugate transpose unit distance 
				}
				//printf("%e %e %e\n",r_ij_hat[0],r_ij_hat[1],r_ij_hat[2]);
			double complex r_i_j_outer_r_i_j_t[3][3];
			for (int i_alpha = 0; i_alpha < 3; i_alpha++)
			{
				for (int j_alpha = 0; j_alpha < 3; j_alpha++)
				{ 
					r_i_j_outer_r_i_j_t[i_alpha][j_alpha] = r_ij_hat[i_alpha]*conj_r_ij_hat[j_alpha];
					//printf("r_i_j_outer_r_i_j_t[%d][%d]= %e+i%e , ",i_alpha,j_alpha, creal(r_i_j_outer_r_i_j_t[i_alpha][j_alpha]), cimag(r_i_j_outer_r_i_j_t[i_alpha][j_alpha]));
				}
				//printf("\n");
			}
			double complex const_1 = cexp(k_0*sqrt(epsilon_ref)*r_ij_mag*I)/(4.*pi*r_ij_mag); 
			//double denom_NF = epsilon_ref*pow(k_0*r_ij_mag,2);
			//double denom_IF = k_0*sqrt(epsilon_ref)*r_ij_mag;
			//double complex const_2 = (1. - 1./denom_NF + 1.*I/denom_IF ) ;
			//double complex const_3 = (1. - 3./denom_NF + 3.*I/denom_IF) ;
			double complex const_2 = (1. - 1./(epsilon_ref*pow(k_0*r_ij_mag,2)) + 1.*I/(k_0*sqrt(epsilon_ref)*r_ij_mag) ) ;
			double complex const_3 = (1. - 3./(epsilon_ref*pow(k_0*r_ij_mag,2)) + 3.*I/(k_0*sqrt(epsilon_ref)*r_ij_mag)) ;
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					G_0[ig_0_2d][jg_0_2d] = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j_t[i_subG_0][j_subG_0]));  
					G_0[jg_0_2d][ig_0_2d] = G_0[ig_0_2d][jg_0_2d];
					//printf("G_0[%d][%d]= %e+i%e , ",ig_0_2d,jg_0_2d, creal(G_0[ig_0_2d][jg_0_2d]), cimag(G_0[ig_0_2d][jg_0_2d]));
				}    
			}
			}//end if (ig_0!=jg_0)
			else //if (ig_0==jg_0) // eq. 26 from Walter et al., PRB 2021 
			{
				double a_j = a_j_function(delta_V_vector[ig_0], pi);
				double part1ii = 1./(3.*delta_V_vector[ig_0]*epsilon_ref*pow(k_0,2));
				double complex part2ii = a_j*k_0*sqrt(epsilon_ref)*I; // com i term 
				double complex part2iiexp = cexp(0. + a_j*k_0*sqrt(epsilon_ref)*I); 
				double complex part3ii = part2iiexp*(1-part2ii) - 1. ; // part3ii is inside brackets
				for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
				{					
					int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
					{
						int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
						G_0[ig_0_2d][jg_0_2d] = eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.); 
						//printf(" G_0[%d][%d]= %e+i%e , ",ig_0_2d,jg_0_2d,creal(G_0[ig_0_2d][jg_0_2d]), cimag(G_0[ig_0_2d][jg_0_2d]));
					}
				} //end i_subG_0
			}//end if (ig_0==jg_0)	
			//printf("\n \n"); 
		}   
	} //end j_subG_0
	
	free(r_ij);
}


void get_G_old_struct_matrix_memory(int tot_sub_vol, double complex G_old[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char multithread){

	// ################### MATRICES STRUCTURE LOOPS ###########################
	// 3N X 3N Matrices structure loops for G^0 and A:
	//double denom_NF, denom_IF ; // used in G^0_ij function
	//double complex const_1, const_2, const_3;
	//double a_j, part1ii;
	// double complex part2ii, part2iiexp,part3ii;
	double eyeG_0[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};
	#pragma omp parallel for // if (multithread == 'Y')	// PARALLELIZE HERE
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++)//
		{
			if (ig_0!=jg_0) // eq. 25 from Walter et al., PRB 2021
			{
			double complex const_1 = cexp(k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0]*I)/(4.*pi*modulo_r_i_j[ig_0][jg_0]); 
			double denom_NF = epsilon_ref*pow(k_0*modulo_r_i_j[ig_0][jg_0],2);
			double denom_IF = k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0];
			double complex const_2 = (1. - 1./denom_NF + 1.*I/denom_IF ) ;
			double complex const_3 = (1. - 3./denom_NF + 3.*I/denom_IF) ;
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					//G_old[ig_0][jg_0][i_subG_0][j_subG_0]=const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
					//G_old[jg_0][ig_0][j_subG_0][i_subG_0]=G_old[ig_0][jg_0][i_subG_0][j_subG_0];
					G_old[ig_0_2d][jg_0_2d] = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
					G_old[jg_0_2d][ig_0_2d] = G_old[ig_0_2d][jg_0_2d];
				}    
			}
			}//end if (ig_0!=jg_0)
			else //if (ig_0==jg_0) // eq. 26 from Walter et al., PRB 2021 
			{
				double a_j = a_j_function(delta_V_vector[ig_0], pi);
				double part1ii = 1./(3.*delta_V_vector[ig_0]*epsilon_ref*pow(k_0,2));
				double complex part2ii = a_j*k_0*sqrt(epsilon_ref)*I; // com i term 
				double complex part2iiexp = cexp(0. + a_j*k_0*sqrt(epsilon_ref)*I); 
				double complex part3ii = part2iiexp*(1-part2ii) - 1. ; // part3ii is inside brackets
				for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
				{
					int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
					{
						int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
						//G_old[ig_0][jg_0][i_subG_0][j_subG_0]=eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.);
						G_old[ig_0_2d][jg_0_2d] = eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.); 
					}
				} //end i_subG_0
			}//end if (ig_0=jg_0)	
		}    
	} //end j_subG_0

}


void get_G_old_struct_matrix_memory_file(int tot_sub_vol, double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char multithread, char* G_old_file_name){

	// ################### MATRICES STRUCTURE LOOPS ###########################
	// 3N X 3N Matrices structure loops for G^0 and A:
	//double denom_NF, denom_IF ; // used in G^0_ij function
	//double complex const_1, const_2, const_3;
	//double complex G_oldValue;
	//double a_j, part1ii;
	//double complex part2ii, part2iiexp,part3ii;
	double eyeG_0[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};
	FILE * G_old_export = fopen(G_old_file_name, "wb"); 
	if (G_old_export == NULL) {
    	perror("Error opening binary file");
    	exit(1); // Exit with an error code
	}
	#pragma omp parallel for //if (multithread == 'Y')	// PARALLELIZE HERE
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++)//
		{
			if (ig_0!=jg_0) // eq. 25 from Walter et al., PRB 2021
			{
			double complex const_1 = cexp(k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0]*I)/(4.*pi*modulo_r_i_j[ig_0][jg_0]); 
			double denom_NF = epsilon_ref*pow(k_0*modulo_r_i_j[ig_0][jg_0],2);
			double denom_IF = k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0];
			double complex const_2 = (1. - 1./denom_NF + 1.*I/denom_IF ) ;
			double complex const_3 = (1. - 3./denom_NF + 3.*I/denom_IF) ;
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					//G_old[ig_0_2d][jg_0_2d] = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
					double complex G_oldValue = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
					int position_old_ij = 9*tot_sub_vol*ig_0+3*tot_sub_vol*i_subG_0+3*jg_0+j_subG_0; //seems to be correct
					int position_old_ji = 9*tot_sub_vol*jg_0+3*tot_sub_vol*j_subG_0+3*ig_0+i_subG_0; //seems to be correct
					//printf("position =%d , ",position); 
					
					fseek(G_old_export, position_old_ij * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
        			size_t elements_ij_read = fwrite(&G_oldValue, sizeof(double complex), 1, G_old_export); // Read the matrix data from the binary file into the struct
					if (elements_ij_read != 1) { // Handle error or add debugging information
            			printf("Error reading data at position %d, when i!=m while writing G_new_ij\n", position_old_ij);
						exit(1); // Exit with an error code
        			}
					//G_old[jg_0_2d][ig_0_2d] = G_old[ig_0_2d][jg_0_2d];
					fseek(G_old_export, position_old_ji * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
        			size_t element_ji_read = fwrite(&G_oldValue, sizeof(double complex), 1, G_old_export); // Read the matrix data from the binary file into the struct
					if (element_ji_read != 1) { // Handle error or add debugging information
            			printf("Error reading data at position %d, when i!=m while writing G_new_ij\n", position_old_ji);
						exit(1); // Exit with an error code
        			}

				}    
			}
			}//end if (ig_0!=jg_0)
			else //if (ig_0==jg_0) // eq. 26 from Walter et al., PRB 2021 
			{
				double a_j = a_j_function(delta_V_vector[ig_0], pi);
				double part1ii = 1./(3.*delta_V_vector[ig_0]*epsilon_ref*pow(k_0,2));
				double complex part2ii = a_j*k_0*sqrt(epsilon_ref)*I; // com i term 
				double complex part2iiexp = cexp(0. + a_j*k_0*sqrt(epsilon_ref)*I); 
				double complex part3ii = part2iiexp*(1-part2ii) - 1. ; // part3ii is inside brackets
				for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
				{
					int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
					{
						int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
						//G_old[ig_0_2d][jg_0_2d] = eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.); 
						int position_old_ij = 9*tot_sub_vol*ig_0+3*tot_sub_vol*i_subG_0+3*jg_0+j_subG_0; //seems to be correct
						double complex G_oldValue = eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.); 
						
						fseek(G_old_export, position_old_ij * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
        				size_t elements_ij_read = fwrite(&G_oldValue, sizeof(double complex), 1, G_old_export); // Read the matrix data from the binary file into the struct
						if (elements_ij_read != 1) { // Handle error or add debugging information
            				printf("Error reading data at position %d, when i!=m while writing G_new_ij\n", position_old_ij);
							exit(1); // Exit with an error code
        				}

					} //end j_subG_0
				} //end i_subG_0
			} // end if (ig_0=jg_0)	
		}    // end jg_0
	} // end ig_0

	fclose(G_old_export);

}

void get_G0_matrix_memory_2D(int tot_sub_vol, double complex G_0[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi, double epsilon_ref, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double delta_V_vector[tot_sub_vol],char multithread){

	// ################### MATRICES STRUCTURE LOOPS ###########################
	// 3N X 3N Matrices structure loops for G^0 and A:
	//double denom_NF, denom_IF ; // used in G^0_ij function
	//double complex const_1, const_2, const_3;
	//double a_j, part1ii;
	//double complex part2ii, part2iiexp,part3ii;
	double eyeG_0[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};
	#pragma omp parallel for //if (multithread == 'Y')	// PARALLELIZE HERE
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	//for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++)//
	//for (int jg_0 = 0; jg_0 < tot_sub_vol-1; jg_0++) //tot_sub_vol
	{
		for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++)//
		//for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++)//
		//for (int ig_0 = jg_0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
		{
			if (ig_0!=jg_0) // eq. 25 from Walter et al., PRB 2021
			{
			double complex const_1 = cexp(k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0]*I)/(4.*pi*modulo_r_i_j[ig_0][jg_0]); 
			double denom_NF = epsilon_ref*pow(k_0*modulo_r_i_j[ig_0][jg_0],2);
			double denom_IF = k_0*sqrt(epsilon_ref)*modulo_r_i_j[ig_0][jg_0];
			double complex const_2 = (1. - 1./denom_NF + 1.*I/denom_IF ) ;
			double complex const_3 = (1. - 3./denom_NF + 3.*I/denom_IF) ;
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					G_0[ig_0_2d][jg_0_2d] = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
					G_0[jg_0_2d][ig_0_2d] = G_0[ig_0_2d][jg_0_2d];
					printf("G_0[%d][%d]= %e+i%e , ",ig_0_2d,jg_0_2d, creal(G_0[ig_0_2d][jg_0_2d]), cimag(G_0[ig_0_2d][jg_0_2d]));
				}    
			}
			}//end if (ig_0!=jg_0)
			else //if (ig_0==jg_0) // eq. 26 from Walter et al., PRB 2021 
			{
				double a_j = a_j_function(delta_V_vector[ig_0], pi);
				double part1ii = 1./(3.*delta_V_vector[ig_0]*epsilon_ref*pow(k_0,2));
				double complex part2ii = a_j*k_0*sqrt(epsilon_ref)*I; // com i term 
				double complex part2iiexp = cexp(0. + a_j*k_0*sqrt(epsilon_ref)*I); 
				double complex part3ii = part2iiexp*(1-part2ii) - 1. ; // part3ii is inside brackets
				for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
				{					
					int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
					{
						int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
						G_0[ig_0_2d][jg_0_2d] = eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.); 
						printf(" G_0[%d][%d]= %e+i%e , ",ig_0_2d,jg_0_2d,creal(G_0[ig_0_2d][jg_0_2d]), cimag(G_0[ig_0_2d][jg_0_2d]));
					}
				} //end i_subG_0
			}//end if (ig_0==jg_0)	
			printf("\n \n"); 
		}   
	} //end j_subG_0

}



void get_A_matrix_2D(int tot_sub_vol, double complex G_0[3*tot_sub_vol][3*tot_sub_vol], double complex A[3*tot_sub_vol][3*tot_sub_vol], double k_0, double complex alpha_0[tot_sub_vol]){ 
//void get_A_matrix_2D(int tot_sub_vol, double complex G_0[3*tot_sub_vol][3*tot_sub_vol], double complex A[3*tot_sub_vol][3*tot_sub_vol], double k_0, double complex alpha_0[tot_sub_vol], char multithread){ 
	
	double complex (*alpha_0_matrix)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*alpha_0_matrix));
	if (alpha_0_matrix == NULL){
		printf("Failure with memory when generating A. Use iterative solver");
		exit(1);
	} 
	#pragma omp parallel for // if (multithread == 'Y')	// PARALLELIZE HERE
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++)//
	{
		for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++)//
		//for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++)//
		{
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					if (ig_0==jg_0 && i_subG_0==j_subG_0)
					{
						alpha_0_matrix[ig_0_2d][jg_0_2d] = alpha_0[ig_0];
						//printf("%e+i%e\n", creal(alpha_0[ig_0]),cimag(alpha_0[ig_0]));
					}
					else
					{
						alpha_0_matrix[ig_0_2d][jg_0_2d] = 0.;
						alpha_0_matrix[jg_0_2d][ig_0_2d] = 0.;
					}
				}
			}
		}
	}
	
	double complex (*prod)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*prod));
	if (prod == NULL){
		printf("Failure with memory when generating A. Use iterative solver");
		exit(1);
	} 

	double complex alpha_parameter = pow(k_0,2) + 0.0 * I;  // Scaling factor for A*B
	double complex beta_parameter = 1.0 + 0.0 * I;   // Scaling factor for C
	
	//cblas_zgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans, 3*tot_sub_vol, 3*tot_sub_vol, 3*tot_sub_vol,&alpha_parameter, G_0_matrix_2D,  3*tot_sub_vol, alpha_0_matrix_2D, 3*tot_sub_vol,&beta_parameter,prod, 3*tot_sub_vol);
	cblas_zgemm(
		CblasRowMajor,CblasNoTrans, CblasNoTrans, 
		//CblasColMajor,CblasNoTrans, CblasNoTrans, 
		3*tot_sub_vol, 3*tot_sub_vol, 3*tot_sub_vol,  	// Dimensions of matrices
		&alpha_parameter, 								// Scaling factor for A*B
		G_0,  3*tot_sub_vol, 					// Matrix A
		alpha_0_matrix, 3*tot_sub_vol, 				// Matrix B
		&beta_parameter, 								// Scaling factor for C
		prod, 3*tot_sub_vol								 // Matrix C
	);
	
	free(alpha_0_matrix);

	#pragma omp parallel for // if (multithread == 'Y')	// PARALLELIZE HERE
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++)//
	{
		for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++)//
		//for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++)//
		{
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					if (ig_0==jg_0 && i_subG_0==j_subG_0) // if i=j:
					{
						A[ig_0_2d][jg_0_2d] = 1. -  prod[ig_0_2d][jg_0_2d]; // eq. 26 from Lindsay's paper
						//printf(" %e+i%e , ",creal(A[ig_0_2d][jg_0_2d]), cimag(A[ig_0_2d][jg_0_2d]));
					
					}	
					else
					{
						A[ig_0_2d][jg_0_2d] = 0. -  prod[ig_0_2d][jg_0_2d]; // eq. 26 from Lindsay's paper
						A[jg_0_2d][ig_0_2d] = 0. -  prod[jg_0_2d][ig_0_2d]; // eq. 26 from Lindsay's paper
						//printf(" %e+i%e , ",creal(A[ig_0_2d][jg_0_2d]), cimag(A[ig_0_2d][jg_0_2d]));
					}
				} //end j_subG_0 
			} //end jg_0  
		}  //end i_subG_0 
			//printf("\n \n");    
	} //end ig_0    

free(prod);

}

