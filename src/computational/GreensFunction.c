#include "computational/GreensFunction.h"
#include <math.h>
#include "functions_DSGF.h"
#include <stdlib.h>
#include "computational/solvers/iterative_solver.h"
#include "mkl.h"

void set_up_get_G0_1D(int tot_sub_vol, double complex G_0[3*tot_sub_vol*3*tot_sub_vol], double k_0, double pi, double epsilon_ref, double delta_V_vector[tot_sub_vol], double R[][3]){
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
				//int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					//int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					long index_ij = j_subG_0+jg_0*3+(long)i_subG_0*3*tot_sub_vol+(long)	ig_0*3*tot_sub_vol*3;
					long index_ji = i_subG_0+ig_0*3+j_subG_0*3*tot_sub_vol+jg_0*3*tot_sub_vol*3;
					G_0[index_ij] = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j_t[i_subG_0][j_subG_0]));  
					G_0[index_ji] = G_0[index_ij];
					//G_0[ig_0_2d][jg_0_2d] = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j_t[i_subG_0][j_subG_0]));  
					//G_0[jg_0_2d][ig_0_2d] = G_0[ig_0_2d][jg_0_2d];
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
					//int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
					{
						//int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
						long index_ij = j_subG_0+jg_0*3+(long)i_subG_0*3*tot_sub_vol+(long)ig_0*3*tot_sub_vol*3;
						G_0[index_ij] = eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.);
						//G_0[ig_0_2d][jg_0_2d] = eyeG_0[i_subG_0][j_subG_0]*part1ii*(2.*part3ii-1.); 
						//printf(" G_0[%d][%d]= %e+i%e , ",ig_0_2d,jg_0_2d,creal(G_0[ig_0_2d][jg_0_2d]), cimag(G_0[ig_0_2d][jg_0_2d]));
					}
				} //end i_subG_0
			}//end if (ig_0==jg_0)	
			//printf("\n \n"); 
		}   
	} //end j_subG_0
	
	free(r_ij);
}


void set_up_get_G_old(int tot_sub_vol, double complex G_old[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi, double epsilon_ref, double delta_V_vector[tot_sub_vol], double R[][3]){
//void set_up_get_G_old(int tot_sub_vol, double complex G_old[3*tot_sub_vol][3*tot_sub_vol], double k_0, double pi, double epsilon_ref, double delta_V_vector[tot_sub_vol],char multithread, double R[][3]){

	//printf("set-up and calculation of free space Green's function\n");
	// ################### MATRICES STRUCTURE LOOPS ###########################
	// 3N X 3N Matrices structure loops for G^0 and A:
	//double denom_NF, denom_IF ; // used in G^0_ij function
	//double complex const_1, const_2, const_3;
	//double a_j, part1ii;
	// double complex part2ii, part2iiexp,part3ii;
	double eyeG_0[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};
	double (*r_ij)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*r_ij));
	if (r_ij == NULL){
		printf("Failure with memory=%ld during setup_G_0_matrices. ",get_mem_usage());
		//return 1;
		exit(1);
	}
	#pragma omp parallel for // if (multithread == 'Y')	// PARALLELIZE HERE
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++)//
		{
			if (ig_0!=jg_0) // eq. 25 from Walter et al., PRB 2021
			{
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
			double denom_NF = epsilon_ref*pow(k_0*r_ij_mag,2);
			double denom_IF = k_0*sqrt(epsilon_ref)*r_ij_mag;
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
					G_old[ig_0_2d][jg_0_2d] = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j_t[i_subG_0][j_subG_0]));  
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

void set_up_get_G_old_file(int tot_sub_vol, double k_0, double pi, double epsilon_ref, double delta_V_vector[tot_sub_vol],char multithread, char* G_old_file_name,double R[][3]){

	//printf("set-up and calculation of free space Green's function\n");
	// ################### MATRICES STRUCTURE LOOPS ###########################
	// 3N X 3N Matrices structure loops for G^0 and A:
	//double denom_NF, denom_IF ; // used in G^0_ij function
	//double complex const_1, const_2, const_3;
	//double complex G_oldValue;
	//double a_j, part1ii;
	//double complex part2ii, part2iiexp,part3ii;
	double eyeG_0[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};
	double (*r_ij)[tot_sub_vol][3] = calloc(tot_sub_vol, sizeof(*r_ij));
	if (r_ij == NULL){
		printf("Failure with memory=%ld during setup_G_0_matrices. ",get_mem_usage());
		//return 1;
		exit(1);
	}
	FILE * G_old_export = fopen(G_old_file_name, "wb"); 
	if (G_old_export == NULL) {
    	perror("Error opening binary file");
    	exit(1); // Exit with an error code
	}
	// #pragma omp parallel for if (multithread == 'Y')	// PARALLELIZE HERE
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{
		for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++)//
		{
			if (ig_0!=jg_0) // eq. 25 from Walter et al., PRB 2021
			{
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
			double denom_NF = epsilon_ref*pow(k_0*r_ij_mag,2);
			double denom_IF = k_0*sqrt(epsilon_ref)*r_ij_mag;
			double complex const_2 = (1. - 1./denom_NF + 1.*I/denom_IF ) ;
			double complex const_3 = (1. - 3./denom_NF + 3.*I/denom_IF) ;
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				//int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
				{
					//int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					//G_old[ig_0_2d][jg_0_2d] = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j[ig_0][jg_0][i_subG_0][j_subG_0]));  
					double complex G_oldValue = const_1*((const_2*eyeG_0[i_subG_0][j_subG_0])-(const_3*r_i_j_outer_r_i_j_t[i_subG_0][j_subG_0]));  
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
					//int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
					{
						//int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
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

void get_A_matrix_1D(int tot_sub_vol, double complex G_0[3*tot_sub_vol*3*tot_sub_vol], double complex A[3*tot_sub_vol*3*tot_sub_vol], double k_0, double complex alpha_0[tot_sub_vol]){ 
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
						long index_ij = j_subG_0+jg_0*3+(long)i_subG_0*3*tot_sub_vol+(long)ig_0*3*tot_sub_vol*3;
						A[index_ij] = 1. -  prod[ig_0_2d][jg_0_2d]; // eq. 26 from Lindsay's paper
						//A[ig_0_2d][jg_0_2d] = 1. -  prod[ig_0_2d][jg_0_2d]; // eq. 26 from Lindsay's paper
						//printf(" %e+i%e , ",creal(A[ig_0_2d][jg_0_2d]), cimag(A[ig_0_2d][jg_0_2d]));
					
					}	
					else
					{
						long index_ij = j_subG_0+jg_0*3+(long)i_subG_0*3*tot_sub_vol+(long)ig_0*3*tot_sub_vol*3;
						long index_ji = i_subG_0+ig_0*3+(long)j_subG_0*3*tot_sub_vol+(long)jg_0*3*tot_sub_vol*3;
						A[index_ij] = 0. -  prod[ig_0_2d][jg_0_2d]; // eq. 26 from Lindsay's paper
						A[index_ji] = 0. -  prod[jg_0_2d][ig_0_2d]; // eq. 26 from Lindsay's paper
						//A[ig_0_2d][jg_0_2d] = 0. -  prod[ig_0_2d][jg_0_2d]; // eq. 26 from Lindsay's paper
						//A[jg_0_2d][ig_0_2d] = 0. -  prod[jg_0_2d][ig_0_2d]; // eq. 26 from Lindsay's paper
						//printf(" %e+i%e , ",creal(A[ig_0_2d][jg_0_2d]), cimag(A[ig_0_2d][jg_0_2d]));
					}
				} //end j_subG_0 
			} //end jg_0  
		}  //end i_subG_0 
			//printf("\n \n");    
	} //end ig_0    

free(prod);
}
