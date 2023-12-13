#include "computational/solvers/iterative_solver.h"
#include <stdio.h>
#include <string.h>
#include "computational/GreensFunction.h"
#include <stdlib.h>
#include <lapacke.h> // for desktop
// #include <mkl_lapacke.h> // for CHPC
#include <cblas.h>
#include <stdbool.h>
#include <complex.h>
#include <math.h>
#include "functions_DSGF.h"
#include "file_utils.h"

//#include "mkl.h" // if uncommented, a series of warning are shown when compiling.

void A2d_solver(double complex epsilon, int mm, int tot_sub_vol, double delta_V, double complex G_sys_old[][3*tot_sub_vol], double complex A_2d[3][3], double k){
	double eyeA_2d[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};

	for(int minor_x = 0; minor_x < 3; minor_x++){
		int new_x = 3*mm + minor_x;
		for (int minor_y = 0; minor_y < 3; minor_y++){
			int new_y = 3*mm + minor_y;
            A_2d[minor_x][minor_y] = eyeA_2d[minor_x][minor_y] - pow(k,2)*delta_V*epsilon*G_sys_old[new_x][new_y];
			//A_2d[minor_x][minor_y] = eyeA_2d[minor_x][minor_y] - pow(k,2)*delta_V*epsilon*upperTriangularMatrix[index];
			//printf(" A_2d[%d,%d]= (%e + i%e) ",minor_x,minor_y, creal(A_2d[minor_x][minor_y]),cimag(A_2d[minor_x][minor_y]));
				
		}
		//printf("\n");
	}
}

void A_b_iterative_populator(int tot_sub_vol, double complex *A_iterative, double complex *b_iterative, double complex A_2d[3][3], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol],int mm, int jg_0){

	//int ipack=0;
	for (int row = 0; row < 3; row++)
	{
		int new_x = (3*mm + row);
		for(int col = 0; col < 3; col++) // 3D coordinate positions 
		{
			int new_y = (3*jg_0 + col);
			int index = col+ 3*row ; //9 * (new_x * tot_sub_vol + new_y) + 3 * row + col;
			A_iterative[index] = A_2d[row][col]; //A_2d[mm_2d][mm_2d]; //A[mm][mm][i_subG_0][j_subG_0];
			b_iterative[index] = G_sys_old[new_x][new_y]; //G_sys_old[mm][jg_0][i_subG_0][j_subG_0];
			//b_iterative[ipack] = upperTriangularMatrix[index];
			//ipack++;
		}    
	}
}

void G_sys_new_populator(int tot_sub_vol, int mm, int jg_0, double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], double complex b_iterative[]){
  	
	for (int mm_sub = 0; mm_sub < 3; mm_sub++) // 3D coordinate positions
	{
		int mm_2d = (3*mm + mm_sub);
		for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
		{
			int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
			int index = j_subG_0 + mm_sub*3; 
			G_sys_new[mm_2d][jg_0_2d] = b_iterative[index]; // stores G^1_11
			G_sys_new[jg_0_2d][mm_2d] = G_sys_new[mm_2d][jg_0_2d]; //reciprocity
		}  
	}
}

void remaining_pertubations(int tot_sub_vol, int mm, double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], double complex A_2d[3][3]){
//void remaining_pertubations(int tot_sub_vol, int mm, double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], double complex A_2d[3][3], char multithread){

	#pragma omp parallel for // if (multithread == 'Y')// PARALELLIZE HERE
	for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) // Only loop through remaining perturbations
	{
		double complex A_iterative[3*3];
		double complex b_iterative[3*3];
		A_b_iterative_populator(tot_sub_vol, A_iterative, b_iterative, A_2d, G_sys_old, mm, jg_0);
		//A_b_iterative_populator(tot_sub_vol, A_iterative, b_iterative, A_2d, size, upperTriangularMatrix, mm, jg_0);
		int info = LAPACKE_zgels(LAPACK_ROW_MAJOR,'N',3,3,3,A_iterative,3,b_iterative,3);
		G_sys_new_populator(tot_sub_vol, mm, jg_0, G_sys_new, b_iterative);
	} // end jg_0					
	
}

void offdiagonal_solver(int tot_sub_vol, int mm, double k, double complex alpha_0[], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol]){
//void offdiagonal_solver(int tot_sub_vol, int mm, double k, double complex alpha_0[], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], char multithread){
			
	//modifications from Martin Cuma		
	#pragma omp parallel for collapse(2) // if (multithread == 'Y')	// PARALLELIZE HERE	
	// Next, solve all systems of equations for ii not equal to mm	
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //lower triangular
	{ 
		for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
		{
			int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
		for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++) //upper triangular matrix
		{
		if(ig_0 != mm  && jg_0 != mm)
		{
			for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
			{
				int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
				double complex G_sys_prod = 0.;
                #pragma omp simd 
				for(int m_sub = 0;  m_sub < 3;  m_sub++)//loop for matricial multiplication
				{
					int mm_2d = (3*mm + m_sub);
					G_sys_prod += G_sys_old[ig_0_2d][mm_2d]*pow(k,2)*alpha_0[mm]*G_sys_new[mm_2d][jg_0_2d];
					//G_sys_prod += G_sys_old[ig_0_2d][mm_2d]*pow(k,2)*alpha_0[mm]*G_sys_new_mm[m_sub][jg_0_2d];
				}
				//int index = 9 * (submatrix_i * tot_sub_vol + submatrix_j) + 3 * row_in_submatrix + col_in_submatrix;
				//G_sys_new[ig_0_2d][jg_0_2d] = upperTriangularMatrix[index] + G_sys_prod; //alpha_0[mm]*[ig_0_2d][jg_0_2d]
				G_sys_new[ig_0_2d][jg_0_2d] = G_sys_old[ig_0_2d][jg_0_2d] + G_sys_prod; //alpha_0[mm]*[ig_0_2d][jg_0_2d]
				G_sys_new[jg_0_2d][ig_0_2d] = G_sys_new[ig_0_2d][jg_0_2d]; // reciprocity
			} // j_subG_0    
		} // jg_0   	
		}// i_subG_0   	        				
		} // if(ig_0 != mm)              	
	} // ig_0    

	/*
	//previous version
	#pragma omp parallel for collapse(2) // if (multithread == 'Y')	// PARALLELIZE HERE	
	// Next, solve all systems of equations for ii not equal to mm		
	//for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //complete system
	for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //upper triangular matrix
	{
		for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
		{
			int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
			//for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //complete system
			for (int ig_0 = jg_0; ig_0 < tot_sub_vol; ig_0++) //lower triangular
			{ 
		if(ig_0 != mm  && jg_0 != mm)
		{
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				
						double complex G_sys_prod = 0.;
						for(int m_sub = 0;  m_sub < 3;  m_sub++)//loop for matricial multiplication
						{
							int mm_2d = (3*mm + m_sub);
							//int  index_mm = 9 * (submatrix_i * mm_2d + submatrix_j) + 3 * row_in_submatrix + col_in_submatrix;
							//G_sys_prod += G_sys_old[ig_0_2d][mm_2d]*G_sys_new_transpose[jg_0_2d][mm_2d];	// old test, do not consider
							G_sys_prod += G_sys_old[ig_0_2d][mm_2d]*pow(k,2)*alpha_0[mm]*G_sys_new[mm_2d][jg_0_2d];
							//G_sys_prod += upperTriangularMatrix[index]*pow(k,2)*alpha_0[mm]*G_sys_new[mm_2d][jg_0_2d];
						}
						//int index = 9 * (submatrix_i * tot_sub_vol + submatrix_j) + 3 * row_in_submatrix + col_in_submatrix;
						//G_sys_new[ig_0_2d][jg_0_2d] = upperTriangularMatrix[index] + G_sys_prod; //alpha_0[mm]*[ig_0_2d][jg_0_2d]
						G_sys_new[ig_0_2d][jg_0_2d] = G_sys_old[ig_0_2d][jg_0_2d] + G_sys_prod; //alpha_0[mm]*[ig_0_2d][jg_0_2d]
						G_sys_new[jg_0_2d][ig_0_2d] = G_sys_new[ig_0_2d][jg_0_2d]; // reciprocity
					} // j_subG_0    
			} // jg_0   	
		}// i_subG_0   	        				
		} // if(ig_0 != mm)              	
	} // ig_0    
	
	*/

}

void update_G_old(int tot_sub_vol, double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol]){

	// new version from Martin Cuma
	#pragma omp parallel for collapse(2) // if (multithread == 'Y')	// PARALLELIZE HERE	
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //lower triangular
	{ 
		for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
		{
			int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
			memcpy(G_sys_old[ig_0_2d],G_sys_new[ig_0_2d],sizeof(double complex)*tot_sub_vol*3);
		}// j_subG_0   	                 	
	} // jg_0 

	/* // old version
	#pragma omp parallel for // if (multithread == 'Y')	// PARALLELIZE HERE	
	for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //upper triangular matrix
	{
		for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
		{
			int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
			for (int ig_0 = jg_0; ig_0 < tot_sub_vol; ig_0++) //lower triangular
			{ 
				for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
				{
					int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
						G_sys_old[ig_0_2d][jg_0_2d] = G_sys_new[ig_0_2d][jg_0_2d]; //alpha_0[mm]*[ig_0_2d][jg_0_2d]
						G_sys_old[jg_0_2d][ig_0_2d] = G_sys_old[ig_0_2d][jg_0_2d]; // reciprocity
					} // j_subG_0    
			} // jg_0   	
		}// i_subG_0   	                 	
	} // ig_0 
	*/	   

}


void core_solver(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[],double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol]){
//void core_solver(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[],double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], char multithread){

	double complex A_2d[3][3];
	for (int mm = 0; mm < tot_sub_vol; mm++) //tot_sub_vol
	{
		printf("%d - ",mm+1);
		double complex epsilon_s = (epsilon - epsilon_ref); // Scattering dielectric function
		A2d_solver(epsilon_s, mm, tot_sub_vol, delta_V_vector[mm], G_sys_old, A_2d, k);	
		remaining_pertubations(tot_sub_vol, mm, G_sys_old, G_sys_new, A_2d);
		
		offdiagonal_solver(tot_sub_vol, mm, k, alpha_0, G_sys_old, G_sys_new);
		update_G_old(tot_sub_vol, G_sys_old, G_sys_new); //new function to update G_old = G_new
	
		//remaining_pertubations(tot_sub_vol, mm, G_sys_old, G_sys_new, A_2d, multithread);
		//offdiagonal_solver(tot_sub_vol, mm, k, alpha_0, G_sys_old, G_sys_new, multithread);
		//memmove(G_sys_old,G_sys_new,3*tot_sub_vol*3*tot_sub_vol*sizeof(double complex));
		//memcpy(G_sys_old,G_sys_new,3*tot_sub_vol*3*tot_sub_vol*sizeof(double complex)); // Update G_old = G_new for next iteration.
	
	}//end mm loop 
}

void iterative_solver(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol],double k_0, double pi, double R[][3]){
//void iterative_solver(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol],double k_0, double pi,char multithread, double R[][3]){
//void iterative_solver(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol],double k_0, double pi,double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3],char multithread){
	double complex (*G_sys_old)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys_old));
	if (G_sys_old == NULL){
			printf("Failure with memory=%ld in iterative solver.",get_mem_usage());
			exit(1);
	} 
	set_up_get_G_old(tot_sub_vol, G_sys_old, k_0, pi, epsilon_ref, delta_V_vector, R);
	core_solver(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_sys, G_sys_old);
	
	//get_G_old_struct_matrix_memory(tot_sub_vol, G_sys_old, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, multithread);
	//set_up_get_G_old(tot_sub_vol, G_sys_old, k_0, pi, epsilon_ref, delta_V_vector, multithread, R);
	//core_solver(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_sys, G_sys_old, multithread);
	free(G_sys_old);

}

void iterative_solver_file_handler(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[],double k_0, double pi,char multithread, char* G_old_file_name, char* G_sys_file_name, double R[][3]){
//void iterative_solver_file_handler(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[],double k_0, double pi,double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3],char multithread, char* G_old_file_name, char* G_sys_file_name){

	//get_G_old_struct_matrix_memory_file(tot_sub_vol, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, multithread,G_old_file_name);
	set_up_get_G_old_file(tot_sub_vol, k_0, pi, epsilon_ref, delta_V_vector, multithread, G_old_file_name, R);

	double complex A_2d[3][3];
	for (int mm = 0; mm < tot_sub_vol; mm++) //tot_sub_vol
	{
		//printf("mm= %d \n",mm+1);
		printf("%d - ",mm+1);

		double complex epsilon_s = (epsilon - epsilon_ref); // Scattering dielectric function
		double eyeA_2d[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};
		double complex G_oldValue; //struct Matrix3x3 G_old_mm; //
		double complex A_2d[3][3]; //struct Matrix3x3 A_2d_test; 
		FILE * G_old_import = fopen(G_old_file_name, "rb"); 
		if (G_old_import == NULL) {
    		perror("Error opening binary file");
    		exit(1); // Exit with an error code
		}
		
		double complex cte = - pow(k, 2) * delta_V_vector[mm] * epsilon_s;
		//printf("%e + i %e ,%e + i %e ,%e + i %e \n %e + i %e ,%e + i %e ,%e + i %e \n %e + i %e ,%e + i %e , %e + i %e \n",creal(G_old_mm.xx),cimag(G_old_mm.xx),creal(G_old_mm.xy),cimag(G_old_mm.xy),creal(G_old_mm.xz),cimag(G_old_mm.xz),creal(G_old_mm.yx),cimag(G_old_mm.yx),creal(G_old_mm.yy),cimag(G_old_mm.yy),creal(G_old_mm.yz),cimag(G_old_mm.yz),creal(G_old_mm.zx),cimag(G_old_mm.zx),creal(G_old_mm.zy),cimag(G_old_mm.zy),creal(G_old_mm.zz),cimag(G_old_mm.zz));
		for (int row = 0; row < 3; row++)
		{
			for (int column = 0; column < 3; column++)
			{
				int position = (9*tot_sub_vol+3)*mm+3*tot_sub_vol*row+column;
				//printf("position =%d , ",position); 
				fseek(G_old_import, position * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
        		size_t elements_read = fread(&G_oldValue, sizeof(double complex), 1, G_old_import); // Read the matrix data from the binary file into the struct
				if (elements_read != 1) { // Handle error or add debugging information
            		printf("Error reading data at position %d, when i!=m while writing G_new_ij\n", position);
					exit(1); // Exit with an error code
        		}
				//printf("G_old[%d][%d] = %e + i %e , ",row, column,creal(G_oldValue),cimag(G_oldValue));
				A_2d[row][column] =  eyeA_2d[row][column] + cte* G_oldValue; //A_2d_test =  &eyeA_2d - pow(k, 2) * delta_V_vector[mm] * epsilon_s * &G_old_mm;
				//printf("A_2d[%d][%d] = %e + i %e , ",row, column,creal(A_2d[row][column]),cimag(A_2d[row][column]));
			}
			//printf("\n");
		}
		
		//This part calculates the first G_new for when i=m. It uses G_old, but we do not need all the terms. 
		double complex A_iterative[3*3];
		double complex b_iterative[3*3];
		#define MAX_LINE_LENGTH 1024 
		
		FILE * G_new_import_export = fopen(G_sys_file_name, "w+");
		if (G_new_import_export == NULL) {
    		perror("Error opening binary file");
    		exit(1); // Exit with an error code
		}
		// #pragma omp parallel for if (multithread == 'Y')// PARALELLIZE HERE
		for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) // Only loop through remaining perturbations
		{
			//int ipack=0;
			//G_old_import = fopen(G_old_file_name, "rb");
			for (int row = 0; row < 3; row++)
			{
				for(int col = 0; col < 3; col++) // 3D coordinate positions 
				{
					int position_mj = 9*tot_sub_vol*mm+3*jg_0+3*tot_sub_vol*row+col; //correct position 
					//int position = (9*tot_sub_vol+3)*mm+3*tot_sub_vol*row+col;
					//printf("position =%d , ",position_mj); 
					fseek(G_old_import, position_mj * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
        			//fseek(G_old_import, 0, position); // 0 bytes from the start of the file
					size_t elements_read = fread(&G_oldValue, sizeof(double complex), 1, G_old_import); // Read the matrix data from the binary file into the struct
					if (elements_read != 1) { // Handle error or add debugging information
            			printf("Error reading data at position %d, when i!=m while writing G_new_ij\n", position_mj);
						exit(1); // Exit with an error code
        			}
					int index = col+3*row; 
					A_iterative[index] = A_2d[row][col]; //A_2d[mm_2d][mm_2d]; //A[mm][mm][i_subG_0][j_subG_0];
					b_iterative[index] = G_oldValue;
					//ipack++;
				}    
			}
			//printf("\n");
			int info = LAPACKE_zgels(LAPACK_ROW_MAJOR,'N',3,3,3,A_iterative,3,b_iterative,3); // G_new using Linear inversion using LAPACK
			
			//b_iterative is stored directly into a file, according to the term's position.
			
			for (int mm_sub = 0; mm_sub < 3; mm_sub++) // 3D coordinate positions
			{
				int mm_2d = (3*mm + mm_sub);
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
				{
					int index = j_subG_0+ 3*mm_sub;
					int position_mj = 9*tot_sub_vol*mm+3*jg_0+3*tot_sub_vol*mm_sub+j_subG_0; //correct position  
					//printf("position_mj =%d , ",position_mj);
                    rewind(G_new_import_export); 
                	fseek(G_new_import_export, position_mj * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
					size_t elements_read_mj = fwrite(&b_iterative[index], sizeof(double complex), 1, G_new_import_export); // Read the matrix data from the binary file into the struct
					if (elements_read_mj != 1) {
            		    printf("Error reading data at position %d\n", position_mj); // Handle error or add debugging information
        			}
					// reciprocity
					int position_jm = 9*tot_sub_vol*jg_0+3*mm+3*tot_sub_vol*j_subG_0+mm_sub; //correct position  
					//printf("position_jm =%d , ",position_jm); 
					rewind(G_new_import_export);
					fseek(G_new_import_export, position_jm * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
					size_t elements_read_jm = fwrite(&b_iterative[index], sizeof(double complex), 1, G_new_import_export); // Read the matrix data from the binary file into the struct
					if (elements_read_jm != 1) {
            		    printf("Error reading data at position %d\n", position_jm); // Handle error or add debugging information
        			}
					//int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					//printf("G_new[%d][%d] = %e + i %e , ",mm_2d, jg_0_2d,creal(b_iterative[gpack]),cimag(b_iterative[gpack]));
				}  
				//printf("\n");
			}
		} // end jg_0
		
		//This part calculates the remaining G_new for when i!=m. It uses G_old, but we do not need all the terms.
		double complex G_new_ij_value, G_old_ij_value;
		double complex G_old_im_value,G_new_mj_value;
		
		// #pragma omp parallel for if (multithread == 'Y') // PARALELLIZE HERE
		for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //upper triangular matrix
		{
			for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
			{
				int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
				for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //complete system
				//for (int ig_0 = jg_0; ig_0 < tot_sub_vol; ig_0++) //lower triangular
				{ 
					//if(ig_0 != mm ) //condition without reciprocity
					if(ig_0 != mm && jg_0 != mm) //using reciprocity, some terms would be calculated twice. The updated condition, each term is calculated one
					{
						for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
						{
							int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
							int position_old_ij = 9*tot_sub_vol*ig_0+3*tot_sub_vol*i_subG_0+3*jg_0+j_subG_0; //seems to be correct

							fseek(G_old_import, position_old_ij * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
        					//fseek(G_old_import, 0, position); // 0 bytes from the start of the file
							size_t elements_read = fread(&G_old_ij_value, sizeof(double complex), 1, G_old_import); // Read the matrix data from the binary file into the struct
							if (elements_read != 1) { // Handle error or add debugging information
            					printf("Error reading data at position %d, when i!=m while reading G_new_ij\n", position_old_ij);
								exit(1); // Exit with an error code
        					}
							double complex G_sys_prod = 0.;
							for(int m_sub = 0;  m_sub < 3;  m_sub++)//loop for matricial multiplication
							{
								int mm_2d = (3*mm + m_sub);
								//printf("G_old_im_value = %e + i %e, ",creal(G_old_im_value),cimag(G_old_im_value));
								//printf("G_new_mj_value = %e + i %e; ",creal(G_new_mj_value[m_sub]),cimag(G_new_mj_value[m_sub]));
								//G_sys_prod += G_old_im_value*pow(k,2)*alpha_0[mm]*G_new_mj_value[m_sub];
								int position_old_im = 9*tot_sub_vol*ig_0+3*tot_sub_vol*i_subG_0+3*mm+m_sub; //seems to be correct
								fseek(G_old_import, position_old_im * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
        						size_t SGF_old_read = fread(&G_old_im_value, sizeof(double complex), 1, G_old_import); // Read the matrix data from the binary file into the struct
								if (SGF_old_read != 1) { // Handle error or add debugging information
            						printf("Error reading data at position %d, when i!=m while reading G_new_im\n", position_old_im);
									exit(1); // Exit with an error code
        						}
								int position_new_mj = 9*tot_sub_vol*mm+3*tot_sub_vol*m_sub+3*jg_0+j_subG_0; // correct
								fseek(G_new_import_export, position_new_mj * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
        						size_t SGF_new_read = fread(&G_new_mj_value, sizeof(double complex), 1, G_new_import_export); // Read the matrix data from the binary file into the struct
								if (SGF_new_read != 1) { // Handle error or add debugging information
            						printf("Error reading data at position %d, when i!=m while reading G_new_mj\n", position_new_mj);
									exit(1); // Exit with an error code
        						}
								//printf("G_new_mj_value = %e + i %e; ",creal(G_new_mj_value),cimag(G_new_mj_value));
								G_sys_prod += G_old_im_value*pow(k,2)*alpha_0[mm]*G_new_mj_value;
							}
							G_new_ij_value = G_old_ij_value + G_sys_prod; //alpha_0[mm]*[ig_0_2d][jg_0_2d]
							//printf("G_sys_new[%d][%d][%d][%d]=%e+i%e\n",ig_0,jg_0,i_subG_0,j_subG_0,creal(G_new_ij_value),cimag(G_new_ij_value));
							fseek(G_new_import_export, position_old_ij * sizeof(double complex), SEEK_SET);
							size_t SGF_new_ij_write = fwrite(&G_new_ij_value, sizeof(double complex), 1, G_new_import_export);
							if (SGF_new_ij_write != 1) { // Handle error or add debugging information
            					printf("Error reading data at position %d, when i!=m while writing G_new_ij\n", position_old_ij);
								exit(1); // Exit with an error code
        					}
						} // i_subG_0    
					} //  if(ig_0 != mm) 	
				}// ig_0 	      
				//printf("\n");  				
			} // j_subG_0                 	
		} // jg_0    
		
		fclose(G_old_import);
		fclose(G_new_import_export); // Close the file
		
		if (mm<tot_sub_vol-1) {
			
    		// Open the old file for writing
    		FILE* oldFile = fopen(G_old_file_name, "wb");
    		if (oldFile == NULL) {
        		perror("Error opening old file for writing");
       			exit(1);
    		}

    		// Open the new file for reading
   			FILE* newFile = fopen(G_sys_file_name, "rb");
    		if (newFile == NULL) {
        		perror("Error opening new file for reading");
        		fclose(oldFile);
        		exit(1);
    		}

    		// Copy the content from new file to old file
    		char buffer[1024]; // You can adjust the buffer size as needed
    		size_t bytesRead;
    		while ((bytesRead = fread(buffer, 1, sizeof(buffer), newFile)) > 0) {
        		fwrite(buffer, 1, bytesRead, oldFile);
    		}

    		// Close both files
    		fclose(oldFile);
    		fclose(newFile);
		}
		
	}
	
}
