#include "computational/solvers/iterative_solver.h"
#include <stdio.h>
#include <string.h>
#include "computational/GreensFunction.h"
#include <stdlib.h>
#include <lapacke.h>
#include <cblas.h>
#include <stdbool.h>
#include <complex.h>
#include <math.h>
#include "functions_DSGF.h"
#include "file_utils.h"

//#include "mkl.h" // if uncommented, a series of warning are shown when compiling.

void matrix_reshape(int inner_size, int outer_size, double complex matrix_2d_1[][3*outer_size], double complex matrix_4d_1[outer_size][outer_size][inner_size][inner_size]){

	for (int major_row = 0; major_row < outer_size; major_row++)
	{
		for (int major_y = 0; major_y < outer_size; major_y++)
		{ 
			for(int minor_x = 0; minor_x < inner_size; minor_x++)
			{
				for(int minor_y = 0; minor_y < inner_size; minor_y++)
				{
					int new_x = 3*major_row + minor_x;
					int new_y = 3*major_y + minor_y;
					matrix_2d_1[new_x][new_y]=matrix_4d_1[major_row][major_y][minor_x][minor_y]; 
				}
			}
		}
	}

}

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

	int ipack=0;

	for (int row = 0; row < 3; row++)
	{
		int new_x = (3*mm + row);
		for(int col = 0; col < 3; col++) // 3D coordinate positions 
		{
			int new_y = (3*jg_0 + col);
			int index = 9 * (new_x * tot_sub_vol + new_y) + 3 * row + col;
			A_iterative[ipack] = A_2d[row][col]; //A_2d[mm_2d][mm_2d]; //A[mm][mm][i_subG_0][j_subG_0];
			b_iterative[ipack] = G_sys_old[new_x][new_y]; //G_sys_old[mm][jg_0][i_subG_0][j_subG_0];
			//b_iterative[ipack] = upperTriangularMatrix[index];
			ipack++;
		}    
	}
}

void G_sys_new_populator(int tot_sub_vol, int mm, int jg_0, double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], double complex b_iterative[]){

	int gpack=0;     	

	for (int mm_sub = 0; mm_sub < 3; mm_sub++) // 3D coordinate positions
	{
		int mm_2d = (3*mm + mm_sub);
		for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
		{
			int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices

			G_sys_new[mm_2d][jg_0_2d] = b_iterative[gpack]; // stores G^1_11
			G_sys_new[jg_0_2d][mm_2d] = G_sys_new[mm_2d][jg_0_2d]; //reciprocity
			gpack++;
		}  
	}
}

void offdiagonal_solver(int tot_sub_vol, int mm, double k, double complex alpha_0[], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol]){
	/*
	double complex G_sys_new_transpose[3*tot_sub_vol][3*tot_sub_vol];
	
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{ 
		for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
		{
			int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
			for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //lower triangular matrix
			{
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
				{
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					G_sys_new_transpose[ig_0_2d][jg_0_2d] = G_sys_new[jg_0_2d][ig_0_2d]; //https://akkadia.org/drepper/cpumemory.pdf
				}
			}
		}
	}		
	*/
		
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
		if(ig_0 != mm)
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

/*
	int N = tot_sub_vol;
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) {
        if (ig_0 != mm) {
            for (int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) {
                int ig_0_2d = (3 * ig_0 + i_subG_0);

                for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) {
                    for (int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) {
                        int jg_0_2d = (3 * jg_0 + j_subG_0);

                        double complex G_sys_prod = 0.0;

                        for (int m_sub = 0; m_sub < 3; m_sub++) {
                            int mm_2d = (3 * mm + m_sub);
                            G_sys_prod += G_sys_old[ig_0_2d][mm_2d] * G_sys_new[mm_2d][jg_0_2d];
                        }

                        // Perform matrix multiplication using CBLAS
                        double complex alpha = pow(k, 2) * alpha_0[mm];
						double complex beta = 0;
                        cblas_zgemm(
                            CblasRowMajor, CblasNoTrans, CblasNoTrans,
                            3*N, 1, 3*N,  // Dimensions of matrices
                            &alpha,   // Scaling factor for G_sys_prod
                            &G_sys_old[ig_0_2d][jg_0_2d], N * 3,  // Matrix G_sys_old
                            &G_sys_prod, 1,  // Matrix G_sys_prod
							&beta,
                            &G_sys_new[ig_0_2d][jg_0_2d], N * 3  // Matrix G_sys_new
                        );
                    }
                }
            }
        }
    }
*/

}


void remaining_pertubations(int tot_sub_vol, int mm, double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], double complex A_2d[3][3]){

	double complex A_iterative[3*3];
	double complex b_iterative[3*3];
	for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) // Only loop through remaining perturbations
	{
		
		A_b_iterative_populator(tot_sub_vol, A_iterative, b_iterative, A_2d, G_sys_old, mm, jg_0);
		//A_b_iterative_populator(tot_sub_vol, A_iterative, b_iterative, A_2d, size, upperTriangularMatrix, mm, jg_0);

		// %%%%%%%%%%% G_new using Linear inversion using LAPACK %%%%%%%%%%%%%%%%%
		int info = LAPACKE_zgels(LAPACK_ROW_MAJOR,'N',3,3,3,A_iterative,3,b_iterative,3);
		G_sys_new_populator(tot_sub_vol, mm, jg_0, G_sys_new, b_iterative);

	} // end jg_0					
	
}

void core_solver(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[],double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol]){


	double complex A_2d[3][3];
	for (int mm = 0; mm < tot_sub_vol; mm++) //tot_sub_vol
	{
		//printf("%d - ",mm+1);

		double complex epsilon_s = (epsilon - epsilon_ref); // Scattering dielectric function

		A2d_solver(epsilon_s, mm, tot_sub_vol, delta_V_vector[mm], G_sys_old, A_2d, k);	
		remaining_pertubations(tot_sub_vol, mm, G_sys_old, G_sys_new, A_2d);
		offdiagonal_solver(tot_sub_vol, mm, k, alpha_0, G_sys_old, G_sys_new);

		memcpy(G_sys_old,G_sys_new,3*tot_sub_vol*3*tot_sub_vol*sizeof(double complex)); // Update G_old = G_new for next iteration.

	}//end mm loop 
}

void core_solver_store(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[],double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], char* G_sys_file_name){


	double complex A_2d[3][3];
	for (int mm = 0; mm < tot_sub_vol; mm++) //tot_sub_vol
	{
		//printf("%d - ",mm+1);

		double complex epsilon_s = (epsilon - epsilon_ref); // Scattering dielectric function

		A2d_solver(epsilon_s, mm, tot_sub_vol, delta_V_vector[mm], G_sys_old, A_2d, k);	
		remaining_pertubations(tot_sub_vol, mm, G_sys_old, G_sys_new, A_2d);
		offdiagonal_solver(tot_sub_vol, mm, k, alpha_0, G_sys_old, G_sys_new);

		memcpy(G_sys_old,G_sys_new,3*tot_sub_vol*3*tot_sub_vol*sizeof(double complex)); // Update G_old = G_new for next iteration.

		FILE* G_sys_file = fopen(G_sys_file_name, "w");
  		 // Write the array elements to the file
		for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //lower triangular
    	{
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions 
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //upper triangular matrix
				{
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
					{
						int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
						fprintf(G_sys_file, "%e + i %e , ", creal(G_sys_new[ig_0_2d][jg_0_2d]), cimag(G_sys_new[ig_0_2d][jg_0_2d]));
    				}
				}
				fprintf(G_sys_file, "\n");
			}
		}			
    	fclose(G_sys_file); // Close the file


	}//end mm loop 
}

void iterative_solver(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol],double k_0, double pi,double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3],char wave_type){

	double complex (*G_0)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*G_0)); 
	//get_G0_matrix(tot_sub_vol, G_0, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, wave_type);
	get_G0_matrix_memory(tot_sub_vol, G_0, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, wave_type);
			
	double complex (*G_sys_old)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys_old));
	matrix_reshape(3, tot_sub_vol, G_sys_old, G_0); // this reshapes 2 4D matrices to 2 2D matrices, where G0 and eyeA are 4D and G_sys_old and eyeA_2d are the respective 2D matrices
	free(G_0);

	//double complex (*G_sys_new)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys_new));
    
	//older version
	//core_solver(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_sys_new, G_sys_old);
	core_solver(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_sys, G_sys_old);
	
	free(G_sys_old);
	//memcpy(G_sys,G_sys_new,3*tot_sub_vol*3*tot_sub_vol*sizeof(double complex)); //Populate G_sys with G_new
	//free(G_sys_new);

}

void iterative_solver_store(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol],double k_0, double pi,double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3],char wave_type, char* G_sys_file_name){

	/*
	double complex (*G_0)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*G_0)); 
	//get_G0_matrix(tot_sub_vol, G_0, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, wave_type);
	get_G0_matrix_memory(tot_sub_vol, G_0, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, wave_type);
			
	double complex (*G_old)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_old));
	matrix_reshape(3, tot_sub_vol, G_old, G_0); // this reshapes 2 4D matrices to 2 2D matrices, where G0 and eyeA are 4D and G_sys_old and eyeA_2d are the respective 2D matrices
	free(G_0);
	*/

	double complex (*G_old)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_old));
	//double complex (*G_old)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*G_old)); 
	if (G_old == NULL){
			printf("Failure with memory in iterative solver.");
			exit(1);
	} 
	//get_G_old_matrix_memory(tot_sub_vol, G_old, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, wave_type);
	get_G_old_struct_matrix_memory(tot_sub_vol, G_old, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, wave_type);

	//double complex (*G_sys_new)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys_new));
    //printf("here\n");
	//older versions
	//core_solver(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_sys_new, G_sys_old);
	//core_solver(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_sys, G_old);
	
	//store G_old and G_sys
	//double complex (*G_new)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_new));
	core_solver_store(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_sys, G_old, G_sys_file_name);
	
	free(G_old);
	//memcpy(G_sys,G_sys_new,3*tot_sub_vol*3*tot_sub_vol*sizeof(double complex)); //Populate G_sys with G_new
	//free(G_sys_new);

}



void iterative_solver_with_data(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[],double k_0, double pi,double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3],char wave_type, char* G_old_file_name, char* G_sys_file_name){

	double complex (*G_old)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_old));
	//double complex (*G_old)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*G_old)); 
	if (G_old == NULL){
			printf("Failure with memory in iterative solver.");
			exit(1);
	} 
	//get_G_old_matrix_memory(tot_sub_vol, G_old, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, wave_type);
	get_G_old_struct_matrix_memory(tot_sub_vol, G_old, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, wave_type);
	
	// Write the array elements to the file
	write_bin(tot_sub_vol, G_old, G_old_file_name); //write_csv(tot_sub_vol, G_old, G_old_file_name);
	
	free(G_old);

	//older versions of core solver:
	//core_solver(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_sys_new, G_sys_old);
	//core_solver(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_sys, G_old);
	
	//core solver storing G_old and G_sys without memory improvement:
	//core_solver_store(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_sys, G_old, G_sys_file_name);
	
	double complex A_2d[3][3];
	double complex A_2d_test[3][3];
	for (int mm = 0; mm < tot_sub_vol; mm++) //tot_sub_vol
	{
		printf("%d - ",mm+1);

		double complex epsilon_s = (epsilon - epsilon_ref); // Scattering dielectric function

		double complex (*G_old)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_old));
		if (G_old == NULL){
			printf("Failure with memory in iterative solver.");
			exit(1);
		}
		// Read file into array elements
		read_bin(tot_sub_vol, G_old, G_old_file_name); //read_csv(tot_sub_vol, G_old, G_old_file_name);

		//This function generates Amm and uses G_old, but we do not need all the terms to generate Amm. 
		//A2d_solver(epsilon_s, mm, tot_sub_vol, delta_V_vector[mm], G_old, A_2d, k);
		
		//I copied the body of the A2d_solver function into the code as an attempt to use the minimum of terms needed (9). 
		double eyeA_2d[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};
		/*
		struct Eye3x3 eyeA_2d = {
        .xx = 1.0, .xy = 0.0, .xz = 0.0,
        .yx = 0.0, .yy = 1.0, .yz = 0.0,
        .zx = 0.0, .zy = 0.0, .zz = 1.0 };
		*/
		double complex G_oldValue; //struct Matrix3x3 G_old_mm; //
		double complex A_2d[3][3]; //struct Matrix3x3 A_2d_test; 
		FILE * G_old_import = fopen(G_old_file_name, "rb"); 
		if (G_old_import == NULL) {
    		perror("Error opening binary file");
    		exit(1); // Exit with an error code
		}
		// Set the file pointer to the start of the matrix data
    	//fseek(G_old_import, 0, SEEK_SET); // 0 bytes from the start of the file
		//size_t elements_read = fread(&G_old_mm, sizeof(struct Matrix3x3), 1, G_old_import); // Read the matrix data from the binary file into the struct
    	double complex cte = - pow(k, 2) * delta_V_vector[mm] * epsilon_s;
		//printf("%e + i %e ,%e + i %e ,%e + i %e \n %e + i %e ,%e + i %e ,%e + i %e \n %e + i %e ,%e + i %e , %e + i %e \n",creal(G_old_mm.xx),cimag(G_old_mm.xx),creal(G_old_mm.xy),cimag(G_old_mm.xy),creal(G_old_mm.xz),cimag(G_old_mm.xz),creal(G_old_mm.yx),cimag(G_old_mm.yx),creal(G_old_mm.yy),cimag(G_old_mm.yy),creal(G_old_mm.yz),cimag(G_old_mm.yz),creal(G_old_mm.zx),cimag(G_old_mm.zx),creal(G_old_mm.zy),cimag(G_old_mm.zy),creal(G_old_mm.zz),cimag(G_old_mm.zz));
		for (int row = 0; row < 3; row++)
		{
			for (int column = 0; column < 3; column++)
			{
				//int position = 9*mm+9*tot_sub_vol*mm+3*tot_sub_vol*row+column;
				int position = (9*tot_sub_vol+3)*mm+3*tot_sub_vol*row+column;
				printf("position =%d , ",position); 
				fseek(G_old_import, position * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
        		//fseek(G_old_import, 0, position); // 0 bytes from the start of the file
				size_t elements_read = fread(&G_oldValue, sizeof(double complex), 1, G_old_import); // Read the matrix data from the binary file into the struct
				printf("G_old[%d][%d] = %e + i %e , ",row, column,creal(G_oldValue),cimag(G_oldValue));
				A_2d[row][column] =  eyeA_2d[row][column] + cte* G_oldValue; //A_2d_test =  &eyeA_2d - pow(k, 2) * delta_V_vector[mm] * epsilon_s * &G_old_mm;
				printf("A_2d[%d][%d] = %e + i %e , ",row, column,creal(A_2d[row][column]),cimag(A_2d[row][column]));
			}
			printf("\n");
		}
    	fclose(G_old_import); // Close the file
		
		
		/*
		A_2d_test.xx = eyeA_2d.xx - cte* G_old_mm.xx;
		A_2d_test.xy = eyeA_2d.xy - cte* G_old_mm.xy;
		A_2d_test.xz = eyeA_2d.xz - cte* G_old_mm.xz;
		*/
		//printf("A2d: %e + i %e ,%e + i %e ,%e + i %e \n ",creal(A_2d_test.xx),cimag(A_2d_test.xx),creal(A_2d_test.xy),cimag(A_2d_test.xy),creal(A_2d_test.xz),cimag(A_2d_test.xz));
		
		/*
		if (elements_read != 1) {
        	perror("Error reading from binary file");
        	fclose(binFile);
        	return 1; // Exit with an error code
    	}
		*/
		//complex double value = realPart + imagPart * I;
		/*
		
		
		
		*/

	
	/*
	// Print values for verification
	 if (mm==0)
	 { 
		//for (int i = 0; i < 3; i++){for (int j = 0; j < 3; j++){printf("A_2d[%d,%d] = %e + %e i , ", i, j, creal(A_2d[i][j]), cimag(A_2d[i][j]));}}

		printf("\nPrint what was imported\n");
		for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //rows
    	{
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions 
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //columns
				{
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
					{
						int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
						printf("G_old[%d,%d] = %e + %e i , ", ig_0_2d, jg_0_2d, creal(G_old[ig_0_2d][jg_0_2d]), cimag(G_old[ig_0_2d][jg_0_2d]));
					}
				}
				printf("\n");
			}
		}	
	 }
	 */
		/*
		// Ideally, we would not need to store this 3Nx3N array.
		double complex (*G_sys_new)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys_new));
		if (G_sys_new == NULL){
			printf("Failure with memory in iterative solver.");
			exit(1);
		} 	
		*/
			
		//This function calculates the first G_new for when i=m. It uses G_old, but we do not need all the terms. 
		//I copied the body of the function into the code as an attempt to use the minimum of terms needed (9).
		//remaining_pertubations(tot_sub_vol, mm, G_old, G_sys_new, A_2d);
	
		double complex A_iterative[3*3];
		double complex b_iterative[3*3];
		#define MAX_LINE_LENGTH 1024
		for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) // Only loop through remaining perturbations
		{
			A_b_iterative_populator(tot_sub_vol, A_iterative, b_iterative, A_2d, G_old, mm, jg_0); // Populate Amm as A_iterative and G_sys_old as b_iterative
			int info = LAPACKE_zgels(LAPACK_ROW_MAJOR,'N',3,3,3,A_iterative,3,b_iterative,3); // G_new using Linear inversion using LAPACK
			
			//Ideally, instead of using this function, we would store b_iterative directly in a file, according to the term's position.
			//G_sys_new_populator(tot_sub_vol, mm, jg_0, G_sys_new, b_iterative);
			int gpack=0;     	
			FILE * G_new_export = fopen(G_sys_file_name, "wb"); 
			if (G_new_export == NULL) {
    			perror("Error opening binary file");
    			exit(1); // Exit with an error code
			}
			for (int mm_sub = 0; mm_sub < 3; mm_sub++) // 3D coordinate positions
			{
				int mm_2d = (3*mm + mm_sub);
				for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
				{
					int position_mj = 9*tot_sub_vol*mm+3*jg_0+3*tot_sub_vol*mm_sub+j_subG_0; //correct position  
					//printf("position_mj =%d , ",position_mj); 
					fseek(G_new_export, position_mj * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
					size_t elements_read_mj = fwrite(&b_iterative[gpack], sizeof(double complex), 1, G_new_export); // Read the matrix data from the binary file into the struct
					
					// reciprocity
					int position_jm = 9*tot_sub_vol*jg_0+3*mm+3*tot_sub_vol*j_subG_0+mm_sub; //correct position  
					//printf("position_jm =%d , ",position_jm); 
					fseek(G_new_export, position_jm * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
					size_t elements_read_jm = fwrite(&b_iterative[gpack], sizeof(double complex), 1, G_new_export); // Read the matrix data from the binary file into the struct
					
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					//G_sys_new[mm_2d][jg_0_2d] = b_iterative[gpack]; // stores G^1_11
					//G_sys_new[jg_0_2d][mm_2d] = G_sys_new[mm_2d][jg_0_2d]; //reciprocity
					printf("G_new[%d][%d] = %e + i %e , ",mm_2d, jg_0_2d,creal(b_iterative[gpack]),cimag(b_iterative[gpack]));
					
					gpack++;
				}  
				printf("\n");
			}
			fclose(G_new_export); // Close the file    	
			
		} // end jg_0

		//This function calculates the remaining G_new for when i!=m. It uses G_old, but we do not need all the terms. 
		//I copied the body of the function into the code as an attempt to use the minimum of terms needed. 
		//Ideally, we would import the three terms needed for the calculations and store the result in the G_sys file.
		//offdiagonal_solver(tot_sub_vol, mm, k, alpha_0, G_old, G_sys_new);
		
		/*
		
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
					//if(ig_0 != mm ) //condition without reciprocity
					if(ig_0 != mm && jg_0 != mm) //using reciprocity, some terms would be calculated twice. The updated condition, each term is calculated one
					{
						for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
						{
							int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
							double complex G_sys_prod = 0.;
							for(int m_sub = 0;  m_sub < 3;  m_sub++)//loop for matricial multiplication
							{
								int mm_2d = (3*mm + m_sub);
								G_sys_prod += G_old[ig_0_2d][mm_2d]*pow(k,2)*alpha_0[mm]*G_sys_new[mm_2d][jg_0_2d];
							}
							G_sys_new[ig_0_2d][jg_0_2d] = G_old[ig_0_2d][jg_0_2d] + G_sys_prod; //alpha_0[mm]*[ig_0_2d][jg_0_2d]
							G_sys_new[jg_0_2d][ig_0_2d] = G_sys_new[ig_0_2d][jg_0_2d]; // reciprocity
						} // i_subG_0    
					} //  if(ig_0 != mm) 	
				}// ig_0 	        				
			} // j_subG_0                 	
		} // jg_0 
		
		*/   

		// data version
		double complex G_new_ij_value, G_old_ij_value;
		double complex G_old_im_value,G_new_mj_value[3];
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
							//int position = 9*mm+9*tot_sub_vol*mm+3*tot_sub_vol*i_subG_0+j_subG_0; // used in Amm
							//int position = 9*tot_sub_vol*mm+3*jg_0+3*tot_sub_vol*i_subG_0+j_subG_0; //used in G_new_m,j  
							int position_old_ij = 9*tot_sub_vol*ig_0+3*tot_sub_vol*i_subG_0+3*jg_0+j_subG_0; //seems to be correct
							//int position_old_im = 9*tot_sub_vol*ig_0+3*tot_sub_vol*i_subG_0+3*mm; //seems to be correct
							int position_old_im = 9*tot_sub_vol*ig_0+3*tot_sub_vol*i_subG_0+3*mm+j_subG_0; //seems to be correct
							int position_new_mj = 9*tot_sub_vol*mm+3*jg_0+3*tot_sub_vol*i_subG_0; //seems to be correct
							printf("position_old_ij =%d , position_old_im = %d, position_new_mj = %d: ",position_old_ij, position_old_im, position_new_mj); 
							
							G_old_import = fopen(G_old_file_name, "rb"); 
							if (G_old_import == NULL) {
    							perror("Error opening binary file");
    							exit(1); // Exit with an error code
							}

							fseek(G_old_import, position_old_ij * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
        					//fseek(G_old_import, 0, position); // 0 bytes from the start of the file
							size_t elements_read = fread(&G_old_ij_value, sizeof(double complex), 1, G_old_import); // Read the matrix data from the binary file into the struct
							
							fseek(G_old_import, position_old_im * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
        					size_t SGF_old_read = fread(&G_old_im_value, sizeof(double complex), 1, G_old_import); // Read the matrix data from the binary file into the struct
							
							fclose(G_old_import);

							FILE * G_new_import = fopen(G_sys_file_name, "r+b");
							if (G_new_import == NULL) {
    						perror("Error opening binary file");
    						exit(1); // Exit with an error code
							}
							fseek(G_new_import, position_new_mj * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
        					size_t SGF_new_read = fread(&G_new_mj_value, sizeof(double complex), 3, G_new_import); // Read the matrix data from the binary file into the struct
							

							double complex G_sys_prod = 0.;
							for(int m_sub = 0;  m_sub < 3;  m_sub++)//loop for matricial multiplication
							{
								int mm_2d = (3*mm + m_sub);
								printf("G_old_im_value = %e + i %e, ",creal(G_old_im_value),cimag(G_old_im_value));
								printf("G_new_mj_value = %e + i %e; ",creal(G_new_mj_value[m_sub]),cimag(G_new_mj_value[m_sub]));
								G_sys_prod += G_old_im_value*pow(k,2)*alpha_0[mm]*G_new_mj_value[m_sub];
							}
							G_new_ij_value = G_old_ij_value + G_sys_prod; //alpha_0[mm]*[ig_0_2d][jg_0_2d]
							printf("G_sys_new[%d][%d][%d][%d]=%e+i%e\n",ig_0,jg_0,i_subG_0,j_subG_0,creal(G_new_ij_value),cimag(G_new_ij_value));
							fseek(G_new_import, position_old_ij * sizeof(double complex), SEEK_SET);
							size_t SGF_new_write = fwrite(&G_new_ij_value, sizeof(double complex), 1, G_new_import);
							
						} // i_subG_0    
					} //  if(ig_0 != mm) 	
				}// ig_0 	      
				printf("\n");  				
			} // j_subG_0                 	
		} // jg_0    
		
		//memcpy(G_old,G_sys_new,3*tot_sub_vol*3*tot_sub_vol*sizeof(double complex)); // Update G_old = G_new for next iteration.
		
	/*
		printf("\nPrint what will be exported\n");
		for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //rows
    	{
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions 
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //columns
				{
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
					{
						int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
						printf("G_old[%d,%d] = %e + %e i , ", ig_0_2d, jg_0_2d, creal(G_old[ig_0_2d][jg_0_2d]), cimag(G_old[ig_0_2d][jg_0_2d]));
					}
				}
				printf("\n");
			}
		}
	*/

		if (mm<tot_sub_vol-1) {
			// Update G_old 
			//write_csv(tot_sub_vol, G_old, G_old_file_name); // 19s for 100 frequencies and 16 subvolumes
			//write_bin(tot_sub_vol, G_old, G_old_file_name); // 12s for 100 frequencies and 16 subvolumes
			
			/*
			//another strategy to update G_old
			char update_G_old[260];
			//sprintf(update_G_old, "cp ./%s/G_sys.csv  ./%s/G_old.csv\n",results_folder,results_folder);
			sprintf(update_G_old, "cp ./%s ./%s",G_sys_file_name,G_old_file_name);
			system(update_G_old);
			*/	

			/*
			FILE *oldFile = fopen(G_old_file_name, "wb");  // Open old binary file for writing (this clears existing data)
			if (oldFile == NULL) {
    			perror("Error opening binary file");
    			exit(1); // Exit with an error code
			}
			fclose(oldFile);
			
			*/

			/*
			FILE *newFile = fopen(G_sys_file_name, "rb");  // Open new binary file for reading
			if (newFile == NULL) {
    			perror("Error opening binary file");
    			exit(1); // Exit with an error code
			}

			FILE *oldFile = fopen(G_old_file_name, "wb");  // Open old binary file for writing (this clears existing data)
			if (oldFile == NULL) {
    			perror("Error opening binary file");
    			exit(1); // Exit with an error code
			}
			size_t dataSize = 9 * tot_sub_vol * tot_sub_vol * sizeof(double complex);
			double complex *dataBuffer = (double complex *)malloc(dataSize); //
			if (dataBuffer == NULL) {
   			// Handle memory allocation error
    		fclose(oldFile);
    		fclose(newFile);
    		exit(1); // Exit with an error code
			}			
			fread(dataBuffer, 1, dataSize, newFile);
			fwrite(dataBuffer, 1, dataSize, oldFile);
			
			fclose(oldFile);
			fclose(newFile);
			free(dataBuffer);
			*/

			//CopyFile(G_sys_file_name, G_old_file_name);
			//system("cp x.exe y.bas -f");

			FILE *oldFile = fopen(G_old_file_name, "wb");  // Open old binary file for writing (this clears existing data)
			if (oldFile == NULL) {
    			perror("Error opening binary file");
    			exit(1); // Exit with an error code
			}
			fclose(oldFile);

			char copy_G_new_to_G_old[260];
			sprintf(copy_G_new_to_G_old, "cp ./%s ./%s\n",G_sys_file_name,G_old_file_name);
			system(copy_G_new_to_G_old);

			FILE *newFile = fopen(G_sys_file_name, "wb");  // Open the new binary file in write mode (this clears existing data)
			if (newFile == NULL) {
    			perror("Error opening binary file");
    			exit(1); // Exit with an error code
			}
			fclose(newFile);	


			/*
			
			FILE *newFile = fopen(G_sys_file_name, "rb");  // Open new binary file for reading
			if (newFile == NULL) {
    			perror("Error opening binary file");
    			exit(1); // Exit with an error code
			}

			FILE *oldFile = fopen(G_old_file_name, "wb");  // Open old binary file for writing (this clears existing data)
			if (oldFile == NULL) {
    			perror("Error opening binary file");
    			exit(1); // Exit with an error code
			}
			char buffer[1024]; // A buffer to hold data
			size_t bytesRead;
			int c;
		
			//while ((bytesRead = fread(buffer, 1, sizeof(buffer), newFile)) > 0) {
			//	//size_t SGF_new_write = fwrite(&G_new_ij_value, sizeof(double complex), 1, G_new_import);
    		//	fwrite(buffer, 1, bytesRead, oldFile); // Write data to the old binary file
			//}
			
			while ((c = fgetc(newFile)) != EOF) fputc(c,oldFile);
			
			fclose(oldFile);
			fclose(newFile);
			*/

			/*	
			FILE *newFile = fopen(G_sys_file_name, "wb");  // Open the new binary file in write mode (this clears existing data)
			fclose(newFile);
			*/			
			
			}

		else {
			// write G_new 
			//write_csv(tot_sub_vol, G_sys_new,G_sys_file_name); 
			
			//write_bin(tot_sub_vol, G_sys_new,G_sys_file_name);
		}

		free(G_old);
		//free(G_sys_new);
		
	}
	
}


// *****************************************************
//   TRIANGULAR SOLUTION TO REDUCE MEMORY REQUIREMENTS 
// *****************************************************

void A2d_solver_memory(double complex epsilon, int mm, int tot_sub_vol, double delta_V, int size, double complex G_0_TriangularMatrix[size], double complex A_2d[3][3], double k){
	printf("\nA2d_solver_memory: ");
	double eyeA_2d[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};
	int index =mm*9;
	for(int minor_x = 0; minor_x < 3; minor_x++){
		int new_x = 3*mm + minor_x;
		for (int minor_y = 0; minor_y < 3; minor_y++){
			int new_y = 3*mm + minor_y;
            //int index = 9 * (new_x * tot_sub_vol + new_y) + 3 * minor_x + minor_y; // Calculate the index in the 1D array
			printf("%d -",index );
			A_2d[minor_x][minor_y] = eyeA_2d[minor_x][minor_y] - pow(k,2)*delta_V*epsilon*G_0_TriangularMatrix[index];
			index++;
		}
	}
}

void remaining_pertubations_memory(int tot_sub_vol, int mm, int size, double complex G_0_TriangularMatrix[size], double complex G_sys_TriangularMatrix[size], double complex A_2d[3][3]){

	double complex A_iterative[3*3];
	double complex b_iterative[3*3];
	for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) // Only loop through remaining perturbations
	{
		
		//A_b_iterative_populator(tot_sub_vol, A_iterative, b_iterative, A_2d, G_sys_old, mm, jg_0);
		A_b_iterative_populator_memory(tot_sub_vol, A_iterative, b_iterative, A_2d, size, G_0_TriangularMatrix, mm, jg_0);

		// %%%%%%%%%%% G_new using Linear inversion using LAPACK %%%%%%%%%%%%%%%%%
		int info = LAPACKE_zgels(LAPACK_ROW_MAJOR,'N',3,3,3,A_iterative,3,b_iterative,3);
		G_sys_new_populator_memory(tot_sub_vol, mm, jg_0, size, G_sys_TriangularMatrix, b_iterative);

	} // end jg_0					
	
}

void A_b_iterative_populator_memory(int tot_sub_vol, double complex *A_iterative, double complex *b_iterative, double complex A_2d[3][3], int size, double complex G_0_TriangularMatrix[size],int mm, int jg_0){

	int ipack=0;

	for (int row = 0; row < 3; row++)
	{
		int new_x = (3*mm + row);
		for(int col = 0; col < 3; col++) // 3D coordinate positions 
		{
			int new_y = (3*jg_0 + col);
			// INDEX!!!!!
			int index = 9 * (new_x * tot_sub_vol + new_y) + 3 * row + col;
			A_iterative[ipack] = A_2d[row][col]; //A_2d[mm_2d][mm_2d]; //A[mm][mm][i_subG_0][j_subG_0];
			b_iterative[ipack] = G_0_TriangularMatrix[index]; 
			ipack++;
		}    
	}
}

void G_sys_new_populator_memory(int tot_sub_vol, int mm, int jg_0, int size, double complex G_sys_TriangularMatrix[size], double complex b_iterative[]){


	int gpack=0;     	

	for (int mm_sub = 0; mm_sub < 3; mm_sub++) // 3D coordinate positions
	{
		int mm_2d = (3*mm + mm_sub);
		for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
		{
			int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
			// INDEX!!!!!!
			int index = index = 9 * (mm_2d * tot_sub_vol + jg_0_2d) + 3 * mm_sub + j_subG_0;
			//int index_transpose = 9 * ( jg_0_2d* tot_sub_vol + mm_2d) + 3 * j_subG_0 + mm_sub;  
			G_sys_TriangularMatrix[index] = b_iterative[gpack]; // stores G^1_11
			//G_sys_TriangularMatrix[jg_0_2d][mm_2d] = G_sys_TriangularMatrix[mm_2d][jg_0_2d]; //reciprocity
			gpack++;
		}  
	}
}

void offdiagonal_solver_memory(int tot_sub_vol, int mm, double k, double complex alpha_0[], int size, double complex G_0_TriangularMatrix[size], double complex G_sys_TriangularMatrix[size]){
		
	// Next, solve all systems of equations for ii not equal to mm		
	//for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //complete system
	for (int ig_0 = 0; ig_0 < 3*tot_sub_vol; ig_0++) //lower triangular
	{ 
		if(ig_0 != mm)
		{
			for (int jg_0 = ig_0; jg_0 < 3* tot_sub_vol; jg_0++) //upper triangular matrix
			{
				int submatrix_i = ig_0 / 3; // Row index of the 3x3 submatrix
            	int submatrix_j = jg_0 / 3; // Column index of the 3x3 submatrix

            	int row_in_submatrix = ig_0 % 3; // Row index within the submatrix
            	int col_in_submatrix = jg_0 % 3; 
				// INDEX!!!!!!
				int index = 9 * (submatrix_i * tot_sub_vol + submatrix_j) + 3 * row_in_submatrix + col_in_submatrix;

				double complex G_sys_prod = 0.;
				for(int m_sub = 0;  m_sub < 3;  m_sub++)//loop for matricial multiplication
				{
					int mm_2d = (3*mm + m_sub);
					int  index_im = 9 * (submatrix_i * tot_sub_vol + mm) + 3 * row_in_submatrix + m_sub;
					int  index_mj = 9 * (mm * tot_sub_vol + submatrix_j) + 3 * m_sub + col_in_submatrix;
					//G_sys_prod += G_sys_old[ig_0_2d][mm_2d]*pow(k,2)*alpha_0[mm]*G_sys_new[mm_2d][jg_0_2d];
					G_sys_prod += G_0_TriangularMatrix[index_im]*pow(k,2)*alpha_0[mm]*G_sys_TriangularMatrix[index_mj];
				}
				//int index = 9 * (submatrix_i * tot_sub_vol + submatrix_j) + 3 * row_in_submatrix + col_in_submatrix;
				//G_sys_new[ig_0_2d][jg_0_2d] = upperTriangularMatrix[index] + G_sys_prod; //alpha_0[mm]*[ig_0_2d][jg_0_2d]
				G_sys_TriangularMatrix[index] = G_0_TriangularMatrix[index] + G_sys_prod; //alpha_0[mm]*[ig_0_2d][jg_0_2d]
				//G_sys_new[jg_0_2d][ig_0_2d] = G_sys_TriangularMatrix[ig_0_2d][jg_0_2d]; // reciprocity   
			} // jg_0    	        				
		} // if(ig_0 != mm)              	
	} // ig_0    
}

void core_solver_memory(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[],int size, double complex G_new_TriangularMatrix[size], double complex G_old_TriangularMatrix[size]){

	double complex A_2d[3][3];
	for (int mm = 0; mm < tot_sub_vol; mm++) //tot_sub_vol
	{
		//printf("%d - ",mm+1);

		double complex epsilon_s = (epsilon - epsilon_ref); // Scattering dielectric function

		A2d_solver_memory(epsilon_s, mm, tot_sub_vol, delta_V_vector[mm], size, G_old_TriangularMatrix, A_2d, k);	//int index, double complex upperTriangularMatrix[index],
		remaining_pertubations_memory(tot_sub_vol, mm, size, G_old_TriangularMatrix, G_new_TriangularMatrix, A_2d);
		offdiagonal_solver_memory(tot_sub_vol, mm, k, alpha_0, size, G_old_TriangularMatrix, G_new_TriangularMatrix);

		memcpy(G_old_TriangularMatrix,G_new_TriangularMatrix,size*sizeof(double complex)); // Update G_old = G_new for next iteration.
	}//end mm loop 
}

void iterative_solver_memory(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[], int size, double complex G_sys_TriangularMatrix[size], double k_0, double pi,double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3],char wave_type){

	double complex* G_old_TriangularMatrix = (double complex*)malloc(size * sizeof(double complex)); // Allocate memory for the 1D array
	get_G0_triangular_matrix(tot_sub_vol, size, G_old_TriangularMatrix, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, wave_type); //
	
	//core_solver_memory(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, size, G_sys_TriangularMatrix, G_0_TriangularMatrix); //
	double complex A_2d[3][3];
	for (int mm = 0; mm < tot_sub_vol; mm++) //tot_sub_vol
	{
		printf("mm=%d - ",mm);

		double complex epsilon_s = (epsilon - epsilon_ref); // Scattering dielectric function

		//A2d_solver_memory(epsilon_s, mm, tot_sub_vol, delta_V_vector[mm], size, G_old_TriangularMatrix, A_2d, k);	//int index, double complex upperTriangularMatrix[index],
		printf("\nCompute Amm - ");
		double eyeA_2d[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};
		int index;
		if (mm==0) {index = 0;}
		else {index =mm*3*3*tot_sub_vol-pow(3,mm);}//{index =3*mm*(3*tot_sub_vol-mm);}
		//account for the -1 index at each row
		//printf("\n");
		for(int minor_x = 0; minor_x < 3; minor_x++)
		{
			int new_x = 3*mm + minor_x;
			//for (int minor_y = 0; minor_y < 3; minor_y++)
			for (int minor_y = minor_x; minor_y < 3; minor_y++)
			{
				int new_y = 3*mm + minor_y;
           	 	//index = 3*mm + 3 * minor_x + minor_y; // Calculate the index in the 1D array
				//index = (tot_sub_vol * minor_x - (minor_x * (minor_x + 1)) / 2) + (minor_y - minor_x);
				printf("[%d] ",index);
				//printf(" G_0[%d]= (%e + i%e) ",index, creal(G_old_TriangularMatrix[index]),cimag(G_old_TriangularMatrix[index]));
				A_2d[minor_x][minor_y] = eyeA_2d[minor_x][minor_y] - pow(k,2)*delta_V_vector[mm]*epsilon_s*G_old_TriangularMatrix[index];
				A_2d[minor_y][minor_x] = A_2d[minor_x][minor_y] ;
				//printf(" A_2d[%d,%d]= (%e + i%e) ",minor_x,minor_y, creal(A_2d[minor_x][minor_y]),cimag(A_2d[minor_x][minor_y]));
				index++;
			}
			index = index + 3*(tot_sub_vol-1) -3*mm;
			//printf("\n");
		}
		printf("\n");
				
		//remaining_pertubations_memory(tot_sub_vol, mm, size, G_old_TriangularMatrix, G_sys_TriangularMatrix, A_2d);
		printf("Pertubations for i=m \n");
		double complex A_iterative[3*3];
		double complex b_iterative[3*3];
		int index_pop=0;
		for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) // Only loop through remaining perturbations
		{
			printf("jg_0 = %d:",jg_0);
			if (mm==0) {index = 0;}
			else {index =mm*3*3*tot_sub_vol-pow(3,mm);}// {index =3*mm*(3*tot_sub_vol-mm);} //{index =mm*3*3*tot_sub_vol-pow(3,mm);}
			
			//A_b_iterative_populator(tot_sub_vol, A_iterative, b_iterative, A_2d, G_sys_old, mm, jg_0);
			//A_b_iterative_populator_memory(tot_sub_vol, A_iterative, b_iterative, A_2d, size, G_old_TriangularMatrix, mm, jg_0);
			int ipack=0;

			for (int row = 0; row < 3; row++)
			{
				//int new_x = (3*mm + row);
				//for(int col = 0; col < 3; col++) // 3D coordinate positions 
				for(int col = row; col < 3; col++) // 3D coordinate positions 
				{
					//int new_y = (3*jg_0 + col);
					if (row==col)
					{
						printf("[%d] ",index); //gpack
						A_iterative[ipack] = A_2d[row][col]; //A_2d[mm_2d][mm_2d]; //A[mm][mm][i_subG_0][j_subG_0];
						b_iterative[ipack] = G_old_TriangularMatrix[index]; 
						ipack++;
					}
					else
					{
						printf("[%d] ",index); //gpack
						A_iterative[ipack] = A_2d[row][col]; //A_2d[mm_2d][mm_2d]; //A[mm][mm][i_subG_0][j_subG_0];
						b_iterative[ipack] = G_old_TriangularMatrix[index]; 
						ipack++;
						printf("[%d] ",index); //gpack
						A_iterative[ipack] = A_2d[col][row]; //A_2d[mm_2d][mm_2d]; //A[mm][mm][i_subG_0][j_subG_0];
						b_iterative[ipack] = G_old_TriangularMatrix[index]; 
						ipack++;
					}
					index++;
				}
				index = index + 3*(tot_sub_vol-1) -3*mm; 
			}
			
			printf("\n");	
			// %%%%%%%%%%% G_new using Linear inversion using LAPACK %%%%%%%%%%%%%%%%%
			int info = LAPACKE_zgels(LAPACK_ROW_MAJOR,'N',3,3,3,A_iterative,3,b_iterative,3);
			
			//G_sys_new_populator_memory(tot_sub_vol, mm, jg_0, size, G_sys_TriangularMatrix, b_iterative);
			int gpack=0;     	
			printf("Populate G_sys mm=%d and jg_0=%d:\n",mm,jg_0);
			//if (jg_0==0) {index_pop = 0;}
			//else{};
			//else if (mm!=0 && jg_0==0){index =mm*3*3*tot_sub_vol-pow(3,mm);}
			//{index =mm*3*3*tot_sub_vol-pow(3,mm);} //{index = index + 3*(tot_sub_vol-1) -3*mm;}
			for (int mm_sub = 0; mm_sub < 3; mm_sub++) // 3D coordinate positions
			{
				if (mm!=jg_0)
				{
					for (int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
					{
						G_sys_TriangularMatrix[index_pop] = b_iterative[gpack]; 
						gpack++;
						printf("[%d] ",index_pop); 
						index_pop++;		
					}  
				}
				else
				{
					for (int j_subG_0 = mm_sub; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
					{
						if (j_subG_0==mm_sub)
						{
							G_sys_TriangularMatrix[index_pop] = b_iterative[gpack]; 
							gpack++;
							printf("[%d] ",index_pop); 
							index_pop++;
						}
						else{
							G_sys_TriangularMatrix[index_pop] = b_iterative[gpack]; 
							gpack++;
							gpack++;
							printf("[%d] ",index_pop); 
							index_pop++;
						}
					}
				}
			printf("\n");	 
			}
			//printf("%d: ",jg_0);
		} // end jg_0	
				
		//printf("\n");		
		//offdiagonal_solver_memory(tot_sub_vol, mm, k, alpha_0, size, G_old_TriangularMatrix, G_sys_TriangularMatrix);
		printf("Run pertubations: \n"); // Next, solve all systems of equations for ii not equal to mm		
		//for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //complete system

		//if (mm==0) {index = 0;}
		//else {index =mm*3*3*tot_sub_vol-pow(3,mm);}
		index = 0 ;
		int gpack=0;
		int decay =0;
		int old_decay = 0;
		for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //lower triangular
		//for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //lower triangular
		{ 
			
			if (ig_0 == mm)
			{
				printf("i==m: \n");
				index_pop= 9*mm;
				for (int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
				{
					for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++) //upper triangular matrix
					{
						if (ig_0!=jg_0)
						{
							for (int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions
							{
								G_sys_TriangularMatrix[index] = G_sys_TriangularMatrix[index_pop] ;
								printf("(%d = %d) ",index, index_pop);
								index_pop++;
								index++;
								gpack++;
							}	
						}
						else //(ig_0==jg_0)
						{
							for(int j_subG_0 = i_subG_0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
							{
								G_sys_TriangularMatrix[index] = G_sys_TriangularMatrix[index_pop] ;
								printf("(%d=%d) ",index,index_pop);
								index_pop++;
								index++;
								gpack++;
							}	
						}
					}
					printf("\n");  	
				}
			}
			
			if(ig_0 != mm)
			{ 
				printf("i!=m: \n");
				for (int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
				{
					for (int jg_0 = ig_0; jg_0 < tot_sub_vol; jg_0++) //upper triangular matrix
					{
						if (ig_0!=jg_0)
						{
							printf("ig_0!=jg_0\n");
							for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
							{
								double complex G_sys_prod = 0.;
								printf("\n{%d} ",index); 
								for(int m_sub = 0;  m_sub < 3;  m_sub++)//loop for matricial multiplication
								{
									printf("mm= %d, ig_0=%d, jg_0=%d, i_subG_0=%d, j_subG_0=%d, m_sub=%d ; ",mm,ig_0, jg_0,i_subG_0,j_subG_0,m_sub );
									int index_im = 3*tot_sub_vol*mm*jg_0 + 3*(jg_0) + 3*tot_sub_vol*i_subG_0 + m_sub -decay; //int index_mj = 3*mm + 3*ig_0 + 3*jg_0 + m_sub ;
									int index_mj = 3*tot_sub_vol*mm*jg_0 + 3*(jg_0) + 3*tot_sub_vol*i_subG_0 + m_sub -decay; //3*mm + 3*ig_0 + 3*jg_0 + m_sub ;
									printf("index_im=%d, ",index_im);
									printf("index_mj=%d ;\n",index_mj);
									//int mm_2d = (3*mm + m_sub);
									//G_sys_prod += G_sys_old[ig_0_2d][mm_2d]*pow(k,2)*alpha_0[mm]*G_sys_new[mm_2d][jg_0_2d];
									//G_sys_prod += G_old_TriangularMatrix[index]*pow(k,2)*alpha_0[mm]*G_sys_TriangularMatrix[index];
									G_sys_prod += G_old_TriangularMatrix[index_im]*pow(k,2)*alpha_0[mm]*G_sys_TriangularMatrix[index_mj];
									
								}
								//int index = 9 * (submatrix_i * tot_sub_vol + submatrix_j) + 3 * row_in_submatrix + col_in_submatrix;
								//G_sys_new[ig_0_2d][jg_0_2d] = upperTriangularMatrix[index] + G_sys_prod; //alpha_0[mm]*[ig_0_2d][jg_0_2d]
								G_sys_TriangularMatrix[index] = G_old_TriangularMatrix[index] + G_sys_prod; //alpha_0[mm]*[ig_0_2d][jg_0_2d]
								//G_sys_new[jg_0_2d][ig_0_2d] = G_sys_TriangularMatrix[ig_0_2d][jg_0_2d]; // reciprocity   
								index++;
								gpack++;
							}	
						}
						else if (ig_0==jg_0)
						{
							printf("ig_0==jg_0\n");
							for(int j_subG_0 = i_subG_0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
							{
								double complex G_sys_prod = 0.;
								printf("\n[%d] \n",index); 
								for(int m_sub = 0;  m_sub < 3;  m_sub++)//loop for matricial multiplication
								{
									printf("mm= %d, ig_0=%d, jg_0=%d, i_subG_0=%d, j_subG_0=%d, m_sub=%d; ",mm,ig_0, jg_0,i_subG_0,j_subG_0,m_sub );
									int index_im = 3*tot_sub_vol*mm*jg_0 + 3*(jg_0) + 3*tot_sub_vol*i_subG_0 + m_sub -decay; //    int index_mj = 3*mm + 3*ig_0 + 3*jg_0 + m_sub ;
									int index_mj = 3*tot_sub_vol*mm*jg_0 + 3*(jg_0) + 3*tot_sub_vol*i_subG_0 + m_sub -decay; //3*mm + 3*ig_0 + 3*jg_0 + m_sub ;
									printf("index_im=%d, ",index_im);
									printf("index_mj=%d ;\n",index_mj);
									//int mm_2d = (3*mm + m_sub);
									//G_sys_prod += G_sys_old[ig_0_2d][mm_2d]*pow(k,2)*alpha_0[mm]*G_sys_new[mm_2d][jg_0_2d];
									//G_sys_prod += G_old_TriangularMatrix[index]*pow(k,2)*alpha_0[mm]*G_sys_TriangularMatrix[index];
									G_sys_prod += G_old_TriangularMatrix[index_im]*pow(k,2)*alpha_0[mm]*G_sys_TriangularMatrix[index_mj];
									
								}
								//int index = 9 * (submatrix_i * tot_sub_vol + submatrix_j) + 3 * row_in_submatrix + col_in_submatrix;
								//G_sys_new[ig_0_2d][jg_0_2d] = upperTriangularMatrix[index] + G_sys_prod; //alpha_0[mm]*[ig_0_2d][jg_0_2d]
								G_sys_TriangularMatrix[index] = G_old_TriangularMatrix[index] + G_sys_prod; //alpha_0[mm]*[ig_0_2d][jg_0_2d]
								//G_sys_new[jg_0_2d][ig_0_2d] = G_sys_TriangularMatrix[ig_0_2d][jg_0_2d]; // reciprocity   
								index++;
								gpack++;
								
								
							}	
						}
						
						printf("jg_0=%d, i_subG_0= %d, decay=%d",jg_0,i_subG_0,decay);
						//old_decay= decay;
						//decay = old_decay + decay + i_subG_0*jg_0; //jg_0;
						decay = decay + i_subG_0+jg_0;
					} // for jg_0
					printf("\n");  
				} // i_subG_0 
				//printf("\n");   	        				
			} // if(ig_0 != mm) 
			//index = index + 3*(tot_sub_vol-1) -3*mm;             	
		} // ig_0    

		printf("terms calculated: %d. Update G_old = G_new\n \n", gpack-1); 

		memcpy(G_old_TriangularMatrix,G_sys_TriangularMatrix,size*sizeof(double complex)); // Update G_old = G_new for next iteration.
	
	}//end mm loop 
		
	free(G_old_TriangularMatrix);

}


