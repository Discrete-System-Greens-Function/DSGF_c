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

void iterative_solver_with_data(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[],double k_0, double pi,double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3],char wave_type, char* G_old_file_name, char* G_sys_file_name){

	double complex (*G_old)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_old));
	if (G_old == NULL){
			printf("Failure with memory in iterative solver.");
			exit(1);
	} 
	//get_G_old_matrix_memory(tot_sub_vol, G_old, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, wave_type);
	get_G_old_struct_matrix_memory(tot_sub_vol, G_old, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, wave_type);
	write_bin(tot_sub_vol, G_old, G_old_file_name); // Write the array elements to the file	
	free(G_old);
	
	double complex A_2d[3][3];
	for (int mm = 0; mm < tot_sub_vol; mm++) //tot_sub_vol
	{
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
		
		//This part calculates the first G_new for when i=m. It uses G_old, but we do not need all the terms. 
		double complex A_iterative[3*3];
		double complex b_iterative[3*3];
		#define MAX_LINE_LENGTH 1024 
		for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) // Only loop through remaining perturbations
		{
			int ipack=0;
			G_old_import = fopen(G_old_file_name, "rb");
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
				
					A_iterative[ipack] = A_2d[row][col]; //A_2d[mm_2d][mm_2d]; //A[mm][mm][i_subG_0][j_subG_0];
					b_iterative[ipack] = G_oldValue;
					ipack++;
				}    
			}
			fclose(G_old_import);
			printf("\n");
			int info = LAPACKE_zgels(LAPACK_ROW_MAJOR,'N',3,3,3,A_iterative,3,b_iterative,3); // G_new using Linear inversion using LAPACK
			
			//Ideally, instead of using this function, we would store b_iterative directly in a file, according to the term's position.
			//G_sys_new_populator(tot_sub_vol, mm, jg_0, G_sys_new, b_iterative);
			int gpack=0;     	
			
			FILE * G_new_export = fopen(G_sys_file_name, "w+");
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
					printf("position_mj =%d , ",position_mj);
                    rewind(G_new_export); 
                	fseek(G_new_export, position_mj * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
					size_t elements_read_mj = fwrite(&b_iterative[gpack], sizeof(double complex), 1, G_new_export); // Read the matrix data from the binary file into the struct
					if (elements_read_mj != 1) {
            		    printf("Error reading data at position %d\n", position_mj); // Handle error or add debugging information
        			}
					// reciprocity
					int position_jm = 9*tot_sub_vol*jg_0+3*mm+3*tot_sub_vol*j_subG_0+mm_sub; //correct position  
					printf("position_jm =%d , ",position_jm); 
					rewind(G_new_export);
					fseek(G_new_export, position_jm * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
					size_t elements_read_jm = fwrite(&b_iterative[gpack], sizeof(double complex), 1, G_new_export); // Read the matrix data from the binary file into the struct
					if (elements_read_jm != 1) {
            		    printf("Error reading data at position %d\n", position_jm); // Handle error or add debugging information
        			}
					int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
					printf("G_new[%d][%d] = %e + i %e , ",mm_2d, jg_0_2d,creal(b_iterative[gpack]),cimag(b_iterative[gpack]));
					gpack++;
				}  
				printf("\n");
			}
			fclose(G_new_export); // Close the file    	
			
		} // end jg_0
		double complex (*G_new)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_new));
	    if (G_new == NULL){
			printf("Failure with memory in iterative solver.");
			exit(1);
	    } 
	    read_bin(tot_sub_vol, G_new, G_sys_file_name);
	
		printf("\nPrint G_new so far\n");
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
						printf("G_new[%d,%d] = %e + %e i , ", ig_0_2d, jg_0_2d, creal(G_new[ig_0_2d][jg_0_2d]), cimag(G_new[ig_0_2d][jg_0_2d]));
					}
				}
				printf("\n");
			}
		}
		printf("\n");
		free(G_new);

		//This function calculates the remaining G_new for when i!=m. It uses G_old, but we do not need all the terms.
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

							FILE * G_new_import = fopen(G_sys_file_name, "rb");
							if (G_new_import == NULL) {
    						perror("Error opening binary file");
    						exit(1); // Exit with an error code
							}
							fseek(G_new_import, position_new_mj * sizeof(double complex), SEEK_SET); // Set the file position to the specified position
        					size_t SGF_new_read = fread(&G_new_mj_value, sizeof(double complex), 3, G_new_import); // Read the matrix data from the binary file into the struct
							fclose(G_new_import);

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
							G_new_import = fopen(G_sys_file_name, "r+b");
							fseek(G_new_import, position_old_ij * sizeof(double complex), SEEK_SET);
							size_t SGF_new_write = fwrite(&G_new_ij_value, sizeof(double complex), 1, G_new_import);
							fclose(G_new_import);
						} // i_subG_0    
					} //  if(ig_0 != mm) 	
				}// ig_0 	      
				printf("\n");  				
			} // j_subG_0                 	
		} // jg_0    
		
		//memcpy(G_old,G_sys_new,3*tot_sub_vol*3*tot_sub_vol*sizeof(double complex)); // Update G_old = G_new for next iteration.
		
	double complex (*G_new_2)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_new_2));
	//double complex (*G_old)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*G_old)); 
	if (G_new_2 == NULL){
			printf("Failure with memory in iterative solver.");
			exit(1);
	} 
	read_bin(tot_sub_vol, G_new_2, G_sys_file_name);
	
		printf("\nPrint G_new complete\n");
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
						printf("G_new[%d,%d] = %e + %e i , ", ig_0_2d, jg_0_2d, creal(G_new_2[ig_0_2d][jg_0_2d]), cimag(G_new_2[ig_0_2d][jg_0_2d]));
					}
				}
				printf("\n");
			}
		}
		printf("\n");
		free(G_new_2);

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

		//free(G_old);
		//free(G_sys_new);
		
	}
	
}