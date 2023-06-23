#include "computational/solvers/iterative_solver.h"
#include <stdio.h>
#include <string.h>
#include "computational/GreensFunction.h"
#include <stdlib.h>
#include <lapacke.h>
#include <cblas.h>
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
		}
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
			A_iterative[ipack] = A_2d[row][col]; //A_2d[mm_2d][mm_2d]; //A[mm][mm][i_subG_0][j_subG_0];
			b_iterative[ipack] = G_sys_old[new_x][new_y]; //G_sys_old[mm][jg_0][i_subG_0][j_subG_0];
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

			gpack++;
		}  
	}
}

void extractSubmatrix(
    double complex *largerMatrix, int largerRows, int largerCols,
    double complex *smallerMatrix, int smallerRows, int smallerCols,
    int startRow, int startCol
) {
    for (int i = 0; i < smallerRows; i++) {
        for (int j = 0; j < smallerCols; j++) {
            int largerIndex = (startRow + i) * largerCols + (startCol + j);
            int smallerIndex = i * smallerCols + j;
            smallerMatrix[smallerIndex] = largerMatrix[largerIndex];
        }
    }
}

void insertSubmatrix(
    double complex* largerMatrix, int largerRows, int largerCols,
    double complex* smallerMatrix, int smallerRows, int smallerCols,
    int startRow, int startCol
) {
    for (int i = 0; i < smallerRows; i++) {
        for (int j = 0; j < smallerCols; j++) {
            int largerIndex = (startRow + i) * largerCols + (startCol + j);
            int smallerIndex = i * smallerCols + j;
            largerMatrix[largerIndex] = smallerMatrix[smallerIndex];
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
	/*
	int smallerRows=3;
	int smallerCols=3;
	
	double complex *G_old_33 = malloc(smallerRows * smallerCols * sizeof(double complex)); //double complex (*G_old_33)[3] = calloc(3, sizeof(*G_old_33));
	double complex *G_new_33 = malloc(smallerRows * smallerCols * sizeof(double complex)); //double complex (*G_new_33)[3] = calloc(3, sizeof(*G_new_33));

	double complex alpha_parameter= pow(k,2)*alpha_0[mm] + 0.0 * I;  // Scaling factor for A*B
	double complex beta_parameter = 1.0 + 0.0 * I;   // Scaling factor for C, beta = 1.0: The existing values in matrix C are preserved, and the result of the matrix multiplication is added to the existing values in C. It performs an accumulation operation where the new result is added to the existing values.
	
	//double complex (*prod)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*prod));

	double complex (*prod)[3] = calloc(3, sizeof(*prod));
	
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{ 
		if(ig_0 != mm)
		{
			for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //lower triangular matrix
			{
				extractSubmatrix( 
				(double complex*)G_sys_old, 3*tot_sub_vol, 3*tot_sub_vol, //[3*tot_sub_vol][3*tot_sub_vol]
        		G_old_33, 3, 3,
        		ig_0*3, jg_0*3  // Starting indices of the submatrix within the larger matrix
    			);
				
				extractSubmatrix(
				(double complex*)G_sys_new, 3*tot_sub_vol, 3*tot_sub_vol, //[3*tot_sub_vol][3*tot_sub_vol]
        		G_new_33, 3, 3,
        		ig_0*3, jg_0*3  // Starting indices of the submatrix within the larger matrix
    			);

				cblas_zgemm(
				CblasRowMajor,CblasNoTrans, CblasNoTrans, 
				3, 3, 3,  	// Dimensions of matrices
				&alpha_parameter, 					// Scaling factor for A*B
				G_old_33, 3, 					// Matrix A
				G_new_33, 3, 				// Matrix B
				&beta_parameter, 					// Scaling factor for C
				prod, 3								// Matrix C
				);

				insertSubmatrix(
        		(double complex*)G_sys_new, 3*tot_sub_vol, 3*tot_sub_vol,
        		(double complex*)prod, 3, 3,
        		ig_0*3, jg_0*3  // Starting indices in the larger matrix to insert the submatrix
    			);

				
				//for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
				//{
				//	int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				//	for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
				//	{
				//		int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
				//		G_sys_new[ig_0_2d][jg_0_2d] = G_sys_old[ig_0_2d][jg_0_2d]+prod[i_subG_0][j_subG_0]; 
				//		//printf("%e+i%e",creal(G_sys_new[ig_0_2d][jg_0_2d]),cimag(G_sys_new[ig_0_2d][jg_0_2d]));
				//	}
				//}
				
			}
		}
	}
	
	free(G_new_33);
	free(G_old_33);
	free(prod);
	*/

	/*
	cblas_zgemm(
							CblasRowMajor,CblasNoTrans, CblasNoTrans, 
							3, 3, 3,  	// Dimensions of matrices
							&alpha_parameter, 					// Scaling factor for A*B
							G_old_33,  3, 					// Matrix A
							G_new_33, 3, 				// Matrix B
							&beta_parameter, 					// Scaling factor for C
							prod, 3								// Matrix C
						);
	*/

	/*
	cblas_zgemm(
							CblasRowMajor,CblasNoTrans, CblasNoTrans, 
							3*tot_sub_vol, 3*tot_sub_vol, 3*tot_sub_vol,  	// Dimensions of matrices
							&alpha_parameter, 					// Scaling factor for A*B
							G_sys_old,  3*tot_sub_vol, 					// Matrix A
							G_sys_new, 3*tot_sub_vol, 				// Matrix B
							&beta_parameter, 					// Scaling factor for C
							prod, 3*tot_sub_vol								// Matrix C
						);
	*/
	
	
	// Next, solve all systems of equations for ii not equal to mm		
	for (int ig_0 = 0; ig_0 < tot_sub_vol; ig_0++) //tot_sub_vol
	{ 
		if(ig_0 != mm)
		{
			for(int i_subG_0 = 0; i_subG_0 < 3; i_subG_0++) // 3D coordinate positions
			{
				int ig_0_2d = (3*ig_0 + i_subG_0); // Set indices
				for (int jg_0 = 0; jg_0 < tot_sub_vol; jg_0++) //lower triangular matrix
				{
					for(int j_subG_0 = 0; j_subG_0 < 3; j_subG_0++) // 3D coordinate positions 
					{
						int jg_0_2d = (3*jg_0 + j_subG_0); // Set indices
						double complex G_sys_prod = 0.;
						for(int m_sub = 0;  m_sub < 3;  m_sub++)//loop for matricial multiplication
						{
							int mm_2d = (3*mm + m_sub);
							G_sys_prod += G_sys_old[ig_0_2d][mm_2d]*G_sys_new[mm_2d][jg_0_2d];
							//G_sys_prod += G_sys_old[ig_0_2d][mm_2d]*G_sys_new_transpose[jg_0_2d][mm_2d];	
						}
						G_sys_new[ig_0_2d][jg_0_2d] = G_sys_old[ig_0_2d][jg_0_2d] + pow(k,2)*alpha_0[mm]*G_sys_prod;
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

		// %%%%%%%%%%% G_new using Linear inversion using LAPACK %%%%%%%%%%%%%%%%%
		int info = LAPACKE_zgels(LAPACK_ROW_MAJOR,'N',3,3,3,A_iterative,3,b_iterative,3);
		G_sys_new_populator(tot_sub_vol, mm, jg_0, G_sys_new, b_iterative);

	} // end jg_0					

}


void core_solver(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[],double complex G_sys_new[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol]){


	double complex A_2d[3][3];
	for (int mm = 0; mm < tot_sub_vol; mm++) //tot_sub_vol
	{
		printf("%d - ",mm+1);

		double complex epsilon_s = (epsilon - epsilon_ref); // Scattering dielectric function

		A2d_solver(epsilon_s, mm, tot_sub_vol, delta_V_vector[mm], G_sys_old, A_2d, k);	
		remaining_pertubations(tot_sub_vol, mm, G_sys_old, G_sys_new, A_2d);
		offdiagonal_solver(tot_sub_vol, mm, k, alpha_0, G_sys_old, G_sys_new);

		memcpy(G_sys_old,G_sys_new,3*tot_sub_vol*3*tot_sub_vol*sizeof(double complex)); // Update G_old = G_new for next iteration.
	}//end mm loop 
}

/*
void iterative_solver(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[],double complex G_0[tot_sub_vol][tot_sub_vol][3][3], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol]){


	double complex (*G_sys_old)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys_old));
	double complex (*G_sys_new)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys_new));

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
	core_solver(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_sys_new, G_sys_old);

	memcpy(G_sys,G_sys_new,3*tot_sub_vol*3*tot_sub_vol*sizeof(double complex)); //Populate G_sys with G_new

	free(G_sys_new);
	free(G_sys_old);

}
*/
/*
void iterative_solver(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[],double complex G_sys_old[3*tot_sub_vol][3*tot_sub_vol], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol]){

	double complex (*G_sys_new)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys_new));
    
	core_solver(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_sys_new, G_sys_old);

	memcpy(G_sys,G_sys_new,3*tot_sub_vol*3*tot_sub_vol*sizeof(double complex)); //Populate G_sys with G_new

	free(G_sys_new);

}
*/

void iterative_solver(int tot_sub_vol, double complex epsilon, double complex epsilon_ref, double k, double delta_V_vector[], double complex alpha_0[], double complex G_sys[3*tot_sub_vol][3*tot_sub_vol],double k_0, double pi,double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3],char wave_type){

	double complex (*G_0)[tot_sub_vol][3][3] = calloc(tot_sub_vol, sizeof(*G_0)); 
	get_G0_matrix(tot_sub_vol, G_0, k_0, pi, epsilon_ref, modulo_r_i_j, r_i_j_outer_r_i_j, delta_V_vector, wave_type);
			
	double complex (*G_sys_old)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys_old));
	matrix_reshape(3, tot_sub_vol, G_sys_old, G_0); // this reshapes 2 4D matrices to 2 2D matrices, where G0 and eyeA are 4D and G_sys_old and eyeA_2d are the respective 2D matrices
	free(G_0);

	double complex (*G_sys_new)[3*tot_sub_vol] = calloc(3*tot_sub_vol, sizeof(*G_sys_new));
    
	core_solver(tot_sub_vol, epsilon, epsilon_ref, k, delta_V_vector, alpha_0, G_sys_new, G_sys_old);
	free(G_sys_old);
	memcpy(G_sys,G_sys_new,3*tot_sub_vol*3*tot_sub_vol*sizeof(double complex)); //Populate G_sys with G_new
	free(G_sys_new);

}