#include "debugging_utils.h"

void print_matrix_complex(int row, int col, double complex array[row][col]){
	
	for (int i = 0; i < row; i++){
		for (int j = 0; j < col; j++){
			printf("%f + i%f | ", creal(array[i][j]), cimag(array[i][j]));
		}
		printf("\n");
		printf("\n");
	}

}

void print_matrix(int row, int col, double array[row][col]){

	for (int i = 0; i < row; i++){
		for (int j = 0; j < col; j++){
			printf("%f | ", creal(array[i][j]));
		}
		printf("\n");
		printf("\n");
	}
}

void print_matrix_range_complex(int row, int col, int start_row, int end_row, int start_col, int end_col, double complex array[row][col]){
	
	for (int i = start_row; i < end_row; i++){
		for (int j = start_col; j < end_col; j++){
			printf("%f + i%f | ", creal(array[i][j]), cimag(array[i][j]));
		}
		printf("\n");
		printf("\n");
	}

}

void print_matrix_range(int row, int col, int start_row, int end_row, int start_col, int end_col, double array[row][col]){
		
	for (int i = start_row; i < end_row; i++){
		for (int j = start_col; j < end_col; j++){
			printf("%f | ", creal(array[i][j]));
		}
		printf("\n");
		printf("\n");
	}

}


void print_4d_internal_complex(int row, int col, int minor_row, int minor_col, int row_index, int col_index, double complex array[row][col][minor_row][minor_col]){
	print_matrix_complex(minor_row, minor_col, array[row_index][col_index]);
}

void p_int(int row_index, int col_index, double complex array[16][16][3][3]){
	print_4d_internal_complex(16,16,3,3,row_index, col_index, array);
}
