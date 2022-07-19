#ifndef __debugging_utils_h__
#define __debugging_utils_h__

#include <complex.h>
#include <stdio.h>

void print_matrix_complex(int row, int col, double complex array[row][col]);

void print_matrix(int row, int col, double array[row][col]);

void print_matrix_range_complex(int row, int col, int start_row, int end_row, int start_col, int end_col, double complex array[row][col]);

void print_matrix_range(int row, int col, int start_row, int end_row, int start_col, int end_col, double array[row][col]);

void print_4d_internal_complex(int row, int col, int minor_row, int minor_col, int row_index, int col_index, double complex array[row][col][minor_row][minor_col]);

void p_int(int row_index, int col_index, double complex array[16][16][3][3]);

#endif
