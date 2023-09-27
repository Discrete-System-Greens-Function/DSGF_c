#include "geometry/shared.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void set_delta_V_vector_T_vector(double T1, double T2, double delta_V_1, double delta_V_2, int tot_sub_vol, int const_N_subvolumes_per_object, double T_vector[], double delta_V_vector[]){


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
}

void set_delta_V_vector(double delta_V_1, double delta_V_2, int tot_sub_vol, int const_N_subvolumes_per_object, double delta_V_vector[]){

	for (int i_vec=0; i_vec<tot_sub_vol; i_vec++)
	{
		if (i_vec < const_N_subvolumes_per_object) // 2-body case
		{
			delta_V_vector[i_vec] = delta_V_1;
		}
		else
		{
			delta_V_vector[i_vec] = delta_V_2;
		}
	}  
}

void set_T_vector(double T1, double T2, int tot_sub_vol, int const_N_subvolumes_per_object, double T_vector[]){
	for (int i_vec=0; i_vec<tot_sub_vol; i_vec++)
	{
		if (i_vec < const_N_subvolumes_per_object) // 2-body case
		{
			T_vector[i_vec] = T1;
		}
		else
		{
			T_vector[i_vec] = T2;
		}
	}  
}




void write_delta_V_vector(double delta_V_1, double delta_V_2, int tot_sub_vol, int const_N_subvolumes_per_object){

	FILE* delta_V_vector_file = fopen("delta_V_vector.bin", "wb");
	// Write the entire array to the file
    //size_t numElements = tot_sub_vol;
	for (int i_vec=0; i_vec<tot_sub_vol; i_vec++)
	{
		if (i_vec < const_N_subvolumes_per_object) // 2-body case
		{
			//T_vector[i_vec] = T1;
			fwrite(&delta_V_1, sizeof(double), 1, delta_V_vector_file);
		}
		else
		{
			//T_vector[i_vec] = T2;
			fwrite(&delta_V_1, sizeof(double), 1, delta_V_vector_file);
		}
	} 	
	fclose(delta_V_vector_file); // Close the file
}
