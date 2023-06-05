#include "geometry/user_defined.h"
#include "file_utils.h"

void populate_R_user_defined(int tot_sub_vol, subvol shape_file[], double R[][3]){

	for (int i=0; i<tot_sub_vol;i++) //tot_sub_vol
	{
		R[i][0] = shape_file[i].x ;
		R[i][1] = shape_file[i].y ;
		R[i][2] = shape_file[i].z ;
	} 

}

void set_up_user_defined_geometry(int tot_sub_vol, int subvol_per_object, int N_bulk_object, double *d, double *T1, double *T2, double R[][3], double delta_V_vector[]){ //char *file_name_ud,
    
    char file_name_ud[256];
	read_geometry_user_defined(d,file_name_ud, T1, T2);
	printf("%s \n", file_name_ud);

	subvol shape_fileud[tot_sub_vol];
    char filename[256];
    sprintf(filename, "library/discretizations/user_defined/%s_discretization.txt", file_name_ud);	
	printf("%s \n",filename);
	printf("line 25 %s \n", file_name_ud);
	populate_subvol_struct(filename, tot_sub_vol, shape_fileud);
	printf("line 27 %s \n", file_name_ud);
	populate_R_user_defined(tot_sub_vol, shape_fileud, R);

	printf("line 28 %s \n", file_name_ud);
    //double delta_v_file_ud[tot_sub_vol];
    char filename_delta_v[256];
    sprintf(filename_delta_v, "library/discretizations/user_defined/%s_delta_V_vector.txt", file_name_ud);	
	

	printf("%s\n",filename_delta_v);
    populate_subvol_delta_v(filename_delta_v, tot_sub_vol, delta_V_vector);
	for (int i = 0; i < tot_sub_vol; i++)
	{
		printf("%le\n",delta_V_vector[i]);
	}
	
}