#include "geometry/sample.h"
#include "file_utils.h"
#include <string.h> // library used to concatenate 2 strings https://stackoverflow.com/questions/46612504/creating-directories-and-files-using-a-loop-in-c


void populate_R_sample(int tot_sub_vol, int subvol_per_object, double origin1[], double origin2[], double delta_V_1, double delta_V_2, subvol shape_file1[], subvol shape_file2[], double R[][3]){

	for (int i=0; i<tot_sub_vol;i++) //tot_sub_vol
	{
		if(i<subvol_per_object)
		{
			R[i][0] = shape_file1[i].x *pow(delta_V_1,1./3) + origin1[0];
			R[i][1] = shape_file1[i].y *pow(delta_V_1,1./3) + origin1[1];
			R[i][2] = shape_file1[i].z *pow(delta_V_1,1./3) + origin1[2];
		}
		else
		{
			R[i][0] = shape_file2[i-subvol_per_object].x* pow(delta_V_2,1./3) + origin2[0]; 
			R[i][1] = shape_file2[i-subvol_per_object].y*pow(delta_V_2,1./3) + origin2[1]; 
			R[i][2] = shape_file2[i-subvol_per_object].z*pow(delta_V_2,1./3) + origin2[2]; 
		}

	}   

}

void set_up_sample_geometry(double pi, int tot_sub_vol, int subvol_per_object, int N_subvolumes_per_object_2, double *T1, double *T2, double *d, double *delta_V_1, double *delta_V_2, double R[][3], char *geometry_1, char *geometry_2){

	double radius1, radius2;
	
	//read_geometry_sphere(d, &radius, T1, T2);
	read_geometry_sample(geometry_1, &radius1, geometry_2, &radius2, d, T1, T2);
	//printf("sample simulation");
	double vol1; 
	if (strcmp(geometry_1, "sphere") == 0){ vol1 = vol_sphere(radius1, pi); } // calls function that calculates the volume for the sphere 1
	if (strcmp(geometry_1, "cube") == 0){ vol1 = pow(radius1, 3); } // calls function that calculates the volume for the sphere 1
	double vol2;
	if (strcmp(geometry_2, "sphere") == 0){ vol2 = vol_sphere(radius2, pi); } // calls function that calculates the volume for the sphere 2
	if (strcmp(geometry_2, "cube") == 0){ vol2 = pow(radius2, 3); } // calls function that calculates the volume for the sphere 1
	*delta_V_1 = vol1/subvol_per_object; // defines the subvolumes' volume for sphere 1
	*delta_V_2 = vol2/N_subvolumes_per_object_2; // defines the subvolumes' volume for sphere 2

	subvol shape_file1[subvol_per_object];
	subvol shape_file2[N_subvolumes_per_object_2];
	
	char file_name1[256];
	sprintf(file_name1, "library/discretizations/%s/%s_%d.txt",geometry_1,geometry_1,subvol_per_object); // path where the file is stored
	populate_subvol_struct(file_name1, subvol_per_object, shape_file1);	

	char file_name2[256];
	sprintf(file_name2, "library/discretizations/%s/%s_%d.txt",geometry_2,geometry_2,N_subvolumes_per_object_2); // path where the file is stored
	populate_subvol_struct(file_name2, N_subvolumes_per_object_2, shape_file2);

	double origin1[3] =  {0,0,0};
	double origin_x; 
	if (strcmp(geometry_2, "sphere") == 0){ 
		origin_x = radius1+*d + radius2; 
	} // calls function that calculates the volume for the sphere 2
	if (strcmp(geometry_2, "cube") == 0){ origin_x = radius1/2+*d+radius2/2; } // calls function that calculates the volume for the sphere 1
	double origin2[3]= {origin_x,0,0};

	populate_R_sample(tot_sub_vol, subvol_per_object, origin1, origin2, *delta_V_1, *delta_V_2, shape_file1, shape_file2 , R);

}
