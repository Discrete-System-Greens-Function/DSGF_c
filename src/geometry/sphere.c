#include "geometry/sphere.h"
#include "file_utils.h"

void populate_R_sphere(int tot_sub_vol, int subvol_per_object, double origin1[], double origin2[], double delta_V_1, double delta_V_2, subvol shape_file[], double R[][3]){

	for (int i=0; i<tot_sub_vol;i++) //tot_sub_vol
	{
		if(i<subvol_per_object)
		{
			R[i][0] = shape_file[i].x *pow(delta_V_1,1./3) + origin1[0];
			R[i][1] = shape_file[i].y *pow(delta_V_1,1./3) + origin1[1];
			R[i][2] = shape_file[i].z *pow(delta_V_1,1./3) + origin1[2];
		}
		else
		{
			R[i][0] = shape_file[i-subvol_per_object].x* pow(delta_V_2,1./3) + origin2[0]; 
			R[i][1] = shape_file[i-subvol_per_object].y*pow(delta_V_2,1./3) + origin2[1]; 
			R[i][2] = shape_file[i-subvol_per_object].z*pow(delta_V_2,1./3) + origin2[2]; 
		}

	}   

}

void set_up_sphere_geometry(double pi, int tot_sub_vol, int subvol_per_object, double *T1, double *T2, double *d, double *delta_V_1, double *delta_V_2, double R[][3]){

	double radius;

	read_geometry_sphere(d, &radius, T1, T2);
	double radius1 = radius; // perfect same-sized spheres
	double radius2 = radius; // perfect same-sized spheres
	double vol1 = vol_sphere(radius1, pi); // calls function that calculates the volume for the sphere 1
	double vol2 = vol_sphere(radius2, pi); // calls function that calculates the volume for the sphere 2
	*delta_V_1 = vol1/subvol_per_object; // defines the subvolumes' volume for sphere 1
	*delta_V_2 = vol2/subvol_per_object; // defines the subvolumes' volume for sphere 2

	subvol shape_file[subvol_per_object];
	
	char file_name[256];
	sprintf(file_name, "library/discretizations/sphere/sphere_subvol_%d.txt",subvol_per_object); // path where the file is stored
	populate_subvol_struct(file_name, subvol_per_object, shape_file);	

	double origin1[3] = {radius1,radius1,radius1};
	double origin2[3]= {origin1[0]+radius1+(*d)+radius2,origin1[1]+radius2-radius1,origin1[2]+radius2-radius1};

	populate_R_sphere(tot_sub_vol, subvol_per_object, origin1, origin2, *delta_V_1, *delta_V_2, shape_file, R);

}
