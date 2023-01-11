#include "geometry/sphere.h"

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
