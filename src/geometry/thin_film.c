#include "geometry/thin_film.h"

void populate_R_thin_film(int tot_sub_vol, subvol shape_file[], double R[][3]){

	for (int i=0; i<tot_sub_vol;i++) //tot_sub_vol
	{
		R[i][0] = shape_file[i].x ;
		R[i][1] = shape_file[i].y ;
		R[i][2] = shape_file[i].z ;
	} 

}
