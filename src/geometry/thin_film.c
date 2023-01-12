#include "geometry/thin_film.h"
#include "file_utils.h"

void populate_R_thin_film(int tot_sub_vol, subvol shape_file[], double R[][3]){

	for (int i=0; i<tot_sub_vol;i++) //tot_sub_vol
	{
		R[i][0] = shape_file[i].x ;
		R[i][1] = shape_file[i].y ;
		R[i][2] = shape_file[i].z ;
	} 

}


void set_up_thin_film_geometry(int tot_sub_vol, int subvol_per_object, int N_bulk_objects, double *T1, double *T2, double *d, double *delta_V_1, double *delta_V_2, double R[][3]){

	double Lx, Ly, Lz;
	read_geometry_thin_films(d, &Lx, &Ly, &Lz, T1, T2);

	vol1 = Lx*Ly*Lz; // calculates the volume for membrane 1
	vol2 = vol1;     // defines the volume of membrane 2 = membrane 1
	*delta_V_1 = vol1/subvol_per_object; // defines the subvolumes' volume for membrane 1
	*delta_V_2 = vol2/subvol_per_object;  // defines the subvolumes' volume for membrane 1

	subvol shape_filetf[tot_sub_vol];

	int Lx_int, Ly_int, Lz_int, d_int;
	Lx_int = Lx*pow(10,9); 
	Ly_int = Ly*pow(10,9); 
	Lz_int = Lz*pow(10,9); 
	d_int = (*d)*pow(10,9);

	char filename[256];

	sprintf(filename, "discretizations/thin-film/%d_thin_films_Lx%dnm_Ly%dnm_Lz%dnm_d%dnm_N%d_discretization.txt", N_bulk_objects, Lx_int, Ly_int, Lz_int, d_int, tot_sub_vol);	
	populate_subvol_struct(filename, tot_sub_vol, shape_filetf);

	populate_R_thin_film(tot_sub_vol, shape_filetf, R);
}
