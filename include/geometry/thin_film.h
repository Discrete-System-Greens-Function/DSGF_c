#ifndef THIN_FILM_H
#define THIN_FILM_H

#include "functions_DSGF.h"
void set_up_thin_film_geometry(int tot_sub_vol, int subvol_per_object, int N_bulk_object, double *T1, double *T2, double *d, double *delta_V_1, double *delta_V_2, double R[][3]);


#endif /* THIN_FILM_H */
