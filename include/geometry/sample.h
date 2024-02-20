#ifndef SAMPLE_H
#define SAMPLE_H

#include <functions_DSGF.h>

void set_up_sample_geometry(double pi, int tot_sub_vol, int subvol_per_object, int N_subvolumes_per_object_2, double *T1, double *T2, double *d, double *delta_V_1, double *delta_V_2, double R[][3], char *geometry_1, char *geometry_2);

#endif /* SAMPLE_H */
