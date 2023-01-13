#ifndef GREENSFUNCTION_H
#define GREENSFUNCTION_H

#include <complex.h>

void setup_G_0_matrices(int tot_sub_vol, double modulo_r_i_j[tot_sub_vol][tot_sub_vol], double complex r_i_j_outer_r_i_j[tot_sub_vol][tot_sub_vol][3][3], double R[][3]);

#endif /* GREENSFUNCTION_H */
