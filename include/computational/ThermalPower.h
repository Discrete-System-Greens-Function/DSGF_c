#ifndef THERMALPOWER_H
#define THERMALPOWER_H
#include <complex.h>

void spectral_post_processing(int tot_sub_vol, int i_omega, int const_N_omega, double k_0, double h_bar, double k_b, double complex epsilon, double omega_value, double T_vector[],double delta_V_vector[], double const_N_subvolumes_per_object, double pi,double complex G_sys[3*tot_sub_vol][3*tot_sub_vol], double *sum_trans_coeff, double Q_subvol[tot_sub_vol][const_N_omega]);

#endif /* THERMALPOWER_H */
