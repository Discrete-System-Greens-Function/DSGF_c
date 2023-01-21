#ifndef DIRECT_SOLVER_H
#define DIRECT_SOLVER_H

#include <complex.h>

void A_b_direct_populator(int tot_sub_vol, double complex A[tot_sub_vol][tot_sub_vol][3][3], double complex G_0[tot_sub_vol][tot_sub_vol][3][3], double complex A_direct[3*3*tot_sub_vol*tot_sub_vol], double complex b_direct[3*3*tot_sub_vol*tot_sub_vol]);


#endif /* DIRECT_SOLVER_H */
