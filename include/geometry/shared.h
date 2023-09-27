#ifndef SHARED_H
#define SHARED_H

void set_delta_V_vector_T_vector(double T1, double T2, double delta_V_1, double delta_V_2, int tot_sub_vol, int const_N_subvolumes_per_object, double T_vector[], double delta_V_vector[]);

void set_delta_V_vector(double delta_V_1, double delta_V_2, int tot_sub_vol, int const_N_subvolumes_per_object, double delta_V_vector[]);

void set_T_vector(double T1, double T2, int tot_sub_vol, int const_N_subvolumes_per_object, double T_vector[]);

void write_T_vector(double T1, double T2, int tot_sub_vol, int const_N_subvolumes_per_object);

#endif /* SHARED_H */
