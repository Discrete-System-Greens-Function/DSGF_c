// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// User inputs definitions in DSGF framework
// Developed by RETL Lab, Department of Mechanical Engineering, The University of Utah, USA.
// LAST UPDATE: JUNE 1, 2022  
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifndef __user_inputs_h //https://stackoverflow.com/questions/28361391/calling-define-from-another-file
#define __user_inputs_h

// %%%%%%%%%%%% Simulations parameters %%%%%%%%%%%%%%%%%%

int N_subvolumes_per_object; // Number of subvolumes per thermal object. Spheres: 1, 8, 32, 136, 280, 552, 912, 1472. Flat surfaces: 320, 640, 1280
int N_bulk_objects; //Number of bulk objects: 2 
int N_omega; //  Number of frequencies to evaluate
char uniform_spectra;

// %%%%%%  Material, geometry and temperature %%%%%%%%%%%

double epsilon_ref; // dielectric function of the background reference medium
char material[]; // Material options: SiC or SiO2
char geometry[]; //Geometry options: sphere, thin-films or user-defined
char discretization_thin_film[260]; //discretization_thin_film = 2_thin_films_Lx200nm_Ly1um_Lz200nm_d150nm_N640_discretization
//char file_name_ud[]; // user_defined

double d; // Separation distance

double radius; // Perfect same-sized spheres

// Dimensions for membranes/thin films 
double Lx; // Depth: 200.e-9
double Ly; // Length: 1.e-6
double Lz; // Thickness: 200.e-9

double T1; // Temperature in K modified to calculate thermal power dissipated
double T2; // Temperature in K
double T_calc; // Temperature at which conductance is calculated [K]
int const N_Tcalc = 5;
double Tcalc_vector[N_Tcalc]; // Multiple temperatures at which conductance is calculated [K]
// %%%%%%%%%%%%%%%%%%% Solution %%%%%%%%%%%%%%%%%%%%%

char solution; // Solution method: D for "direct" or I for "iterative"

char single_analysis ='y'; // Initial input for spectrum range analysis, used to compute memory_usage.
char single_spectrum_analysis; // Analysis of one or a range of frequencies

double G_12_total_SGF_from_omega; // total conductance
//double G_12_total_SGF_from_lambda;


// %%%%%%%%%%%%%%%%%% Save data %%%%%%%%%%%%%%%%%%%%

char save_A_matrix;
char save_G0_matrix;
char save_SGF_matrix;
char save_spectral_conductance;
char save_total_conductance;
char save_power_dissipated_spectral_subvolumes;
char save_power_dissipated_total_subvolumes;
char save_power_dissipated_spectral_bulk;
char save_power_dissipated_total_bulk;
char save_power_density_total_subvolumes;
char save_spectral_transmissivity;
char save_power_dissipated;
char multiple_conductance_temperatures;

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#endif // __user_inputs_h
