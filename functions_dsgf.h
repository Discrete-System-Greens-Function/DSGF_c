// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Functions used in the Discrete System Green's Function 
// Developed by RETL Lab, Department of Mechanical Engineering, The University of Utah, USA.
// LAST UPDATE: MAY 31, 2022
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifndef __functions_dsgf_h //https://stackoverflow.com/questions/28361391/calling-define-from-another-file
#define __functions_dsgf_h

// ########## dsgf_functions.h ##########

// Constants: 
const double pi = 3.14159265359;          // pi number   
const double q = 1.602176634e-19;         // number of Joules per eV
const double h_bar = 1.054571817e-34;    // Planck's constant [J*s]
const double k_b = 1.38064852e-23;       // Boltzmann constant [J/K]
const double epsilon_0 = 8.8542e-12;     // Permittivity of free space [F/m]
const double c_0 = 299792458;            // Speed of light in vacuum [m/s]
double mu_0;                            // Permeability of free space [H/m]
double epsilon_ref = 1. ;             // dielectric function of the background reference medium

double complex epsilon;  // Definition of dielectric function

// Function for wave vector in free space calculation:
double k_0_function(double omega)         // function definition
{
    double k_0 = omega*sqrt(epsilon_0*mu_0);
    return k_0;                  // return statement
}

// Function for wave vector in background reference medium calculation:
double k_function(double omega)         // function definition
{
    double k = omega*sqrt(epsilon_ref*epsilon_0*mu_0);
    return k;                  // return statement
}

// Import user input as a const: https://stackoverflow.com/questions/8048203/is-it-possible-to-set-const-using-a-user-input
/*
// Instead of using a function, the .txt user_inputs.txt file is imported in the main code. 
// This way, the code will not need to be recompiled for every user modification.
// The function remains commented if we decide to modify it.
int get_user_inputs(char* material, char* geometry, double* d, double* radius, double* Lx, double* Ly, double* Lz, double* T1, double* T2, double* T_calc, char* solution, char* single_spectrum_analysis, char* save_A_matrix, char* save_G0_matrix, char* save_SGF_matrix, char* save_spectral_conductance, char* save_spectral_transmissivity)
{
    FILE *import_inputs; 
    char dirPathUserInputs[50];
    sprintf(dirPathUserInputs, "user_inputs.txt");
    import_inputs=fopen(dirPathUserInputs, "r"); 
    char buffer[50];
//    fscanf(import_inputs,"%s = %d\n",buffer, &N_bulk_spheres);
//    fscanf(import_inputs,"%s = %d\n",buffer, &N_sphere);
//    fscanf(import_inputs,"%s = %d\n",buffer, &N_omega); 
    fscanf(import_inputs,"%s = %s\n",buffer, material);
    fscanf(import_inputs,"%s = %s\n",buffer, geometry);
    fscanf(import_inputs,"%s = %lf\n",buffer, &d);
    fscanf(import_inputs,"%s = %lf\n",buffer, &radius);
    fscanf(import_inputs,"%s = %lf\n",buffer, &Lx);
    fscanf(import_inputs,"%s = %lf\n",buffer, &Ly);
    fscanf(import_inputs,"%s = %lf\n",buffer, &Lz);
    fscanf(import_inputs,"%s = %lf\n",buffer, &T1);
    fscanf(import_inputs,"%s = %lf\n",buffer, &T2);
    fscanf(import_inputs,"%s = %lf\n",buffer, &T_calc);
    fscanf(import_inputs,"%s = %s\n",buffer, solution);
    fscanf(import_inputs,"%s = %c\n",buffer, &single_spectrum_analysis);
    fscanf(import_inputs,"%s = %c\n",buffer, &save_A_matrix);
    fscanf(import_inputs,"%s = %c\n",buffer, &save_G0_matrix);
    fscanf(import_inputs,"%s = %c\n",buffer, &save_SGF_matrix);
    fscanf(import_inputs,"%s = %c\n",buffer, &save_spectral_conductance);
    fscanf(import_inputs,"%s = %c",buffer, &save_spectral_transmissivity);
    fclose(import_inputs);
//    const int const_user_input_N_sphere = N_sphere;
    //return const_user_input_N_sphere;
}
*/

// The parameters N_subvolumes_per_object, N_bulk_objects, and N_omega need to treated differently in the code.
// The reason is because they are used to structure the solution matrices and are required to be non-variable parameters.  
int get_N_subvolumes_per_object()
{
    FILE *import_N_subvolumes_per_object; 
    char dirPathN_subvolumes_per_object[50];
    sprintf(dirPathN_subvolumes_per_object, "N_subvolumes_per_object.txt");
    import_N_subvolumes_per_object=fopen(dirPathN_subvolumes_per_object, "r"); 
    fscanf(import_N_subvolumes_per_object,"%d", &N_subvolumes_per_object);
    fclose(import_N_subvolumes_per_object);
    const int const_user_input_N_subvolumes_per_object = N_subvolumes_per_object;
    return const_user_input_N_subvolumes_per_object;
}

int get_N_bulk_objects()
{
    FILE *import_N_bulk_objects; 
    char dirPathN_bulk_objects[50];
    sprintf(dirPathN_bulk_objects, "N_bulk_objects.txt");
    import_N_bulk_objects=fopen(dirPathN_bulk_objects, "r"); 
    fscanf(import_N_bulk_objects,"%d", &N_bulk_objects);
    fclose(import_N_bulk_objects);
    const int const_user_input_N_bulk_objects = N_bulk_objects;
    return const_user_input_N_bulk_objects;
}

int get_N_omega()
{
    FILE *import_N_omega; 
    char dirPathN_omega[50];
    sprintf(dirPathN_omega, "N_omega.txt");
    import_N_omega=fopen(dirPathN_omega, "r"); 
    fscanf(import_N_omega,"%d", &N_omega);
    fclose(import_N_omega);
    const int const_user_input_N_omega = N_omega;
    return const_user_input_N_omega;
}


// Function for the thermal volumes calculation:
double vol_sphere(double radius)               // function definition, radius is a variable inside the function   
{
    //float diameter= 2*radius;
    double volume = (4./3)*pi*pow(radius,3);  // volume for sphere 
    return volume;                            // return statement
}

// Definition used to import the discretization positions
typedef struct node {
    float x;
    float y;
    float z;
} subvol;

// Function a for free space green's function calculation:
double a_j_function(double delta_v_vector)         
{
    double a = pow( (3.* delta_v_vector)/(4*pi), 1./3);  // volume for sphere 
    return a;                  // return statement
}

// Function for the mean energy of an electromagnetic state calculation:
double theta_function(double omega, double T)    
{
    double theta = (h_bar*omega)/(exp(h_bar*omega/(k_b*T))-1) ;   
    return theta;                  // return statement
}

// Function for the calculation of the derivative of the mean energy of an electromagnetic state with respect to the temperature
double  dtheta_dT_function(double omega,double T_calc)
{
    double psi = h_bar*(omega)/(k_b*T_calc);     // Input for exponential dtheta_dT function[dimensionless]
    double dtheta_dT = (k_b*pow(psi,2)*exp(psi))/(pow(exp(psi)-1.,2)); // [J/K]
    return dtheta_dT;           // return statement
}

/*          
double Q_omega_subvol_function(double theta_1,double theta_2, double trans)    
{
    double Q_omega_subvol = (1./(2.*pi))*(fabs(theta_2-theta_1)*trans) ;   
    return Q_omega_subvol;                  // return statement
}
*/

// How to measure memory usage inside my program? (getrusage): https://youtu.be/Os5cK0H8EOA
long get_mem_usage() 
{
    struct rusage myusage;
    getrusage(RUSAGE_SELF,&myusage); // getrusage() https://man7.org/linux/man-pages/man2/getrusage.2.html
    return myusage.ru_maxrss;           // return statement
   
}

// Geometry parameters' definitions
double radius1; 
double radius2; 
double vol1;
double vol2;
double delta_V_1; 
double delta_V_2;

// Frequency parameters' definitions
int N_range;
int omega_range;
double omega_value; 

int bulk; //Indice used to thermal power dissipated 

double denom1, denom2 ; // used in G^0_ij function

int ipack, gpack; //lapacke counters
int ig_0_2d,jg_0_2d,mm_2d,mm_2d_n, mm_sub, mm_sub_n; // Set indices used in "iterative" solution

//double theta_1, theta_2; // Definition for the mean energy of an electromagnetic state of the thermal objects
//double dtheta_dT; // function used to calculate conductance

double trapz; // Definition for trapezoidal integration. Used in total conductance
//double trapz_Q; // Definition for trapezoidal integration. Used in total power dissipated


// Code's output files definitions
char matrices_folder[100];
char frequency_folder[100];
char sep_distance_folder[100];
char spectral_transmissivity_folder[100];
char spectral_conductance_folder[100];
char dirPathpos_processing_summary_FileName[260];

struct stat st = {0}; // used in the condition to create new directories

// ########## end dsgf_functions.h ##########

#endif // __functions_dsgf_h

