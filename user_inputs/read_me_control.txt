The user can modify inputs and select outputs in the control.txt file.

geometry = user_defined             Choose between sphere or user_defined
number_subvolumes_per_object = 36   Number of subvolumes per thermal object, restricted by the .txt discretization files
number_bulk_objects = 2             Number of thermal objects
material = SiO2                     Choose between SiO2, SiC, SiN, or u-SiC
epsilon_ref = 1                     Background reference medium
number_omega = 20                   Define a number of frequencies
uniform_spectra = S                 Choose between Y for uniform, N for non-uniform or S for split
solution = D                        Choose between D for direct or I for iterative
multithread = N                     Choose between N for run in serial or Y for run in parallel 
single_spectrum_analysis = N        Choose Y if a single frquency will be analyzed. Otherwise, keep as N
wave_type = T                       Choose T for total, P for propagating, and E for evanescent

Outputs are selected by Y or N to indicate yes or no.

save_spectral_conductance = Y
save_total_conductance = Y
save_power_dissipated_spectral_subvolumes = N
save_power_dissipated_total_subvolumes = N
save_power_dissipated_spectral_bulk = N
save_power_dissipated_total_bulk = N
save_power_density_total_subvolumes = N
save_spectral_transmissivity = N

For the geometry file,
    - inputs in user_inputs/Geometry/sphere.txt for spheres	
		d = 150.e-9			--- distance between the bulk objects
		radius = 50.e-9		--- radius of the bulk objects
		T1 = 300			--- Temperature for bulk object 1
		T2 = 400			--- Temperature for bulk object 2
    - inputs in user_inputs/Geometry/user_defined.txt 
        d = 500.e-9         --- distance between the bulk objects
        file_name = 2_thin_films_Lx500nm_Ly500nm_Lz500nm_d500nm_N72     --- file name
        T1 = 400            --- Temperature for bulk object 1
        T2 = 300            --- Temperature for bulk object 2