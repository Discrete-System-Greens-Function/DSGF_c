
%%%%%%%%%%%%%
% Constants %
%%%%%%%%%%%%%

q = 1.602176634e-19;        % Number of joules per eV [J/eV]
h_bar = 1.054571817e-34;    % Planck's constant [J*s]
k_b = 1.38064852e-23;       % Boltzmann constant [J/K]
epsilon_0 = 8.8542e-12;     % Permittivity of free space [F/m]
mu_0 = (4*pi)*(10^-7);      % Permeability of free space [H/m]
c_0 = 299792458;            % Speed of light in vacuum [m/s]

material = 'SiC';

%Non-uniform spectrum for SiC
%New spectra
omega_i = 1.e14;%7.5e13;
N_ref_1 = 3;
omega_ref_1 = 1.2e14;
N_ref_2 = 8; %35;
omega_ref_2 = 1.3e14 ;
N_ref_3 = 18; %35;
omega_ref_3 = 1.5e14 ;
N_ref_4 = 28; %35;
omega_ref_4 = 1.6e14 ;
N_ref_5 = 25; %40;
omega_ref_5 = 1.7e14;%1.82e14 ;
N_ref_6 = 18;
omega_ref_6 = 1.8e14 ;%2.3e14 ;
N_ref_7 = 6;
omega_f = 2.e14 ;%3.85e14;
 
%{
% old spectra
omega_i = 7.5e13;
N_ref_1 = 5;
omega_ref_1 = 1.3e14;
N_ref_2 = 18; %35;
omega_ref_2 = 1.47e14 ;
N_ref_3 = 6; %35;
omega_ref_3 = 1.51e14 ;
N_ref_4 = 10; %35;
omega_ref_4 = 1.64e14 ;
N_ref_5 = 28; %40;
omega_ref_5 = 1.82e14 ;
N_ref_6 = 4;
omega_ref_6 = 2.3e14 ;
N_ref_7 = 3;
omega_f = 3.85e14;
%}

N_ref_total = N_ref_1+N_ref_2+N_ref_3+N_ref_4+N_ref_5+N_ref_6+N_ref_7;
omega_range_ref_1 = linspace(omega_i, omega_ref_1, N_ref_1);
omega_range_ref_2 = linspace(omega_ref_1, omega_ref_2, N_ref_2);
omega_range_ref_3 = linspace(omega_ref_2, omega_ref_3, N_ref_3);
omega_range_ref_4 = linspace(omega_ref_3, omega_ref_4, N_ref_4);
omega_range_ref_5 = linspace(omega_ref_4, omega_ref_5, N_ref_5);
omega_range_ref_6 = linspace(omega_ref_5, omega_ref_6, N_ref_6);
omega_range_ref_7 = linspace(omega_ref_6, omega_f, N_ref_7);
%size(omega_range_1);

omega = [omega_range_ref_1,omega_range_ref_2,omega_range_ref_3,omega_range_ref_4,omega_range_ref_5,omega_range_ref_6,omega_range_ref_7];
omega = unique (omega);
N_omega = length(omega);
  
%N_omega = N_ref_total-6;
%lambda = 2*pi*c_0./omega;
 
% Save spectra vector to a .csv file
back = cd;
saveDir = fullfile(back,'../../library/Non_uniform_spectra');
spectra = [saveDir '/' material '_non_uniform_spectra_' num2str(N_omega)  '.csv'];
writematrix(omega', spectra);
        
        