%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Non-uniform spectrum generation for the c version of DSGF solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%
% Constants %
%%%%%%%%%%%%%

q = 1.602176634e-19;        % Number of joules per eV [J/eV]
h_bar = 1.054571817e-34;    % Planck's constant [J*s]
k_b = 1.38064852e-23;       % Boltzmann constant [J/K]
epsilon_0 = 8.8542e-12;     % Permittivity of free space [F/m]
mu_0 = (4*pi)*(10^-7);      % Permeability of free space [H/m]
c_0 = 299792458;            % Speed of light in vacuum [m/s]

%%%%%%%%%%%%%%%%%%%%%%%
% Select the material %
%%%%%%%%%%%%%%%%%%%%%%%
%{
material = 'SiO2'; % in terms of wavelength
initial_lambda = 5.e-6; 
final_lambda = 25.e-6;
Nref = 100; 
lambda = linspace(start_lambda, end_lambda, Nref); % [m]
c_0 = 299792458;            % Speed of light in vacuum [m/s]
omega = (2*pi*c_0./lambda); % [rad/s]
omega = sort(omega, 'ascend');
N_omega = length(omega);
lambda_i = initial_lambda*10^(6);
lambda_f = final_lambda*10^(6);
spectra = [ material '_' num2str(N_omega) '_uniform_' num2str(lambda_i) '_' num2str(lambda_f) '_um.csv']; %if local
%}
%{
material = 'SiO2'; % in terms of radial frequency
initial = 7.5e13; 
final = 3.5e14;
Nref = 100;
omega = linspace(initial, final, Nref);
N_omega = length(omega);
wi = initial*10^(-12);
wf = final*10^(-12);
spectra = [ material '_' num2str(N_omega) '_uniform_' num2str(wi) '_' num2str(wf) '_Trad_s.csv']; %if local
%}

material = 'SiC'; % in terms of radial frequency
initial = 1.4e14; 
final = 1.9e14;
Nref = 100;
omega = linspace(initial, final, Nref);
N_omega = length(omega);
wi = initial*10^(-12);
wf = final*10^(-12);
spectra = [ material '_' num2str(N_omega) '_uniform_' num2str(wi) '_' num2str(wf) '_Trad_s.csv']; %if local

%{
material = 'Si3N4'; % in terms of radial frequency
initial = 2.e13; 
final = 3.e14;
Nref = 100;
omega = linspace(initial, final, Nref);
N_omega = length(omega);
wi = initial*10^(-12);
wf = final*10^(-12);
spectra = [ material '_' num2str(N_omega) '_uniform_' num2str(wi) '_' num2str(wf) '_Trad_s.csv']; %if local
%}

%lambda = 2*pi*c_0./omega;
 
% Save spectra vector to a .csv file
back = cd;
%saveDir = fullfile(back,'../../library/Non_uniform_spectra');
%spectra = [saveDir '/' material '_non_uniform_spectra_' num2str(N_omega)  '.csv']; %save in DSGF library folder
writematrix(omega', spectra);

%{
% Split spectra into multiple simulations. 
% Designed to 100 frequecies. Modify the code according to the range to be split. 
spectra_1 = [ material '_' num2str(N_omega/5) '_uniform_split_1.csv']; %if local
writematrix(omega(1:20)', spectra_1);
spectra_2 = [ material '_' num2str(N_omega/5) '_uniform_split_2.csv']; %if local
writematrix(omega(21:40)', spectra_2);
spectra_3 = [ material '_' num2str(N_omega/5) '_uniform_split_3.csv']; %if local
writematrix(omega(41:60)', spectra_3);
spectra_4 = [ material '_' num2str(N_omega/5) '_uniform_split_4.csv']; %if local
writematrix(omega(61:80)', spectra_4);
spectra_5 = [ material '_' num2str(N_omega/5) '_uniform_split_5.csv']; %if local
writematrix(omega(81:100)', spectra_5);
%}

if strcmp('SiO2',material)
    
    % Dielectric function approximated with a Lorentz model
    epsilon_inf = 2.03843;                          % High-frequency limit of permittivity
    h_bar_omega_0 = [0.05624, 0.09952, 0.13355];    % [eV]
    omega_0 = (h_bar_omega_0*q)/(h_bar);            % [rad/s]
    %lambda_0 = [22.04432, 12.45818, 9.28364]*1e-6;  % [m]
    S = [0.93752, 0.05050, 0.60642];                % [-]
    Gamma = [0.09906, 0.05511, 0.05246];            % [-]

    summation = zeros(length(omega),1); % Preallocate vector
    for k = 1:length(S)
        num = S(k);
        term_1 = (omega'./omega_0(k)).^2;
        term_2 = 1i*Gamma(k).*(omega'./omega_0(k));
        denom = 1 - term_1 - term_2;
        summation = summation + (num./denom);
    end
    epsilon = epsilon_inf + summation;
    
elseif strcmp('SiC',material)
    
    % Dielectric function approximated with a Lorentz model
    omega_TO = 1.494e14;              % Transverse optical phonon frequency [rad/s]
    omega_LO = 1.825e14;              % Longitudinal optical phonon frequency [rad/s]
    Gamma_D = 8.966e11;               % Damping constant [rad/s]
    epsilon_inf = 6.7;                % High-frequency limit of permittivity
    epsilon = epsilon_inf*( ((omega.^2) - (omega_LO^2) + 1i*Gamma_D*omega)./ ((omega.^2) - (omega_TO^2) + 1i*Gamma_D*omega) );
    
    
elseif strcmp('SiC-poly',material)
    
    % Dielectric function approximated with a Lorentz model
    omega_TO = 1.486e14; %1.494e14;              % Transverse optical phonon frequency [rad/s]
    omega_LO = 1.801e14; %1.825e14;              % Longitudinal optical phonon frequency [rad/s]
    Gamma_D = 3.767e12; %8.966e11;               % Damping constant [rad/s]
    epsilon_inf = 8; %6.7;                % High-frequency limit of permittivity
    epsilon = epsilon_inf*( ((omega.^2) - (omega_LO^2) + 1i*Gamma_D*omega)./ ((omega.^2) - (omega_TO^2) + 1i*Gamma_D*omega) );
        
elseif strcmp('Si3N4',material)
   
    M = 5; % Number of oscillations
    epsilon_real = [7.582, 6.754, 6.601, 5.430, 4.601, 4.562];
    epsilon_imag = [0, 0.3759, 0.0041, 0.1179, 0.2073, 0.0124];
    epsilon = epsilon_real + 1i.*epsilon_imag;
    epsilon_inf = epsilon(M+1);
    omega_T = (2*pi).*(1e12).*[13.913, 15.053, 24.521, 26.440, 31.724];
    Gamma = (2*pi).*(1e12).*[5.810, 6.436, 2.751, 3.482, 5.948];
    alpha = [0.0001, 0.3427, 0.0006, 0.0002, 0.0080];

    summation = 0;

    for jj = 1:M
        
        num_1 = ((omega_T(jj))^2) - omega'.^2;
        denom_1 = omega'.*Gamma(jj);
        Gamma_prime = Gamma(jj).*exp(-alpha(jj).*((num_1./denom_1).^2));
    
        delta_epsilon = epsilon(jj) - epsilon(jj+1);
    
        num_2 = delta_epsilon*((omega_T(jj))^2);
        denom_2 = ((omega_T(jj))^2) - omega'.^2 - 1i.*omega'.*Gamma_prime;
    
        summation = summation + num_2./denom_2;
    
    end

    epsilon = epsilon_inf + summation;
    
end

epsilon_fig = figure(1);
semilogx(omega, real(epsilon), '-o', 'linewidth', 2) % '--'
hold on
semilogx(omega, imag(epsilon), '-x', 'linewidth', 2)
%plot(lambda*1e6, real(epsilon), lambda*1e6, imag(epsilon), '--', 'linewidth', 2)
xlabel('Frequency [rad/s]', 'fontsize', 12)
%xlabel('Wavelength, \lambda [\mum]', 'fontsize', 12)
ylabel('Dielectric function, \epsilon', 'fontsize', 12)
title(['Dielectric function of ' material ', N_o_m_e_g_a = ' num2str(N_omega)], 'fontsize', 16)

%title(['Dielectric function of Si_3N_4, N_o_m_e_g_a = ' num2str(N_omega)], 'fontsize', 16)
legend('Real part', 'Imaginary part', 'location', 'best')
set(gca, 'fontsize', 16)
axis tight
grid on
hold off        
%saveas(epsilon_fig, 'fig_SiC_dieletric_function.fig');        