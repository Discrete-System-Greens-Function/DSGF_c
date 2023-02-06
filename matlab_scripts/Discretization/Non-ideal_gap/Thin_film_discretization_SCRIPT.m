%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lindsay Walter
% u0928979
% Thin film discretization
% Updated 11/5/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all
format long
profile on
disp(['Running MATLAB script ' mfilename])

% This code creates discretizations for a system of two thin films
% separated by a gap.  Results are output in either .txt format or in an
% Excel file based on user input.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******************************USER INPUTS********************************%

%%%%%%%%%%%%%%%%%%%%%%
% Set results export %
%%%%%%%%%%%%%%%%%%%%%%

% Directory where results files will be saved (string)
back = cd;
saveDir = fullfile(back, '../../../library/discretizations/thin-film/');

% Save discretization in an Excel file? (0 = "no", 1 = "yes")
save_Excel = 0;

% Save discretization in an .txt file? (0 = "no", 1 = "yes")
save_txt = 1;

% Save figures? (0 = "no", 1 = "yes")
save_fig = 1;


%%%%%%%%%%%%%%%%%%
% Figure options %
%%%%%%%%%%%%%%%%%%

% Show figure axes? (0 = "no", 1 = "yes")
show_axes = 1;

% Show figure titles? (0 = "no", 1 = "yes")
show_titles = 1;


%%%%%%%%%%%%%%%%%%%%
% Dimension inputs %
%%%%%%%%%%%%%%%%%%%%

d = 100e-9;  % Distance between the two films [m]

% Thin film #1
Lx_1 = 1000e-9;         % Length of film #1   -- Previously tested: 200e-9, 400e-9, 800e-9
Ly_1 = 0.5e-6;           % Width of film #1
Lz_1 = 20e-9;         % Thickness of film #1 -- Previously tested: 200e-9
origin_1 = [0,0,0];    % Point of back, left, bottom of cube [x,y,z]
cord_1 = 'min';          % String of coordinate along which direction mesh is to be refined 
mesh_1 = 1;            % Number of subvolumes across the dimension given in the input 'cord' -- Previously tested:4, 5

% Thin film #2
Lx_2 = Lx_1;              % Length of film #2
Ly_2 = Ly_1;              % Width of film #2
Lz_2 = Lz_1;              % Thickness of film #2
origin_2 = [Lx_1+d,0,0];  % Point of back, left, bottom of cube [x,y,z]
cord_2 = 'min';             % String of coordinate along which direction mesh is to be refined 
mesh_2 = 1;               % Number of subvolumes across the dimension given in the input 'cord' -- Previously tested:4, 5



%****************************END USER INPUTS******************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% Create discretization %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Combine inputs
L_1 = [Lx_1, Ly_1, Lz_1]; % Vector containing length of each side of film #1 [m]
L_2 = [Lx_1, Ly_1, Lz_1]; % Vector containing length of each side of film #1 [m]

% Create discretization for thin film #1
[ R_1, delta_V_1, N_1 ] = thin_film_discretization_slope_t_130_1(L_1, origin_1, mesh_1, cord_1);

% Create discretization for thin film #2
[ R_2, delta_V_2, N_2 ] = thin_film_discretization_t_130_2(L_2, origin_2, mesh_2, cord_2);

% Combine discretization for both films in one matrix
r = [R_1; R_2];

% Total number of subvolumes in both films
N = N_1 + N_2;

% Bulk object start index
ind_bulk = [1, N_1+1];

% Vector of volume of each subvolume [m^3]
delta_V_vector = [delta_V_1, delta_V_2];

% Vector of length of a side of each cubic subvolumes
L_sub = delta_V_vector.^(1/3);


%%%%%%%%%%%%%
% File name %
%%%%%%%%%%%%%

% File name for saved discretizations
description = '2_thin_films';  % Very short description of results
%file_name_saved = [description '_Lx' num2str((1e6)*Lx_1) 'um_Ly' num2str((1e6)*Ly_1) 'um_Lz' num2str((1e9)*Lz_1) 'nm_d' num2str((1e9)*d) 'nm_N'  num2str(N)]; % File name where results will be saved
file_name_saved = [description '_Lx' num2str((1e9)*Lx_1) 'nm_Ly' num2str((1e9)*Ly_1) 'nm_Lz' num2str((1e9)*Lz_1) 'nm_d' num2str((1e9)*d) 'nm_N'  num2str(N)]; % File name where results will be saved


%%%%%%%%%%%%%%%%%%%%%%%
% Plot discretization %
%%%%%%%%%%%%%%%%%%%%%%%

% Figure title
title_string = {['Location of each subvolume for N = ', num2str(N) ' total subvolumes'],['L_x = ' num2str((1e6)*Lx_1) ' \mum, L_y = ' num2str((1e6)*Ly_1) ' \mum, L_z = ' num2str((1e9)*Lz_1) ' nm, d = ' num2str((1e9)*d) ' nm']};

% Visualize discretized lattice
FIG_discretization = figure(1);
plot3(r(:,1), r(:,2), r(:,3), 'x')
if show_titles == 1
    title(title_string, 'fontsize', 14)
end
xlabel('x-axis [m]')
ylabel('y-axis [m]')
zlabel('z-axis [m]')
if show_axes == 0
    grid off
    axis off
end
set(gca, 'Fontsize', 20)
set(gca,'DataAspectRatio',[1 1 1])
view(3)
grid on

% Plot voxel image of discretization
FIG_voxel = figure(2);
[vert, fac] = voxel_image( r, L_sub(1), [], [], [], [], 'on', [], [] );
if show_titles == 1
    title(title_string, 'fontsize', 14)
end
xlabel('x-axis (m)');
ylabel('y-axis (m)');
zlabel('z-axis (m)');
if show_axes == 0
    grid off
    axis off
end
set(gca, 'fontsize', 20)
view(35,20)
%view(5,20)
%view(2)


%%%%%%%%%%%%%%%%
% Save figures %
%%%%%%%%%%%%%%%%

if save_fig == 1
    % Center of subvolume discretization
    fig_path_1 = [saveDir 'figure_voxel/' file_name_saved '_discretization.fig'];
    saveas(FIG_discretization, fig_path_1)
    clear FIG_discretization % Remove previous plot handles

    % Voxel discretization
    fig_path_2 = [saveDir 'figure_voxel/' file_name_saved '_voxel.fig'];
    saveas(FIG_voxel, fig_path_2)
    clear FIG_voxel % Remove previous plot handles
end


%%%%%%%%%%%%%%%%%%%%%%%%
% Save discretizations %
%%%%%%%%%%%%%%%%%%%%%%%%

% Save discretization matrix to an Excel file
if save_Excel == 1
    disc_path_1 = [saveDir 'Excel_files/' file_name_saved '_discretization.xlsx'];
    writematrix(r, disc_path_1)
end

% Save discretization matrix to a .txt file
if save_txt == 1
    disc_path_2 = [saveDir file_name_saved '_discretization.txt'];
    writematrix(r, disc_path_2,'Delimiter','tab')
end
