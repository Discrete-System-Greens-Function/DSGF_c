

shape_file = ['sphere_subvol_8.xlsx'];
%shape_file = ['2_thin_films_Lx0.2um_Ly1um_Lz200nm_d150nm_N640_discretization.xlsx'];
%shape_file = ['2_thin_films_Lx0.2um_Ly1um_Lz200nm_d150nm_N640_discretization.txt'];
%pow_dissipated_file = ['power_dissipated_subvol.xlsx'];         
pow_dissipated_file = ['power_dissipated_subvol.csv']; %'teste.xlsx'

r = xlsread([shape_file]);  % Complete discretized lattice [m]
%r = readtable([shape_file]);
%Q_total_subvol = readtable([pow_dissipated_file],'Range','B1:B639');
%Q_total_subvol = xlsread([pow_dissipated_file]); %power_dissipated_subvol.xls  teste.xlsx
Q_total_subvol = readtable([pow_dissipated_file]);

delta_V_1=1.25*10^(-22);
delta_V_2=delta_V_1;
N1=length(r);
N2=N1;

% Vector of all subvolume sizes [m^3]
delta_V_vector = [delta_V_1*ones(1,N1), delta_V_2*ones(1,N2)]; % Volume of each subvolume (uniform subvolume size)
L_sub = delta_V_vector.^(1/3);                                 % Length of side of a cubic subvolume


% Set heatmap color axis limits
abs_limit = max(abs(Q_total_subvol)); %max(abs(1.540425E-26)); %max(abs(Q_total_subvol));
c_limits = [-abs_limit, abs_limit];

% Subvolume heat map for full particles (VIEW 1)
FIG_1 = figure(1);
%subplot(1,2,1)
%[vert, fac] = voxel_image( r, L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol.', c_limits );
[vert, fac] = voxel_image( r, L_sub(1), [], [], [], [], 'heatmap',Q_total_subvol, c_limits );

xlabel('x-axis (m)');
ylabel('y-axis (m)');
zlabel('z-axis (m)');
if show_axes == 0
    grid off
    axis off
    colorbar off
end
%view(2)
view(-30,35)