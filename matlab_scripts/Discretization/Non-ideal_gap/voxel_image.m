function [vert, fac] = voxel_image( pts, vox_sz, color, alpha, edgec, w_line, light_option, vox_mag, c_limits )
%VOXEL_IMAGE Creates a 3D voxel image
%   Parameters:
%   pts    - n x 3 matrix with 3D points
%   vox_sz - 1 x 3 vector with voxel size - if vox_sz is a scalar, all
%            edges will have the same length
%   color  - face color
%   alpha  - face alpha (opacity)
%   edgec  - edge color
%   w_line - linewidth in plot
%   light_option - "on" or "off" for external lighting or "heatmap" for a heatmap
%   vox_mag - assigned magnitude to each voxel
%   c_limits - limits of vox_mag values to plot in colorbar and heat map
%
%   Return values:
%   vert   - 8n x 3 matrix containing the vertices of the voxels
%   fac    - 6n x 4 matrixes containing the indexes of vert that form a
%            face

% Example:
%
% pts = [1,1,1; 2,2,2];
% vs = [0.5, 0.5, 0.5];
% voxel_image(pts, vs);
% view([-37.5, 30]);
%
% This example creates two voxels (at (1,1,1) and (2,2,2)) with edges of
% length 0.5
%
% Author: Stefan Schalk, 11 Feb 2011


% Number of points
np = size(pts,1); 

% Extract colors
CM = get(gca,'ColorOrder');

if (nargin < 1)
    error('No input arguments given');
end
if (nargin < 2)
    vox_sz = [1,1,1];
end
if (nargin < 3)  || isempty(color)
    %color = 'b';
    %color = CM(1,:);
    %color = rgb('navy');
    %color = rgb('peru');
    %color = rgb('royalblue');
    color = rgb('red');
end
if (nargin < 4) || isempty(alpha)
    alpha = 1;
end
if (nargin < 5) || isempty(edgec)
    %edgec = 'k';
    %edgec = CM(6,:);
    %edgec = rgb('royalblue');
    %edgec = rgb('peachpuff');
    %edgec = rgb('whitesmoke');
    edgec = rgb('black');
end
if (nargin < 6) || isempty(w_line)
    w_line = 0.25;
end
if (nargin < 7) || isempty(light_option)
    light_option = 'on';
end
if (nargin < 8)  || isempty(vox_mag)
    vox_mag = ones(np,1);
end
if (nargin < 9)  || isempty(c_limits)
    c_limits = [min(vox_mag), max(vox_mag)];
end


if (size(pts,2) ~= 3)
    error('pts should be an n x 3 matrix');
end
if (isscalar(vox_sz))
    vox_sz = vox_sz*ones(1,3);
end
if (size(vox_sz,1) ~= 1 || size(vox_sz,2) ~= 3)
    error('vox_sz should be an 1 x 3 vector');
end


%%%%%%%%%%%%%
% Color map %
%%%%%%%%%%%%%

if strcmp(light_option, 'heatmap')
    % Convert vector of values of heat dissipated in each subvolume to an RGB colormap value
    color_vox_mag = vals2colormap(vox_mag, 'fireice', c_limits);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%



vert = zeros(8*np,3);
%vert_color = zeros(8*np,3);
fac = zeros(6*np,4,'uint32');
vert_bas = [...
    -0.5,-0.5,-0.5;
    0.5,-0.5,-0.5;
    0.5,0.5,-0.5;
    -0.5,0.5,-0.5;
    -0.5,-0.5,0.5;
    0.5,-0.5,0.5;
    0.5,0.5,0.5;
    -0.5,0.5,0.5];
vert_bas = vert_bas.*([vox_sz(1).*ones(8,1), vox_sz(2).*ones(8,1), vox_sz(3).*ones(8,1)]);
fac_bas = [...
    1,2,3,4;
    1,2,6,5;
    2,3,7,6;
    3,4,8,7;
    4,1,5,8;
    5,6,7,8];



face_color = zeros(6*np,3); % Preallocate
for vx = 1:np % Loop through all points
    a = ((vx-1)*8+1):vx*8;
    for dim = 1:3
        vert( a,dim ) = vert_bas(:,dim) + pts(vx,dim);
    end
    fac ( ((vx-1)*6+1):vx*6,: ) = (vx - 1)*8*ones(6,4) + fac_bas;
    if strcmp(light_option, 'heatmap')
        face_color( ((vx-1)*6+1):vx*6,: ) = repmat(color_vox_mag(vx,:), [6,1]); % Set voxel face color based on heat dissipated in that subvolume
    end
end






%%%%%%%%%%%%%%%
% Create plot %
%%%%%%%%%%%%%%%

if strcmp(light_option, 'on')
    h = patch('Vertices',vert,'Faces',fac,'FaceColor',color,'FaceAlpha',alpha,'Edgecolor',edgec, 'linewidth', w_line,...
        'FaceLighting', 'gouraud', 'EdgeLighting', 'gouraud', 'AmbientStrength', 0.45, 'SpecularStrength', 1);
    
    N_lights = 8; % Max at 8
    x_light = -100*abs(min(pts(:,1)));
    y_light = -500*abs(min(pts(:,2)));
    z_light = 0*abs(min(pts(:,3)));
    for ii = 1:N_lights
        %light('Position',[-1 -5 0],'Style','local')
        light('Position',[x_light y_light z_light],'Style','local')
    end
    
elseif strcmp(light_option, 'off')
    h = patch('Vertices',vert,'Faces',fac,'FaceColor',color,'FaceAlpha',alpha,'Edgecolor',edgec, 'linewidth', w_line);
    
elseif strcmp(light_option, 'heatmap')
    h = patch('Vertices',vert,'Faces',fac,'FaceColor','flat','FaceAlpha',alpha,'Edgecolor',edgec, 'linewidth', w_line,...
        'FaceVertexCData', face_color, 'CDataMapping', 'direct');
    CB = colorbar; % Define colorbar as an object
    CB.Label.String = 'Heat dissipated in a subvolume, Q [W]'; % Label colorbar
    caxis(c_limits)  % Set colorbar limits
    colormap fireice % Set colormap color scheme
end

%h = patch('Vertices',vert,'Faces',fac, 'FaceLighting', 'gouraud', 'FaceColor', 'interp');
view(3)
axis equal
%title('Objects discretized into subvolumes on a cubic lattice')
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
set(gca, 'fontsize', 20)
grid on










end

