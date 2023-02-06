function [ R, delta_V, N_disc ] = thin_film_discretization( L, origin, mesh, cord )

% This function makes a cubic lattice that consists of N_disc points.  
% The lattice locations along with the size of each subvolume is output.


% INPUTS:  L = vector containing length of each side of total cube [m]
%            = [L_x, L_y, L_z]
%               L_x = length of side in x-direction
%               L_y = length of side in y-direction
%               L_z = length of side in z-direction
%          origin = point of back, bottom, left of cube [x,y,z]
%          mesh = number of subvolumes across the dimension given in the input 'cord'
%                 (i.e., degree of mesh refinement (1-10), where low values correspond 
%                 to coarse meshing and large values correspond to more refined meshing)
%          cord = string of coordinate along which direction mesh is to be refined 
%                 (Takes inputs 'x', 'y', 'z', or 'min')
%                 'x' = x-coordinate direction
%                 'y' = y-coordinate direction
%                 'z' = z-coordinate direction
%                 'min' = coordinate direction of minimum length min(L)
%
%
% OUTPUTS: R = cubic lattice (matrix of size N_disc x 3)
%          delta_V = (1 x N_disc) vector of size of discretized subvolumes [m^3]
%          N_disc = total number of discretized subvolumes


% Preallocate vector of number of discretizations for each side
N = zeros(1,3); 

% if strcmp(cord, 'min')           % Discretize along minimum side length.
%     L_min = min(L);              % Find length of smallest side
%     delta_L = L_min/mesh;        % Size of discretized length
%     N = round(L./delta_L);       % Vector of number of discretizations for each side
% elseif strcmp(cord, 'x')         % Discretize along x-direction.
%     delta_L = L(1)/mesh;         % Size of discretized length
%     N(1) = round(L(1)./delta_L); % Vector of number of discretizations for each side
%     N(2) = round(L(2)./delta_L); % Vector of number of discretizations for each side
%     N(3) = round(L(3)./delta_L); % Vector of number of discretizations for each side
% elseif strcmp(cord, 'y')         % Discretize along y-direction.
%     delta_L = L(2)/mesh;         % Size of discretized length
%     N(1) = round(L(1)./delta_L); % Vector of number of discretizations for each side
%     N(2) = round(L(2)./delta_L); % Vector of number of discretizations for each side
%     N(3) = round(L(3)./delta_L); % Vector of number of discretizations for each side
% elseif strcmp(cord, 'z')         % Discretize along z-direction.
%     delta_L = L(3)/mesh;         % Size of discretized length
%     N(1) = round(L(1)./delta_L); % Vector of number of discretizations for each side
%     N(2) = round(L(2)./delta_L); % Vector of number of discretizations for each side
%     N(3) = round(L(3)./delta_L); % Vector of number of discretizations for each side
% end

if strcmp(cord, 'min')           % Discretize along minimum side length.
    L_min = min(L);              % Find length of smallest side
    delta_L = L_min/mesh;        % Size of discretized length
elseif strcmp(cord, 'x')         % Discretize along x-direction.
    delta_L = L(1)/mesh;         % Size of discretized length
    
elseif strcmp(cord, 'y')         % Discretize along y-direction.
    delta_L = L(2)/mesh;         % Size of discretized length
    
elseif strcmp(cord, 'z')         % Discretize along z-direction.
    delta_L = L(3)/mesh;         % Size of discretized length
  
end
    N = round(L./delta_L);       % Vector of number of discretizations for each side




% Cut cube into subvolumes
x1 = linspace(origin(1), origin(1)+L(1), N(1)+1);  % Location of discretized edges along x direction
y1 = linspace(origin(2), origin(2)+L(2), N(2)+1);  % Location of discretized edges along y direction
z1 = linspace(origin(3), origin(3)+L(3), N(3)+1);  % Location of discretized edges along z direction

% Total number of subvolumes in rectangular prism
N_disc = N(1)*N(2)*N(3);

% Discretized volume (volume of a single subvolume)
delta_V = (delta_L^3)*ones(1,N_disc);

% Make discretized lattice for rectangular prism
R = zeros(N_disc, 3);  % Preallocate location vector matrix
n = 1;                 % Set number of iterations to 1
for k = 1:N(3)         % Loop through discretized z axis
    for j = 1:N(2)     % Loop through discretized y axis
        for h = 1:N(1) % Loop through discretized x axis
            x = x1(h) + delta_L/2; % Location of discretized volume centers along x direction
            y = y1(j) + delta_L/2; % Location of discretized volume centers along y direction
            z = z1(k) + delta_L/2; % Location of discretized volume centers along z direction
            R(n,1) = x; % Populate location matrix
            R(n,2) = y; % Populate location matrix
            R(n,3) = z; % Populate location matrix
            n = n+1;    % Increase count of loop iteration
        end
    end
end







% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Lowest common multiple approach (exact volume and geometry maintained) %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Find length of smallest side
% L_min = min(L);
% 
% % Find greatest common divisor of all side lengths
% N_smallest_side = mesh*double(lcm(sym(L)));
% 
% 
% % Size of discretized length
% delta_L = L_min/N_smallest_side;
% 
% % Vector of number of discretizations for each side
% N = round(L./delta_L);
% 
% % Cut cube into subvolumes
% x1 = linspace(origin(1), origin(1)+L(1), N(1)+1);  % Location of discretized edges along x direction
% y1 = linspace(origin(2), origin(2)+L(2), N(2)+1);  % Location of discretized edges along y direction
% z1 = linspace(origin(3), origin(2)+L(3), N(3)+1);  % Location of discretized edges along z direction
% 
% % Discretized volume (volume of a single subvolume)
% delta_V = delta_L^3;
% 
% % Total number of subvolumes in rectangular prism
% N_disc = N(1)*N(2)*N(3);
% 
% % Make discretized lattice for rectangular prism
% R = zeros(N_disc, 3);  % Preallocate location vector matrix
% n = 1;                 % Set number of iterations to 1
% for k = 1:N(3)         % Loop through discretized z axis
%     for j = 1:N(2)     % Loop through discretized y axis
%         for h = 1:N(1) % Loop through discretized x axis
%             x = x1(h) + delta_L/2; % Location of discretized volume centers along x direction
%             y = y1(j) + delta_L/2; % Location of discretized volume centers along y direction
%             z = z1(k) + delta_L/2; % Location of discretized volume centers along z direction
%             R(n,1) = x; % Populate location matrix
%             R(n,2) = y; % Populate location matrix
%             R(n,3) = z; % Populate location matrix
%             n = n+1;    % Increase count of loop iteration
%         end
%     end
% end