%
% Setup array and sampling grid.
%



if ~exist('foc_dep')
  disp('Using foc_dep = 50 mm'); 
  foc_dep = 50;                            % Focus depth.
end

fprintf('Array focused at  = %1.1f [mm]\n',foc_dep);


% Focusing parameters.
if ~exist('theta')
  disp('Using theta = 0 ');
  theta = 0; 			        % Steering angle.
end

fprintf('Array steered = %1.1f [deg]\n',theta);

% Soundspeed in water (in [m/s]).
cp = 1500;

% Length of A-scans.
K = 1000;

%
% Define the array.
%

% Spacing between array elements.
%arr_samp = 0.25; % [mm] 
arr_samp = 1; % [mm] 

a = 0.9;
%a = 0.9/2;
%a = 0.4;
%a = 0.2;

% Transmit.
%nt = 2*8;
nt = 8;
at = a; % [mm] Transmit element size(s).
xt =  arr_samp*(-nt+0.5:nt-0.5);
trans = [xt; at*ones(1,2*nt);];

fprintf('Array aperture = %1.1f [mm]\n',xt(length(xt))-xt(1));

% Receive.
%nr = 2*8;
nr = 8;
ar = a;  % [mm] Receive element size(s).
xr = arr_samp*(-nr+0.5:nr-0.5);
rec = [xr; ar*ones(1,2*nr);];

% Number of transducer elements.
L = size(xt,2);

% Load data just to get us vector.
%dp = '/disk2/users/staff/fl/projects/articles/opt_foc/data/pulse_echo_data';
%[D,xd,yd,us] = ndt2mat([dp '/16foc_st20_1_2_scatt_50mm_32av.ndt']);    
%Ts = us(2) - us(1);
%idx_z_min = 150;
%idx_z_max = 450;

%
% For simulations.
%
Ts = 1/25; % Sampling freq. in [us].
us = (0:Ts:Ts*40000);


%
% Spatial sampling stuff.
%
if strcmp(grid_size,'large')
  x_max = 40;
  x_max = 25;
  
  %idx_z_min = 1501;
  %idx_z_max = 2001-166;
  idx_z_min = 1;
  idx_z_max = 2001-166;
  
  %spat_samp = 0.125; % 1 mm.
  %spat_samp_z = 0.05; % Samples / mm.
  spat_samp = 0.3; % 1 mm.
  spat_samp_z = 0.3; % Samples / mm.
end
  
if strcmp(grid_size,'normal')
  x_max = 25;
  
  idx_z_min = 1501;
  idx_z_max = 2001-166;
  
  spat_samp = 1; % 1 mm.
  spat_samp_z = 1/10; % Samples / mm.
end

if strcmp(grid_size,'all')
  x_max = 25;
  
  idx_z_min = 1;
  idx_z_max = 2678; % 80 mm
  
  spat_samp = 0.5; % 1 mm.
  spat_samp_z = 1; % Samples / mm.
end

%
% Horizontal grid.
%
% Spatial sampling (in mm).
x = -x_max:spat_samp:x_max;
N = length(x)

%
% Vertical grid size.
%
% spat_samp_z = Ts*cp/1e6; 		% [m].
z_min = us(idx_z_min)*cp/1e3/2
z_max = us(idx_z_max)*cp/1e3/2
z = (z_min:spat_samp_z:z_max);
M = length(z)

%
% Observation grid matrix.
%
[X,Z] = meshgrid(x,z);
Y = zeros(size(X));
xo = X(:);
yo = Y(:);
zo = Z(:);
Ro = [xo yo zo];

%fprintf('z_min = %1.1f\n',z_min);
%fprintf('z_max = %1.1f\n',z_max);

fprintf('KL x MN = %d x %d \n',K*L,M*N);
