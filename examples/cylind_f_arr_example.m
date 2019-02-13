%
% Cylindrical Focused Array Transducer Example 1.
%

% $Revision: 565 $ $Date: 2009-09-17 22:24:06 +0200 (Thu, 17 Sep 2009) $ $LastChangedBy: dream $

Fs = 10;   % Sampling freq. in MHz.
Ts = 1/Fs;

%
% Observation point(s).
%
z = 50; % [mm]

% 1 - one point, 0 - 50 points along x-axis.
if 1
  % One point.
  xo = 0;
  yo = 0;
  zo = z;
  ro = [xo yo zo];
  Ro = ro;
  disp(['Observation point (x,y,z) = ' num2str(Ro)])
else
  %  Points along x-axis.
  d  = 1;                               % [mm]
  xo = (0:d:50);                        % 0-50 mm.
  yo = zeros(length(xo),1);
  zo = z*ones(length(xo),1);
  Ro = [xo(:) yo(:) zo(:)];
end

% Descretization parameters.
dx = 0.01;                              % [mm].
dy = 0.01;                              % [mm]
dt = Ts;                                % [us].
nt = 500;                               % Length of spatial impulse response vector.
s_par = [dx dy dt nt];

% Material parameters.
v     = 1.0;                            % Normal velocity.
cp    = 1000;                           % Sound speed [m/s].
alfa  = 0;                              % Absorbtion [dB/(cm MHz)].
m_par = [v cp alfa];

% Delay.
t_z = z*1e3/cp;
%delay = 0;                             % Start at 0 [us].
delay = t_z;                            % Start at t_z [us].

% Element size [mm].
a = 1;                                  % x-size.
b = 20;                                 % y-size.
R = 100;				% Curvature radius.
geom_par = [a b R];

% Grid function (position vectors of the elements).
x = -10:1:10;
[gx,gy] = meshgrid(x);
gx = gx(:);
gy = zeros(length(gx),1);
gz = zeros(length(gx),1);
G = [gx gy gz];

%
% Focusing parameters.
%
%foc_met = 'off';
%foc_met = 'x';
%foc_met = 'y';
foc_met = 'xy';
%foc_met = 'x+y';
focal = 50;                             % Focus radius [mm].

% User defined focusing.
%foc_met = 'ud';
%
% Focusing vector for 'ud'. Delays in [us].
%focal = zeros(length(gx),1);           % unfocused.

% Beam steering.
steer_met = 'off';
%steer_met = 'x';
%steer_met = 'y';
%steer_met = 'xy';

theta  = 0;				% Angle in x-direction.
phi    = 0;                             % Angle in y-direction.
steer_par = [theta phi];

% Apodization.
apod_met = 'off';
%apod_met = 'ud';                       % User defined.
%apod_met = 'triangle';
%apod_met = 'gauss';
%apod_met = 'raised';                   % Raised cosine.
%apod_met = 'simply';                   % Simply supported.
%apod_met = 'clamped';
apod = ones(length(gx),1);              % Apodization weights for 'ud'.
win_par = 1;                            % Parameter for raised cos and Gaussian apodization functions.

H = dream_arr_cylind_f(Ro,geom_par,G,s_par,delay,m_par,...
    foc_met,focal,steer_met,steer_par,apod_met,apod,win_par,'stop');


figure(1)
clf

t = 0:Ts:Ts*(nt-1);

if size(H,2)>1
  mesh(xo,t+delay,H);
  xlabel('x [mm]')
  ylabel('t [\mus]')
  axis tight
  view([135 32])
else
  plot(t+delay,H);
  axis tight
  ax = axis;
  %axis([0 50 ax(3) ax(4)])
  xlabel('t [\mus]')
  grid
end
title(['SIR for a Linear Array Transducer with Cylindrical ' ...
      'Elements'])
