%
% Rectangular Transducer Example 1.
%

% $Revision: 565 $ $Date: 2009-09-17 22:24:06 +0200 (Thu, 17 Sep 2009) $ $LastChangedBy: dream $

Fs = 10;   % Sampling freq. in MHz.
Ts = 1/Fs;

%
% Observation point(s).
%
z = 10; % [mm]

% 1 : one point, 0 : 50 points along x-axis.
if 0
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
dx = 0.03;                              % [mm].
dy = 0.03;                              % [mm]
dt = Ts;                                % [us].
nt = 400;                               % Length of spatial impulse response vector.
s_par = [dx dy dt nt];

% Material parameters.
v     = 1.0;                            % Normal velocity.
cp    = 1500;                           % Sound speed.
alfa  = 0;                              % Absorbtion [dB/(cm MHz)].
m_par = [v cp alfa];

% Delay.
t_z = z*1e3/cp;
%delay = 0;                             % Start at 0 [us].
delay = t_z;                            % Start at t_z [us].

% Geometrical parameters.
a = 10;					% x-size [mm].
b = 15;					% y-size [mm].
geom_par = [a b];

H = dreamrect(Ro,geom_par,s_par,delay,m_par,'stop');


figure(1)
clf

t = 0:Ts:Ts*(nt-1);

if size(H,2)>1
  mesh(xo,t,H);
  xlabel('x [mm]')
  ylabel('t [\mus]')
  axis tight
  view([135 32])
else
  plot(t,H);
  axis tight
  ax = axis;
  %axis([0 50 ax(3) ax(4)])
  xlabel('t [\mus]')
  grid
end
title('Spatial Impulse Response for a Rectangular Transducer')
