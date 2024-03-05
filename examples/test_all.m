%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Run all transducer functions.
%
% Copyright (C) 2005,2006,2008,2009,2019,2021,2024 Fredrik Lingvall
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ONE_POINT = 1;
%ONE_POINT = 0; %Many points.

figure(1);
clf

Fs = 4;   % Sampling freq. [MHz].
Ts = 1/Fs; % [us].

%
% Observation point(s).
%

z = 10;

if ONE_POINT
  % One point.
  xo = 0;
  yo = 0;
  zo = z;
  ro = [xo yo zo];
  Ro = ro;
  disp(['Observation point (x,y,z) = ' num2str(Ro)]);
else
  %  Points along x-axis.
  d  = 1;
  xo = (0:d:50);
  yo = zeros(length(xo),1);
  zo = z*ones(length(xo),1);
  Ro = [xo(:) yo(:) zo(:)];
end

% Descretization parameters.
dx = 0.05;                              % [mm].
dy = 0.05;                              % [mm]
dt = Ts;                                % [us].
nt = 300;                               % Length of spatial impulse response vector.
s_par = [dx dy dt nt];

t = 0:Ts:Ts*(nt-1);

% Material parameters.
v     = 1.0;                            % Normal velocity.
cp    = 1000;                           % Sound speed.
alfa  = 0;                              % Absorbtion (dB/cm Hz).
m_par = [v cp alfa];

t_z = z*1e-3/cp * 1e6;          % us

delay = 0;
%delay = t_z;

% ------------- Line Strip Transducer --------------------------

% Geometrical parameters.
a = 10;				% x-size [mm].
                                %width = dy;				% Strip size (= dy).
geom_par = [a];

[H,err] = dreamline(Ro,geom_par,s_par,delay,m_par,'stop');

subplot(3,2,1);

if size(H,2)>1
  mesh(xo,t,H);
  axis('tight');
  view(135,32);
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)]);
  %xlabel('t [\mus]')
  grid('on');
end
title('Line strip transducer')
fprintf('dreamline\n');

% ------------- Rectangular Transducer --------------------------

% Geometrical parameters.
a = 10;				% x-size.
b = 15;				% y-size.
geom_par = [a b];

[H,err] = dreamrect(Ro,geom_par,s_par,delay,m_par,'stop');

subplot(3,2,2);
if size(H,2)>1
  mesh(xo,t,H);
  axis('tight');
  view(135,32)
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  %xlabel('t [\mus]')
  grid('on');
end
fprintf('dreamrect\n');

% Analytical SIR
[Ha] = rect_sir(Ro,geom_par,s_par(3:4),delay,m_par(1:2));
hold on;
if size(Ha,2)>1
  mesh(xo,t,Ha);
  axis('tight');
  view(135,32)
else
  plot(t,Ha);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  %xlabel('t [\mus]')
  grid('on');
end
fprintf('rect_sir\n');

title('Rectangular transducer')
legend('DREAM','Analytical');


% -------------Focused Rectangular Transducer --------------------------

% Geometrical parameters.
a = 10;				% x-size.
b = 15;				% y-size.
geom_par = [a b];

% Focusing parameters.
%foc_met = 'off';
%foc_met = 'x';
%foc_met = 'y';
foc_met = 'xy';
%foc_met = 'x+y';
focal = 15;                     % Focus distance.

[H,err] = dreamrect_f(Ro,geom_par,s_par,delay,m_par,foc_met,focal,'stop');

subplot(3,2,3);
if size(H,2)>1
  mesh(xo,t,H);
  axis('tight');
  view(135,32);
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  %xlabel('t [\mus]')
  grid('on');
end
title('Focused Rectangular transducer')
fprintf('dreamrect_f\n');

% ------------- Circular Transducer --------------------------

% Geometrical parameters.
a  = 10;				% Radius of the transdecer.
geom_par = [a];

[H,err] = dreamcirc(Ro,geom_par,s_par,delay,m_par,'stop');

subplot(3,2,4);
if size(H,2)>1
  mesh(xo,t,H);
  axis('tight');
  view(135,32);
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  %xlabel('t [\mus]')
  grid('on');
end
title('Circular transducer')

fprintf('dreamcirc\n');

% ------------- Focused Circular Transducer --------------------------

% Geometrical parameters.
a  = 10;                        % Radius of the transdecer.
geom_par = [a];

%% Focusing parameters.
%foc_met = 'off';
%foc_met = 'x';
%foc_met = 'y';
foc_met = 'xy';
%foc_met = 'x+y';
focal = 15;                     % Focus distance.

[H,err] = dreamcirc_f(Ro,geom_par,s_par,delay,m_par,foc_met,focal,'stop');

subplot(3,2,5);
if size(H,2)>1
  mesh(H);
  axis('tight');
  view(135,32);
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  xlabel('t [\mus]')
  grid('on');
end
title('Focused Circular transducer')
fprintf('dreamcirc_f\n');

% ------------- Cylindrical Transducer --------------------------

% Geometrical parameters.
a = 10;				% x-dim size of the transducer.
b = 20;				% y-dim size of the transducer.
Rcurv = 15;                     % Radius of the curvature.
geom_par = [a b Rcurv];

[H,err] = dreamcylind(Ro,geom_par,s_par,delay,m_par,'stop');

subplot(3,2,6);
if size(H,2)>1
  mesh(H);
  axis('tight');
  view(135,32);
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  xlabel('t [\mus]')
  grid('on');
end
title('Cylindrical transducer focused')
fprintf('dreamcylind (concave)\n');

figure(2);
clf;

% Geometrical parameters.
a = 10;				% x-dim size of the transducer.
b = 20;				% y-dim size of the transducer.
Rcurv = -15;                    % Radius of the curvature.
geom_par = [a b Rcurv];

[H,err] = dreamcylind(Ro,geom_par,s_par,delay,m_par,'stop');

subplot(3,2,1);
if size(H,2)>1
  mesh(H);
  axis('tight');
  view(135,32);
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  xlabel('t [\mus]')
  grid('on');
end
title('Cylindrical transducer defocused')
fprintf('dreamcylind (convex)\n');


% ------------- Focused Spherical Transducer --------------------------

% Geometrical parameters.
R = 10;                         % Radius of the transducer.
Rcurv = 25;                     % Radius of the curvature.
geom_par = [R Rcurv];

[H,err] = dreamsphere(Ro,geom_par,s_par,delay,m_par,'stop');

subplot(3,2,2);
if size(H,2)>1
  mesh(xo,t,H);
  axis('tight');
  view(135,32);
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  %xlabel('t [\mus]')
  grid('on');
end
title('Focused spherical transducer')
fprintf('dreamsphere (focused)\n');

% ------------- Defocused Spherical Transducer --------------------------

% Geometrical parameters.
R = 10;              % Radius of the transducer.
Rcurv = -25;         % Radius of the curvature (negative for defocused
geom_par = [R Rcurv];

[H,err] = dreamsphere(Ro,geom_par,s_par,delay,m_par,'stop');

subplot(3,2,3);
if size(H,2)>1
  mesh(xo,t,H);
  axis('tight');
  view(135,32);
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  %xlabel('t [\mus]')
  grid('on');
end
title('Defocused spherical transducer')
fprintf('dreamsphere (defocused)\n');

% ------------- Array with circular elements --------------------------

% Element diameter [mm].
r = 0.5;

% Grid function (position vectors of the elements).
x = -10:1:10;
[gx,gy] = meshgrid(x);
gx = gx(:);
gy = gy(:);
gz = zeros(length(gx),1);
G = [gx gy gz];

% Focusing parameters.
foc_met = 'off';
%foc_met = 'x';
%foc_met = 'y';
%foc_met = 'xy';
%foc_met = 'x+y';
focal = 100;                            % Focus radius

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
%apod_met = 'hann';                     % Hann window.
%apod_met = 'hamming';                  % Hamming window.
%apod_met = 'simply';                   % Simply supported.
%apod_met = 'clamped';
apod = ones(length(gx),1);              % Apodization weights for 'ud'.
win_par = 1;                            % Parameter for raised cos and Gaussian apodization functions.

[H,err] = dream_arr_circ(Ro,r,G,s_par,delay,m_par,foc_met,focal,...
                         steer_met,steer_par,apod_met,apod,win_par,'stop');

subplot(3,2,4);
if size(H,2)>1
  mesh(xo,t,H);
  axis('tight');
  view(135,32);
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  %xlabel('t [\mus]')
  grid('on');
end
title('Array with circular elements')
fprintf('dream_arr_circ\n');

% ------------- Array with rectangular elements --------------------------

% Element size [mm].
a = 1;
b = 1;
geom_par = [a b];

% Grid function (position vectors of the elements).
x = -10:1:10;
[gx,gy] = meshgrid(x);
gx = gx(:);
gy = gy(:);
gz = zeros(length(gx),1);
G = [gx gy gz];

% Focusing parameters.
foc_met = 'off';
%foc_met = 'x';
%foc_met = 'y';
%foc_met = 'xy';
%foc_met = 'x+y';
focal = 100;                            % Focus radius

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
%apod_met = 'hann';                     % Hann window.
%apod_met = 'hamming';                  % Hamming window.
%apod_met = 'simply';                   % Simply supported.
%apod_met = 'clamped';
apod = ones(length(gx),1);              % Apodization weights for 'ud'.
win_par = 1;                            % Parameter for raised cos and Gaussian apodization functions.

[H,err] = dream_arr_rect(Ro,geom_par,G,s_par,delay,...
                         m_par,foc_met,focal,steer_met,steer_par,apod_met,apod,win_par,'stop');

subplot(3,2,5);
if size(H,2)>1
  mesh(xo,t,H);
  axis('tight');
  view(135,32);
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  xlabel('t [\mus]')
  grid('on');
end
title('Array with rectangular elements')
fprintf('dream_arr_rect\n');

% ------------- Array with (concave) cylindrical elements --------------------------

% Element size [mm].
a = 1;                          % x-width.
b = 20;                         % y-width.
Rcurv = 100;                    % Radius.
geom_par = [a b Rcurv];

% Grid function (position vectors of the elements).
gx = -9.5:1:9.5;
gy = zeros(length(gx),1);
gz = zeros(length(gx),1);
gx=gx(:); gy=gy(:); gz=gz(:);
G = [gx gy gz];

% Focusing parameters.
foc_met = 'off';
%foc_met = 'x';
%foc_met = 'y';
%foc_met = 'xy';
%foc_met = 'x+y';
focal = 100;                            % Focus radius

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
%apod_met = 'hann';                     % Hann window.
%apod_met = 'hamming';                  % Hamming window.
%apod_met = 'simply';                   % Simply supported.
%apod_met = 'clamped';
apod = ones(length(gx),1);              % Apodization weights for 'ud'.
win_par = 1;                            % Parameter for raised cos and Gaussian apodization functions.

[H,err] = dream_arr_cylind(Ro,geom_par,G,s_par,delay,...
                           m_par,foc_met,focal,steer_met,steer_par,apod_met,apod,win_par,'stop');

subplot(3,2,6);
if size(H,2)>1
  mesh(xo,t,H);
  axis('tight');
  view(135,32);
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  xlabel('t [\mus]')
  grid('on');
end
title('Array with cylindrical concave  elements')
fprintf('dream_arr_cylind (concave)\n');

figure(3);
clf

% ------------- Array with (convex) cylindrical elements --------------------------

% Element size [mm].
a = 1;                          % x-width.
b = 20;                         % y-width.
Rcurv = -100;                   % Radius.
geom_par = [a b Rcurv];

% Grid function (position vectors of the elements).
gx = -10:1:10;
gy = zeros(length(gx),1);
gz = zeros(length(gx),1);
gx=gx(:); gy=gy(:); gz=gz(:);
G = [gx gy gz];

% Focusing parameters.
foc_met = 'off';
%foc_met = 'x';
%foc_met = 'y';
%foc_met = 'xy';
%foc_met = 'x+y';
focal = 100;                            % Focus radius

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
%apod_met = 'hann';                     % Hann window.
%apod_met = 'hamming';                  % Hamming window.
%apod_met = 'simply';                   % Simply supported.
%apod_met = 'clamped';
apod = ones(length(gx),1);              % Apodization weights for 'ud'.
win_par = 1;                            % Parameter for raised cos and Gaussian apodization functions.

[H,err] = dream_arr_cylind(Ro,geom_par,G,s_par,delay,...
                           m_par,foc_met,focal,steer_met,steer_par,apod_met,apod,win_par,'stop');

subplot(3,2,1);

if size(H,2)>1
  mesh(H);
  axis('tight');
  view(135,32);
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  xlabel('t [\mus]')
  grid('on');
end
title('Array with cylindrical convex elements')
fprintf('dream_arr_cylind (convex)\n');

% ------------- Annualar Array  --------------------------

% Grid function (position vector of the elements).
%gr = 0:1:10;
%gr = 5:5:10;
%gr = [5 5 10];
gr = [5 5 6 6 7 7 8 8 9 9 10];
G = gr(:);

% Focusing parameters.
foc_met = 'on';
%foc_met = 'off';
focal = 100;                            % Focus radius

% Apodization.
apod_met = 'off';
%apod_met = 'ud';                       % User defined.
%apod_met = 'triangle';
%apod_met = 'gauss';
%apod_met = 'raised';                   % Raised cosine.
%apod_met = 'hann';                     % Hann window.
%apod_met = 'hamming';                  % Hamming window.
%apod_met = 'simply';                   % Simply supported.
%apod_met = 'clamped';
apod = ones(length(gr),1);              % Apodization weights for 'ud'.
win_par = 1;                            % Parameter for raised cos and Gaussian apodization functions.

[H,err] = dream_arr_annu(Ro,G,s_par,delay,...
                         m_par,foc_met,focal,apod_met,apod,win_par,'stop');

subplot(3,2,2);
if size(H,2)>1
  mesh(xo,t,H);
  axis('tight');
  view(135,32);
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  xlabel('t [\mus]')
  grid('on');
end
title('Annular array');
fprintf('dream_arr_annu\n');

% ------------- Attenuation response --------------------------

alpha = 3.0;
s_par = [dt nt];
m_par = [cp alpha];
[H] = dream_att(Ro,s_par,delay,m_par);

subplot(3,2,3);
if size(H,2)>1
  mesh(xo,t,H);
  axis('tight');
  view(135,32);
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)]);
  %xlabel('t [\mus]')
  grid('on');
end
title('Attenuation Response')
fprintf('dream_att\n');
