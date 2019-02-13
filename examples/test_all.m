%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Run all transducer functions.
%
% Copyright (C) 2005,2006,2008,2009 Fredrik Lingvall
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% $Revision: 565 $ $Date: 2009-09-17 22:24:06 +0200 (Thu, 17 Sep 2009) $ $LastChangedBy: dream $

ONE_POINT = 1;
%ONE_POINT = 0; %Many points.

% Number of threads to use in the "_p" functions.
% n_cpus = 1;

% Are we running on a Linux machine?
tmp_str = computer;
if size(strfind(tmp_str,'linux')) > 0 | ...
      size(strfind(tmp_str,'GLN')) > 0

  % How many CPUs do we have?
  [dummy,tmp_str]= system('cat /proc/cpuinfo | grep model | grep name');

  n_cpus = size(strfind(tmp_str,'model'),2);

  fprintf('\n*** Detected a %d cpu system ***\n\n',n_cpus);

else
  % Default to one cpu.
  n_cpus = 1;
end



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
cp    = 1000;                   % Sound speed.
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

subplot(2,2,1); % Needed for due to bug in Octave.
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
title('Rectangular transducer')
fprintf('dreamrect\n');

% -------------Focused Rectangular Transducer --------------------------

% Geometrical parameters.
a = 10;				% x-size.
b = 15;				% y-size.
geom_par = [a b];

% Focusing parameters.
foc_met = 'off';
%foc_met = 'x';
%foc_met = 'y';
foc_met = 'xy';
%foc_met = 'x+y';
focal = 100;                            % Focus radius

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

title('Circular transducer')
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

fprintf('dreamcirc\n');

% ------------- Focused Circular Transducer --------------------------

% Geometrical parameters.
a  = 10;				% Radius of the transdecer.
geom_par = [a];

% Focusing parameters.
foc_met = 'off';
%foc_met = 'x';
%foc_met = 'y';
foc_met = 'xy';
%foc_met = 'x+y';
focal = 100;                            % Focus radius

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

% ------------- Focused Cylindrical Transducer --------------------------

% Geometrical parameters.
a = 10;				% x-size of the transducer.
b = 20;				% y-size of the transducer.
R = 100;				% Radius of the curvature.
geom_par = [a b R];

[H,err] = dreamcylind_f(Ro,geom_par,s_par,delay,m_par,'stop');

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
title('Focused cylindrical transducer')
fprintf('dreamcylind_f\n');

figure(2);
clf;

% ------------- Defocused Cylindrical Transducer --------------------------

% Geometrical parameters.
a = 10;				% x-size of the transducer.
b = 20;				% y-size of the transducer.
R = 100;				% Radius of the curvature.
geom_par = [a b R];

[H,err] = dreamcylind_d(Ro,geom_par,s_par,delay,m_par,'stop');

subplot(2,2,1); % Needed for due to bug in Octave.
subplot(3,2,1);

if size(H,2)>1
  mesh(xo,t,H);
  axis('tight');
  view(135,32);
else
  plot(t,H);
  ax = axis;
  %axis([0 50 ax(3) ax(4)])
  %xlabel('t [\mus]')
  grid('on');
end
title('Defocused cylindrical transducer')
fprintf('dreamcylind_d\n');

% ------------- Focused Spherical Transducer --------------------------

% Geometrical parameters.
r = 10;                                 % Radius of the transducer.
R = 100;				% Radius of the curvature.
geom_par = [r R];

[H,err] = dreamsphere_f(Ro,geom_par,s_par,delay,m_par,'stop');

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
fprintf('dreamsphere_f\n');

% ------------- Defocused Spherical Transducer --------------------------

% Geometrical parameters.
r = 10;				% Radius of the transducer.
R = 100;				% Radius of the curvature.
geom_par = [r R];

[H,err] = dreamsphere_d(Ro,geom_par,s_par,delay,m_par,'stop');

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
fprintf('dreamsphere_d_d\n');

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
b = 20;                                 % y-width.
R = 100;				% Radius.
geom_par = [a b R];

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
%apod_met = 'simply';                   % Simply supported.
%apod_met = 'clamped';
apod = ones(length(gx),1);              % Apodization weights for 'ud'.
win_par = 1;                            % Parameter for raised cos and Gaussian apodization functions.

[H,err] = dream_arr_cylind_f(Ro,geom_par,G,s_par,delay,...
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
fprintf('dream_arr_cylind_f\n');

figure(3);
clf

% ------------- Array with (convex) cylindrical elements --------------------------

% Element size [mm].
a = 1;                          % x-width.
b = 20;                                 % y-width.
R = 100;				% Radius.
geom_par = [a b R];

% Grid function (position vectors of the elements).
gx = -10:1:10;
gy = zeros(length(gx),1);
gz = zeros(length(gx),1);
gx=gx(:);  gy=gy(:);  gz=gz(:);
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
%apod_met = 'simply';                   % Simply supported.
%apod_met = 'clamped';
apod = ones(length(gx),1);              % Apodization weights for 'ud'.
win_par = 1;                            % Parameter for raised cos and Gaussian apodization functions.

[H,err] = dream_arr_cylind_d(Ro,geom_par,G,s_par,delay,...
                             m_par,foc_met,focal,steer_met,steer_par,apod_met,apod,win_par,'stop');

subplot(2,2,1); % Needed for due to bug in Octave.
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
fprintf('dream_arr_cylind_d\n');

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Run all transducer functions with parallel computing support.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('Using %d threads.\n',n_cpus);

figure(3);
clf


Fs = 4;   % Sampling freq. [MHz].
Ts = 1/Fs; % [us].

%
% Observation point.
%

z = 10;

if ONE_POINT
  % One point.
  xo = 0;
  yo = 0;
  zo = z;
  ro = [xo yo zo];
  Ro = ro;
  disp(['Observation point (x,y,z) = ' num2str(Ro)])
else
  %  Points along x-axis.
  d  = 10;
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
cp    = 1000;                   % Sound speed.
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

[H,err] = dreamline_p(Ro,geom_par,s_par,delay,m_par,n_cpus,'stop');

subplot(3,2,1)
if size(H,2)>1
  mesh(xo,t,H);
  axis('tight')
  view(135,32)
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  %xlabel('t [\mus]')
  grid('on');
end
title('Line strip transducer')
fprintf('dreamline_p\n');

% ------------- Rectangular Transducer --------------------------

% Geometrical parameters.
a = 10;				% x-size.
b = 15;				% y-size.
geom_par = [a b];

[H,err] = dreamrect_p(Ro,geom_par,s_par,delay,m_par,n_cpus,'stop');

subplot(3,2,2)
if size(H,2)>1
  mesh(xo,t,H);
  axis('tight')
  view(135,32)
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  %xlabel('t [\mus]')
  grid('on');
end
title('Rectangular transducer')
fprintf('dreamrect_p\n');

% -------------Focused Rectangular Transducer --------------------------

% Geometrical parameters.
a = 10;				% x-size.
b = 15;				% y-size.
geom_par = [a b];

% Focusing parameters.
foc_met = 'off';
%foc_met = 'x';
%foc_met = 'y';
foc_met = 'xy';
%foc_met = 'x+y';
focal = 100;                            % Focus radius

[H,err] = dreamrect_f_p(Ro,geom_par,s_par,delay,m_par,foc_met,focal,n_cpus,'stop');

subplot(3,2,3);
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
title('Focused Rectangular transducer')
fprintf('dreamrect_f_p\n');

% ------------- Circular Transducer --------------------------

% Geometrical parameters.
a  = 10;				% Radius of the transdecer.
geom_par = [a];

[H,err] = dreamcirc_p(Ro,geom_par,s_par,delay,m_par,n_cpus,'stop');

subplot(3,2,4);
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
title('Circular transducer')
fprintf('dreamcirc_p\n');

% ------------- Focused Circular Transducer --------------------------

% Geometrical parameters.
a  = 10;				% Radius of the transdecer.
geom_par = [a];

% Focusing parameters.
foc_met = 'off';
%foc_met = 'x';
%foc_met = 'y';
foc_met = 'xy';
%foc_met = 'x+y';
focal = 100;                            % Focus radius

[H,err] = dreamcirc_f_p(Ro,geom_par,s_par,delay,m_par,foc_met,focal,n_cpus,'stop');

subplot(3,2,5);
if size(H,2)>1
  mesh(H);
  axis('tight');
  view(135,32)
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  xlabel('t [\mus]')
  grid('on');
end
title('Focused Circular transducer')
fprintf('dreamcirc_f_p\n');

% ------------- Focused Cylindrical Transducer --------------------------

% Geometrical parameters.
a = 10;				% x-size of the transducer.
b = 20;				% y-size of the transducer.
R = 100;				% Radius of the curvature.
geom_par = [a b R];

[H,err] = dreamcylind_f_p(Ro,geom_par,s_par,delay,m_par,n_cpus,'stop');

subplot(3,2,6);
if size(H,2)>1
  mesh(H);
  axis('tight');
  view(135,32)
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  xlabel('t [\mus]')
  grid('on');
end
title('Focused cylindrical transducer')
fprintf('dreamcylind_f_p\n');

figure(4);
clf

% ------------- Defocused Cylindrical Transducer --------------------------

% Geometrical parameters.
a = 10;				% x-size of the transducer.
b = 20;				% y-size of the transducer.
R = 100;				% Radius of the curvature.
geom_par = [a b R];

[H,err] = dreamcylind_d_p(Ro,geom_par,s_par,delay,m_par,n_cpus,'stop');

subplot(2,2,1); % Needed for due to bug in Octave.
subplot(3,2,1)

if size(H,2)>1
  mesh(xo,t,H);
  axis('tight');
  view(135,32)
else
  plot(t,H);
  ax = axis;
  %axis([0 50 ax(3) ax(4)])
  %xlabel('t [\mus]')
  grid('on');
end
title('Defocused cylindrical transducer')
fprintf('dreamcylind_d_p\n');

% ------------- Focused Spherical Transducer --------------------------

% Geometrical parameters.
r = 10;				% Radius of the transducer.
R = 100;				% Radius of the curvature.
geom_par = [r R];

[H,err] = dreamsphere_f_p(Ro,geom_par,s_par,delay,m_par,n_cpus,'stop');

subplot(3,2,2)
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
title('Focused spherical transducer')
fprintf('dreamsphere_f_p\n');

% ------------- Defocused Spherical Transducer --------------------------

% Geometrical parameters.
r = 10;				% Radius of the transducer.
R = 100;				% Radius of the curvature.
geom_par = [r R];

[H,err] = dreamsphere_d_p(Ro,geom_par,s_par,delay,m_par,n_cpus,'stop');

subplot(3,2,3);
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
title('Defocused spherical transducer')
fprintf('dreamsphere_d_p\n');

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
%apod_met = 'simply';                   % Simply supported.
%apod_met = 'clamped';
apod = ones(length(gx),1);              % Apodization weights for 'ud'.
win_par = 1;                            % Parameter for raised cos and Gaussian apodization functions.

[H,err] = dream_arr_circ_p(Ro,r,G,s_par,delay,m_par,foc_met,focal,...
                           steer_met,steer_par,apod_met,apod,win_par,n_cpus,'stop');

subplot(3,2,4);
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
title('Array with circular elements')
fprintf('dream_arr_circ_p\n');

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
%apod_met = 'simply';                   % Simply supported.
%apod_met = 'clamped';
apod = ones(length(gx),1);              % Apodization weights for 'ud'.
win_par = 1;                            % Parameter for raised cos and Gaussian apodization functions.

[H,err] = dream_arr_rect_p(Ro,geom_par,G,s_par,delay,...
                           m_par,foc_met,focal,steer_met,steer_par,apod_met,apod,win_par,n_cpus,'stop');

subplot(3,2,5);
if size(H,2)>1
  mesh(xo,t,H);
  axis('tight');
  view(135,32)
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  xlabel('t [\mus]')
  grid('on');
end
title('Array with rectangular elements')
fprintf('dream_arr_rect_p\n');

% ------------- Array with (concave) cylindrical elements --------------------------

% Element size [mm].
a = 1;                          % x-width.
b = 20;                                 % y-width.
R = 100;				% Radius.
geom_par = [a b R];

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
%apod_met = 'simply';                   % Simply supported.
%apod_met = 'clamped';
apod = ones(length(gx),1);              % Apodization weights for 'ud'.
win_par = 1;                            % Parameter for raised cos and Gaussian apodization functions.

[H,err] = dream_arr_cylind_f_p(Ro,geom_par,G,s_par,delay,...
                               m_par,foc_met,focal,steer_met,steer_par,apod_met,apod,win_par,n_cpus,'stop');

subplot(3,2,6);
if size(H,2)>1
  mesh(xo,t,H);
  axis('tight');
  view(135,32)
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  xlabel('t [\mus]')
  grid('on');
end
title('Array with cylindrical concave  elements')
fprintf('dream_arr_cylind_f_p\n');

figure(5);
clf

% ------------- Array with (convex) cylindrical elements --------------------------

% Element size [mm].
a = 1;                          % x-width.
b = 20;                                 % y-width.
R = 100;				% Radius.
geom_par = [a b R];

% Grid function (position vectors of the elements).
gx = -10:1:10;
gy = zeros(length(gx),1);
gz = zeros(length(gx),1);
gx=gx(:);  gy=gy(:);  gz=gz(:);
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
%apod_met = 'simply';                   % Simply supported.
%apod_met = 'clamped';
apod = ones(length(gx),1);              % Apodization weights for 'ud'.
win_par = 1;                            % Parameter for raised cos and Gaussian apodization functions.

[H,err] = dream_arr_cylind_d_p(Ro,geom_par,G,s_par,delay,...
                               m_par,foc_met,focal,steer_met,steer_par,apod_met,apod,win_par,n_cpus,'stop');

subplot(2,2,1); % Needed for due to bug in Octave.
subplot(3,2,1);

if size(H,2)>1
  mesh(H);
  axis('tight');
  view(135,32)
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  xlabel('t [\mus]')
  grid('on');
end
title('Array with cylindrical convex elements')
fprintf('dream_arr_cylind_d_p\n');

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
%apod_met = 'simply';                   % Simply supported.
%apod_met = 'clamped';
apod = ones(length(gr),1);              % Apodization weights for 'ud'.
win_par = 1;                            % Parameter for raised cos and Gaussian apodization functions.

[H,err] = dream_arr_annu_p(Ro,G,s_par,delay,...
                           m_par,foc_met,focal,apod_met,apod,win_par,n_cpus,'stop');

subplot(3,2,2);
if size(H,2)>1
  mesh(xo,t,H);
  axis('tight');
  view(135,32)
else
  plot(t,H);
  ax = axis;
  axis([0 50 ax(3) ax(4)])
  xlabel('t [\mus]')
  grid('on');
end
title('Annular array');
fprintf('dream_arr_annu_p\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Misc functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Using %d threads.\n',n_cpus);

L1 = 5000;
L2 = 5000;
N=100;
X = randn(L1,N);
Y = randn(L2,N);

eval('tic')
Z = conv_p(X,Y,n_cpus);
t1= toc;
fprintf('conv_p\n');


eval('tic')
Z2 = fftconv_p(X,Y,n_cpus);
t2= toc;
fprintf('fftconv_p\n\n');

fprintf('Elapsed time conv_p: %f, fftconv_p: %f\n\n',t1,t2);
