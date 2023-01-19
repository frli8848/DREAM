%%
%%
%%

Fs = 4;   % Sampling freq. [MHz].
Ts = 1/Fs; % [us].

%% NB. The number of obs points must be >= 64 (which is the OpenCL
%% work group size) and be of size N * 64 (ie., divisible by 64)!

z = 10.0;

%% Points along x-axis.
d  = 50.0/255;
xo = (-25:d:25);
yo = zeros(length(xo),1);
zo = z*ones(length(xo),1);
Ro = [xo(:) yo(:) zo(:)];

%% Descretization parameters.
dx = 0.05;                % [mm].
dy = 0.05;                % [mm]
dt = Ts;                  % [us].
nt = 300;                 % Length of spatial impulse response vector.
s_par = [dx dy dt nt];

t = 0:Ts:Ts*(nt-1);

%% Material parameters.
v     = 1.0;                    % Normal velocity.
cp    = 1000;                   % Sound speed.
alpha  = 0.0;                   % Absorbtion (dB/cm Hz).
m_par = [v cp alpha];

t_z = z*1e-3/cp * 1e6;          % us

delay = 0;


%% dreamrect
a = 10;				% x-size.
b = 15;				% y-size.
geom_par = [a b];
Hrect = dreamrect(Ro,geom_par,s_par,delay,m_par,'stop','gpu');

%% dreamcirc
R = 10;				% Radius.
geom_par = [R];
Hrect = dreamcirc(Ro,geom_par,s_par,delay,m_par,'stop','gpu');

%% rect_sir
a = 10;				% x-size.
b = 15;				% y-size.
geom_par = [a b];
Hrect_sir = rect_sir(Ro,geom_par,s_par(3:4),delay,m_par(1:2),'gpu');
