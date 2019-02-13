%
% Example: The Synthetic Aperture Focusing Technique (SAFT).
%


Fs = 10;   % Sampling freq. in MHz.
Ts = 1/Fs;

n_cpus = 2;

z = 50; % Point target depth [mm].

d  = 0.125;                             % Array pitch [mm].
xo = (-50:d:50);                        % Horizontal scanning
                                        % distance  [mm].
yo = zeros(length(xo),1);
zo = z*ones(length(xo),1);
Ro = [xo(:) yo(:) zo(:)];


% Descretization parameters.
dx = 0.03;                              % [mm].
dy = 0.03;                              % [mm]
dt = Ts;                                % [us].
nt = 700;                               % Length of spatial impulse response vector.
s_par = [dx dy dt nt];

% Material parameters.
v     = 1.0;                            % Normal velocity.
cp    = 1500;                           % Sound speed.
alfa  = 0;                              % Absorbtion [dB/(cm MHz)].
m_par = [v cp alfa];

% Delay.
delay = 0;                              % Start at 0 [us].
%t_z = z*1e3/cp;
%delay = t_z;                           % Start at t_z [us].

% Geometrical parameters.
r = 1.0;                                % Radius [mm].
geom_par = [r];

H = dreamcirc(Ro,geom_par,s_par,delay,m_par,'stop');

%
% Simulate a 3.5 MHz transducer,
%
he_nt = 50;                             % Length of the electrical impulse response.
t = (0:((he_nt-1)))*Ts;

f0 = 3.5;                               % Center frequency [MHz].
fprintf('\nf = %1.2f [MHz]\n',f0);
lambda = cp/f0/1e3;                     % [mm].
fprintf('lambda = %1.2f [mm]\n',lambda);
t0 = 0.55;                              % Time delay to max amplitude [mus].
a_n = 20;                               % Envelop parameter.
h_e = -exp(-a_n.*(t-t0).^2).*cos(2.*pi.*f0.*t);
%h_e = h_e - mean(h_e);

P = fftconv_p(H,h_e,n_cpus);            % Electrical immpulse response.
B = fftconv_p(P,P,n_cpus);              % Double-path response.

figure(1)
clf

N = size(B,1);
t = 0:Ts:Ts*(N-1);

imagesc(xo,t,abs(B));
xlabel('x [mm]','FontSize',16)
ylabel('t [\mus]','FontSize',16)
axis([-50 50 62 76])
title('B-scan data','FontSize',16);


%
% Positions of the scanned transducer.
%
xo = (-50:d:50);                        % -50 -> 50 mm.
yo = zeros(length(xo),1);
zo = zeros(length(xo),1);
To = [xo(:) yo(:) zo(:)];

%
% Observation grid matrix.
%
x = -50:0.25:50;
z = 40:.25:60;
N = length(x);
M = length(z);

[X,Z] = meshgrid(x,z);
Y = zeros(size(X));
xo = X(:);
yo = Y(:);
zo = Z(:);
Ro2 = [xo yo zo];

%
% SAFT processing.
%
s_par = dt;
m_par = cp;
a = 50; % Synthetic aperture.

if (n_cpus == 1)
  y = saft(B,To,delay,s_par,m_par,Ro2,a);
else
  y = saft_p(B,To,delay,s_par,m_par,Ro2,a,n_cpus);
end

% Rshape the vector y to an MxN image.
Y = reshape(y,M,N);

figure(2)
clf
imagesc(x,z,abs(Y));
xlabel('x [mm]','FontSize',16)
ylabel('z [mm]','FontSize',16)
title('SAFT Processed Image','FontSize',16);
