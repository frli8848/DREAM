
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


%
% 2D Rectangular Array Transducer Example 1.
%

Fs = 10;   % Sampling freq. in MHz.
Ts = 1/Fs;

% Snapshot point
z = 50;
dz = 2; % Average field over dz mm.

%
% Observation point(s).
%

d  = 0.5;                               % [mm]

% X-Z plane;
[xo,zo] = meshgrid(-50:d:50,(0:d:80));
yo = zeros(size(xo));

% Y-Z plane;
[yo2,zo2] = meshgrid(-50:d:50,(0:d:80));
xo2 = zeros(size(yo2));

Ro = [[xo(:)' xo2(:)']'  [yo(:)' yo2(:)']' [zo(:)' zo2(:)']'];

fprintf('Number of observation points = %d\n\n',size(Ro,1));

% Material parameters.
v     = 1.0;                            % Normal velocity.
cp    = 1500;                           % Sound speed.
alfa  = 0;                              % Absorbtion [dB/(cm MHz)].
m_par = [v cp alfa];

% Delay.
t_z = (z-dz/2)*1e3/cp;
delay = t_z;                            % Start at t_z [us].
%delay = 0;                             % Start at 0 [us].


% Descretization parameters.
dx = 0.03;                              % [mm].
dy = 0.03;                              % [mm]
dt = Ts;                                % [us].
nt = round(dz*1e3/cp/Ts);               % Length of spatial impulse response vector.
s_par = [dx dy dt nt];

% Element size [mm].
a = 0.9;
%b = 10; % For 1D array.
b = 0.9; % For 2D array.
geom_par = [a b];

%
% Grid function (position vectors of the elements).
%

x = -10:1:10;


% 1D array
%gx = x(:);
%gy = zeros(length(gx),1);
%gz = zeros(length(gx),1);

% 2D array
[gx,gy] = meshgrid(x);
gx = gx(:);
gy = gy(:);
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
%steer_met = 'off';
steer_met = 'x';
%steer_met = 'y';
%steer_met = 'xy';

theta  = -25;				% Angle in x-direction.
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

fprintf('Computing array responses...');
H = dream_arr_rect_p(Ro,geom_par,G,s_par,delay,m_par,foc_met,...
                     focal,steer_met,steer_par,apod_met,apod,win_par,n_cpus,'ignore');
fprintf('done\n');

figure(1);
clf;

%title(['Spatial Impulse Response for a 2D Array Transducer with Rectangular ' ...
%    'Elements'])



%
% Simulated electrical impulse response of the array elements.
%

he_nt=50;
t = (0:((he_nt-1)))*Ts;

f0 = 3.5;                             % Center frequency [MHz].
t0 = 0.55;                            % Time delay to max amplitude [mus].
a_n = 20;                             % Envelop parameter.

h_e = -exp(-a_n.*(t-t0).^2).*cos(2.*pi.*f0.*t);

fprintf('\nf = %1.2f [MHz]\n',f0);
lambda = cp/f0/1e3; % [mm].
fprintf('lambda = %1.2f [mm]\n',lambda);

h_e = h_e - mean(h_e);

% Number of samples to max amplitude of pulse.
[mx,p_delay] = max(h_e);

f_tmp = abs(freqz(h_e,1,1024));
h_e = h_e/max(f_tmp); % Unity gain at center freq.
f_e = abs(freqz(h_e,1,1024));
f = (0:1023)/1045/Ts/2;


%
% Input signal.
%

u = 1; % A unit pulse.
%u = ones(20,1); % A square pulse.


%
% Simulated pressure response.
%

h = conv(h_e,u);

n_cpus = 1;
Ht = fftconv_p(H,h,n_cpus);


%
% Plot Wavefield snaphot.
%

figure(1);
clf;

h_z = sum(abs(hilbert(Ht)));

n_x = size(xo,1)*size(xo,2);
%n_y = size(yo,1)*size(yo,2);

Hz_x = reshape(h_z(1:n_x),size(xo,1),size(xo,2));
Hz_y = reshape(h_z(n_x+1:end),size(yo,1),size(yo,2));

[I,J] = find(Hz_x == 0);
Hz_x2 = Hz_x;
for n=1:length(I)
  Hz_x2(I(n),J(n)) = nan;
end
%Hz_x2(I,J) = Hz_x(I,J);

[I,J] = find(Hz_y == 0);
Hz_y2 = Hz_y;
for n=1:length(I)
  Hz_y2(I(n),J(n)) = nan;
end
%Hz_y2(I,J) = Hz_y(I,J);

h1 = surf(xo,yo,-zo,Hz_x2);
%h1 = surf(xo,yo,-zo,Hz_x);
set(h1,'EdgeAlpha',0);
set(h1,'EdgeColor','none');

hold on;
h2 = surf(yo,xo,-zo,Hz_y2);
%h2 = surf(yo,xo,-zo,Hz_y);
set(h2,'EdgeAlpha',0);
set(h2,'EdgeColor','none');

xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');

view([50 40]);

title('Wavefiled snapshot for a 2D Array Transducer with Rectangular Elements')
