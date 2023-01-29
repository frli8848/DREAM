%%
%% ------------- Delay-and-sum --------------------------
%%

Fs = 50.0;                      % Sampling freq. [MHz].
Ts = 1/Fs;                      % [us].

%% Descretization parameters.
dx = 0.05;                % [mm].
dy = 0.05;                % [mm]
dt = Ts;                  % [us].
nt = 2000;                % Length of spatial impulse response vector.
s_par = [dx dy dt nt];

t = 0:Ts:Ts*(nt-1);

% Material parameters.
v     = 1.0;                    % Normal velocity.
cp    = 1500;                   % Sound speed.
alpha  = 0.0;                   % Absorbtion (dB/cm Hz).
m_par = [v cp alpha];

%% Simulate a single point scatterer at (x=0,z=10)
z_pt = 10;

%%
%% Generate simulated data
%%

%% Simulated electrical impulse response.

nt_he = 150;
t = (0:((nt_he-1)))*Ts;

f0 = 2.5;                             % Center frequency [MHz].
t0 = 0.55;                            % Time delay to max amplitude [us].
a_n = 10;                             % Envelop parameter.

system_delay = t0+0.21; % Delay to the max of the pulse.

h_e = -exp(-a_n.*(t-t0).^2).*cos(2.*pi.*f0.*t);

fprintf('\nf = %1.2f [MHz]\n',f0);
lambda = cp/f0/1e3; % [mm].
fprintf('lambda = %1.2f [mm]\n',lambda);

f_e = abs(freqz(h_e,1,1024));
h_e = h_e/max(f_e); % Unity gain at center freq.

if (exist('DO_PLOTTING'))
  figure(1);
  clf;
  subplot(211)
  plot(t, h_e);
  xlabel('t [{\mu}s]')
  title('System impulse response')

  subplot(212)
  f = (0:1023)/1024/Ts/2;
  plot(f,20*log10(f_e/max(f_e)));
  xlabel('f [MHz]')
end

%%
%% RCA TFM data (full matrix capture - a.k.a FMC)
%%

d  = 0.5;                       % Array pitch
xo = (-25:d:25);
yo = zeros(length(xo),1);
zo = z_pt*ones(length(xo),1);
Ro_t = [xo(:) yo(:) zo(:)];

%% Crossed transmit and receve elemets.

% Geometrical parameters.
a = 0.4;                        % x-size.
b = 50;				% y-size.
geom_par_t = [a b];

delay = 0.0;
[Ht,err] = dreamrect(Ro_t,geom_par_t,s_par,delay,m_par,'stop');

geom_par_r = [b a];
Ro_r = [yo(:) xo(:) zo(:)];
[Hr,err] = dreamrect(Ro_r,geom_par_r,s_par,delay,m_par,'stop');

L = length(xo);
Yfmc = zeros(nt+nt-1+nt_he-1,L^2);

%% Loop over all transmit elements
n_t=1;
for n=1:L:L^2
  Hdp = fftconv_p(Hr,Ht(:,n_t)); % Double-path SIRs for the n_t:th transmit
  Yfmc(:,n:(n+L-1)) = fftconv_p(Hdp,h_e);
  n_t = n_t+1;
end

Yfmc = Yfmc/max(max(abs(Yfmc))); % Normalize amplitudes

num_elements = size(xo,2);
Gt = [xo(:) zeros(num_elements,1) zeros(num_elements,1)];
Gr = [zeros(num_elements,1) xo(:) zeros(num_elements,1)];


%%
%% 3D Observation points for DAS RCA
%%

x = -25:0.5:25;
y = -25:0.5:25;
Nz = 64*4; % Make sure its a factor of 64 (the OpenCL work group size).
z = (0:(Nz-1))/Nz*20;
[X,Y,Z] = meshgrid(x,y,z);
Ro_tfm = [X(:) Y(:) Z(:)];

Nx = length(x);
Ny = length(y);

delay = system_delay; % Compensate for the pulse/system (transducer) delay.
tic
Im_tfm = das(Yfmc, Gt, Gr, Ro_tfm, dt, delay, cp,'rca');
toc;

if (exist('DO_PLOTTING'))

  figure(4);
  clf;

  O_cpu = reshape(Im_tfm, Nx*Ny, Nz)';
  c_scan_cpu = reshape(max(abs(O_cpu)), Nx, Ny);
  mx = max(max(c_scan_cpu));

  imagesc(x, y, 20.0*log10(c_scan_cpu/mx));
  h_cb = colorbar;
  h_cb_title = get(h_cb,'Title');
  set(h_cb_title,'String','Normalized Amplitude [dB]')
  axis square;
  xlabel('x [mm]');
  ylabel('y [mm]');
  title('C-scan RCA DAS beamformed data');
end

tic
Im_tfm_gpu = das(Yfmc, Gt, Gr, Ro_tfm, dt, delay, cp, 'rca','ignore','gpu');
toc

if (exist('DO_PLOTTING'))
  figure(5);
  clf;

  O_gpu = reshape(Im_tfm_gpu, Nx*Ny, Nz)';
  c_scan_gpu = reshape(max(abs(O_gpu)), Nx, Ny);
  mx = max(max(c_scan_gpu));

  imagesc(x, y, 20.0*log10(c_scan_gpu/mx));
  h_cb = colorbar;
  h_cb_title = get(h_cb,'Title');
  set(h_cb_title,'String','Normalized Amplitude [dB]')
  axis square;
  xlabel('x [mm]');
  ylabel('y [mm]');
  title('C-scan RCA GPU DAS beamformed data');
end
