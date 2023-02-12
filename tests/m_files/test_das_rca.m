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

if (exist('DO_PLOTTING'))
  figure(2);
  clf;
  t_dp = 0:Ts:Ts*(size(Yfmc,1)-1);
  imagesc(1:L^2,t_dp,Yfmc)
  title('FMC RCA B-scan')
  xlabel('A-scan index')
  ylabel('t [{\mu}s]')

  figure(3);
  clf;
  t_dp = 0:Ts:Ts*(size(Yfmc,1)-1);
  imagesc(xo,t_dp,Yfmc(:,1:L:end))
  title('FMC RCA B-scan when transmit element = receive element')
  xlabel('x [mm]')
  ylabel('t [{\mu}s]')

end

num_elements = size(xo,2);
Gt = [xo(:) zeros(num_elements,1) zeros(num_elements,1)];
Gr = [zeros(num_elements,1) xo(:) zeros(num_elements,1)];

%% Observation points for DAS
x = -25:0.5:25;
z = (0:63)/64*20; % Make sure its a factor of 64 (the OpenCL work group size).
[X,Z] = meshgrid(x,z);
Y = zeros(size(X));
Ro_rca = [X(:) Y(:) Z(:)];

delay = system_delay; % Compensate for the pulse/system (transducer) delay.
Im_rca = das(Yfmc, Gt, Gr, Ro_rca, dt, delay, cp,'rca');

if (exist('DO_PLOTTING'))
  figure(4);
  clf;
  imagesc(x,z,reshape(Im_rca,length(z),length(x)))
  title('RCA Reconstruction')
  xlabel('x [mm]')
  ylabel('z [mm]')
end

Im_rca_gpu = das(Yfmc, Gt, Gr, Ro_rca, dt, delay, cp, 'rca', 'ignore','gpu');

if (exist('DO_PLOTTING'))
  figure(5);
  clf;
  imagesc(x,z,reshape(Im_rca_gpu,length(z),length(x)))
  title('RCA Reconstruction')
  xlabel('x [mm]')
  ylabel('z [mm]')
end

%%
%% Single precision
%%

disp('Single precision');

Yfmc_f = single(Yfmc);
Gt_f = single(Gt);
Gr_f = single(Gr);
Ro_rca_f = single(Ro_rca);

Im_rca_f = das(Yfmc_f, Gt_f, Gr_f, Ro_rca_f, single(dt), single(delay), single(cp), 'rca');

if (exist('DO_PLOTTING'))
  figure(5);
  clf;
  imagesc(x,z,reshape(Im_rca_f,length(z),length(x)))
  title('RCA CPU Single Reconstruction')
  xlabel('x [mm]')
  ylabel('z [mm]')
end

Im_rca_f_gpu = das(Yfmc_f, Gt_f, Gr_f, Ro_rca_f, single(dt), single(delay), single(cp), 'rca', 'ignore', 'gpu');

if (exist('DO_PLOTTING'))
  figure(6);
  clf;
  imagesc(x,z,reshape(Im_rca_f_gpu,length(z),length(x)))
  title('RCA GPU Single Reconstruction')
  xlabel('x [mm]')
  ylabel('z [mm]')
end
