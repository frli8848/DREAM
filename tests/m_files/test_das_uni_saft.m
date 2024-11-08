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

% Geometrical parameters.
a = 0.4;                        % x-size.
b = 15;				% y-size.
geom_par = [a b];

%%
%% SAFT data
%%

min_tr = -25.0;
pitch_tr = 0.5;
max_tr = 25.0;

xo = (min_tr:pitch_tr:max_tr);
yo = zeros(length(xo),1);
zo = z_pt*ones(length(xo),1);
Ro = [xo(:) yo(:) zo(:)];

delay = 0.0;
[H,err] = dreamrect(Ro,geom_par,s_par,delay,m_par,'stop');

Hdp = fftconv_p(H,H); % Double-path SIRs
Ysaft = fftconv_p(Hdp,h_e);
Ysaft = Ysaft/max(max(abs(Ysaft))); % Normalize amplitudes

if (exist('DO_PLOTTING'))
  figure(2);
  clf;
  t_dp = 0:Ts:Ts*(size(Ysaft,1)-1);
  imagesc(xo,t_dp,Ysaft)
  title('SAFT B-scan')
  xlabel('x [mm]')
  ylabel('t [{\mu}s]')
end

gt = [min_tr, pitch_tr max_tr];
gr = [];

%% Observation points for DAS
min_Rx = -25.0;
dx = 0.5;
max_Rx = 25.0;
x = min_Rx:dx:max_Rx;
z = (0:63)/64*20; % Use the same as for the DAS TFM test.

min_Rz = min(z);
dz = z(2)-z(1);
max_Rz = max(z);

ro_saft = [min_Rx, dx,  max_Rx;
          0.0,    dx, 0.0;
          min_Rz, dz,  max_Rz;];

delay = system_delay; % Compensate for the pulse/system (transducer) delay.
Im_saft_gpu = das_uni(Ysaft, gt, gr, ro_saft, dt, delay, cp, 'saft'); % NB. There is no CPU version

if (exist('DO_PLOTTING'))
  figure(3);
  clf;
  imagesc(x,z,reshape(Im_saft_gpu,length(z),length(x)))
  title('SAFT GPU Uniform Grid Reconstruction')
  xlabel('x [mm]')
  ylabel('z [mm]')
end

%%
%% Single precision
%%

disp('Single precision');

Ysaft_f = single(Ysaft);
gt_f = single(gt);
gr_f = single(gr);
ro_saft_f = single(ro_saft);

delay = system_delay; % Compensate for the pulse/system (transducer) delay.
Im_saft_f_gpu = das_uni(Ysaft_f, gt_f, gr_f, ro_saft_f, single(dt), single(delay), single(cp), 'saft');

if (exist('DO_PLOTTING'))
  figure(4);
  clf;
  imagesc(x,z,reshape(Im_saft_f_gpu,length(z),length(x)))
  title('SAFT GPU Uniform Grid Single Reconstruction')
  xlabel('x [mm]')
  ylabel('z [mm]')
end
