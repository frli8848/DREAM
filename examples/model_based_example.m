
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A model based array imaging example.
%
%
% This example illustrates model based imaging with a
% 32-element linear (concave) array. The array is focused
% at 50 mm and reconstructed images for 7 point targets,
% also at 50 mm depth, is shown. Note: you probably need 
% >= 2 GB of RAM to run this script using the default 
% parameter values.
%
% Outline of the script:
%
% 1) Build a model for a 3.5 MHz, 32 element linear array 
%    (with concave elements).
%
% 2) Simulate the system (show B-scan and el. impulse resp.) 
%
% 3) Optimal linear beamforming.

% 4) Delay-and-sum imaging (optional).
%
% Copyright (C) 2008,2009 Fredrik Lingvall.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% $Revision: 565 $ $Date: 2009-09-17 22:24:06 +0200 (Thu, 17 Sep 2009) $ $LastChangedBy: dream $

%
% Flags
% 

%DAS = 0; % Set to 1 to enable delay-and-sum processing
DAS = 1; 


% Number of cpus/cores on your system.
n_cpus = 2;

%
% Temporal sampling parameters.
%

Ts = 1/25; % Sampling freq. in [us].
us = (0:Ts:Ts*40000);


% Soundspeed [m/s].
cp = 1500;

%
% Simulate a 3.5 MHz array.
%
he_nt = 50;                             % Length of the electrical impulse response.
t = (0:((he_nt-1)))*Ts;

f0 = 3.5;                               % Center frequency [MHz].
fprintf('\nf = %1.2f [MHz]\n',f0);
lambda = cp/f0/1e3;                     % [mm].
fprintf('lambda = %1.2f [mm]\n',lambda);
t0 = 0.55;                              % Time delay to max amplitude [mus].
a_n = 20;                               % Envelop parameter.
A = 1e-5;                               % Amplitude.
h_e = -A*exp(-a_n.*(t-t0).^2).*cos(2.*pi.*f0.*t);
%h_e = h_e - mean(h_e);


% Set length of A-scans. The max length of 
% K is 2*nt-1 + lenght(h_e) - 1.
nt = 410;
K = 2*nt + length(h_e) - 2;

%
% Define the observation points.
%

x_max = 25;

z_min = 47.5;
z_max = 52.5;

spat_samp = 0.5; % 1 mm.
spat_samp_z = 0.25; % Samples / mm. !!!!!

% Horizontal grid.
x = -x_max:spat_samp:x_max;
N = length(x);

% Vertical grid size.
z = z_min:spat_samp_z:z_max;
M = length(z);

% Observation grid matrix.
[X,Z] = meshgrid(x,z);
Y = zeros(size(X));
xo = X(:);
yo = Y(:);
zo = Z(:);
Ro = [xo yo zo];

% Number of observation points
MN = size(Ro,1);

fprintf('\n');
disp('********************************************')
disp('* Propagation matrix computation started *')
disp('********************************************')
disp(' ')

%
% DREAM parameters.
%

% Descretization parameters.
dx = 0.05; 				% [mm]
dy = 0.05; 				% [mm]
dt = Ts; 				% [us] 
                                        %nt = 400;   				% Length of single-path 
                                        % spatial impulse response vector.
s_par = [dx dy dt nt];

fprintf('\nDREAM length nt = %d\n\n',nt);

% Material parameters.
v     = 1.0; 				% Normal velocity.
alfa  = 0; 				% Absorbtion (dB/cm Hz).
m_par = [v cp alfa];

% Element size [mm].
a = 0.9;
b = 10;
R = 50; % Focus depth.
geom_par = [a b R];

fprintf('Element size: a=%1.2f x b=%1.2f, R=%1.2f\n',a,b,R);

% Grid function (position vectors of the transmit elements).
gx = -15.5:1:15.5; % A 32 element array.
gx = gx(:);
gy = zeros(length(gx),1);
gz = zeros(length(gx),1);
G = [gx gy gz];

% Number of array elements.
L = length(gx);

% Focusing parameters.

foc = 50; % [mm]

%foc_met = 'ud';
%focal = foc;                          % foc is a vector for 'ud'.
%foc_met = 'off';
%focal = 0;
foc_met = 'x';
%foc_met = 'y';
%foc_met = 'xy';
%foc_met = 'x+y';
focal = foc; 				% Focus radius [mm].


% Beam steering.
steer_met = 'off';
%steer_met = 'x';
%steer_met = 'y';
%steer_met = 'xy';

theta = 0;                              % Angle in x-direction.
phi = 0; 				% Angle in y-direction.
steer_par = [theta phi];

% Apodization.
apod_met = 'off';
%apod_met = 'ud'; 			% User defined.
%apod_met = 'triangle';
%apod_met = 'gauss';
%apod_met = 'raised'; 			% Raised cosine.
%apod_met = 'simply'; 			% Simply supported.
%apod_met = 'clamped';
apod = ones(length(gx),1); 		% Apodization weights for 'ud'.
win_par = 1; 				% Parameter for raised cos and Gaussian apodization functions.


% Print problem size.
fprintf('\nz_min = %1.1f [mm], z_max = %1.1f [mm]\n',z_min,z_max);
fprintf('K = %d, L = %d, M = %d, N = %d.\n',K,L,M,N);
fprintf('P is of size: (K*L = %d) x (M*N = %d).\n\n',K*L,M*N);

z_min = min(Ro(:,3));                   % Depth to closest observation point.
delay = z_min/cp*1000.0;                % Single path delay.

fprintf('\n----------------------------------------------------------------------\n');

% Allacate space for the propagation matrix.
P = zeros(K*L,MN);

%
% Transmit SIRs (same for all receive elements)
%

fprintf(' Computing transmit SIRs using DREAM: ');

eval('tic');
if n_cpus == 1
  % Serial. 
  H_t = dream_arr_cylind_f(Ro,geom_par,G,s_par,delay,m_par,foc_met,...
                           focal,steer_met,steer_par,apod_met,apod,win_par,'stop');
else
  % Parallel.
  H_t = dream_arr_cylind_f_p(Ro,geom_par,G,s_par,delay,m_par,foc_met,...
                             focal,steer_met,steer_par,apod_met,apod,win_par,n_cpus,'stop');
end

t = toc;
if t > 60
  fprintf('elapsed time = %1.1f [min]\n',t/60);
else
  fprintf('elapsed time = %1.1f [s]\n',t);
end

fprintf('----------------------------------------------------------------------\n');

% Loop over all receive elements.
for l=1:L
  
  % Compute the SIRs for a new element.
  fprintf(' Computing SIRs for receive element %d using DREAM: ', l);
  eval('tic')
  
  Ro_r = [Ro(:,1)-gx(l) Ro(:,2) Ro(:,3)];
  
  if n_cpus == 1
    H_r = dreamcylind_f(Ro_r, geom_par, s_par, delay, m_par,'stop');
  else
    H_r = dreamcylind_f_p(Ro_r, geom_par, s_par, delay, m_par,n_cpus,'stop');
  end

  %
  % Compute pulse-echo responses.
  %

  if (l == 1) % Compute the fftw wisdoms only once. 
    [H_tmp,wisdom_1]  = fftconv_p(H_t,H_r,n_cpus); % Double-path.
    clear H_r;
    [H_tmp2,wisdom_2] = fftconv_p(H_tmp,h_e,n_cpus); % Electrical impulse response.
    clear  H_tmp;
  else
    H_tmp  = fftconv_p(H_t,H_r,n_cpus,wisdom_1); % Double-path.
    clear H_r;
    H_tmp2 = fftconv_p(H_tmp,h_e,n_cpus,wisdom_2); % Electrical impulse response.
    clear  H_tmp;
  end
  
  % Add responses for the l:th element.
  %P(K*(l-1)+1:l*K,:) = H_tmp2(1:K,:);
  copy_p(P, [K*(l-1)+1 l*K], [1 size(Ro,1)], H_tmp2, n_cpus); 
  
  clear H_tmp2;
  
  t = toc;
  
  if t > 60
    fprintf('elapsed time = %1.1f [min]\n',t/60);
  else
    fprintf('elapsed time = %1.1f [s]\n',t);
  end

end
fprintf('----------------------------------------------------------------------\n');


%
% Simulate the system.
%

fig_n = 1;

figure(fig_n);
clf;
fig_n = fig_n+1;


subplot(211);
plot(us(1:length(h_e)),h_e/A);
xlabel('Time [{\mu}s]');

title('The Double-path Eletrical Impulse Response')

f_e = abs(fft(h_e/A,1024));
f = (0:1023)/1024 * 1/Ts;

subplot(212);
plot(f(1:512),f_e(1:512));
xlabel('Frequency [MHz]');

figure(fig_n);
clf;
fig_n = fig_n+1;

% Seven point targets at z=50 mm.
o = zeros(M*N,1);
o(round(M*N/2),1) = 1; 
o(round(M*N/2-M*10),1) = 1; 
o(round(M*N/2+M*10),1) = 1; 
o(round(M*N/2-M*20),1) = 1; 
o(round(M*N/2+M*20),1) = 1; 
o(round(M*N/2-M*30),1) = 1; 
o(round(M*N/2+M*30),1) = 1; 

O = reshape(o,M,N);
mesh(x,z,O);
xlabel('x [mm]')
ylabel('z [mm]')
title('True image')

figure(fig_n);
clf;
fig_n = fig_n+1;

sigma_e = 0.1;
e = randn(K*L,1);

y = P*o + e;

% A B-scan image.
B = reshape(y,K,L);
imagesc(gx,delay+us(1:K),abs(B));
xlabel('x [mm]')
ylabel('t [{\mu}s]')
title('B-scan image')

%
% The matched filter.
%

figure(fig_n);
clf;
fig_n = fig_n+1;

o_hat_mf = P'*y;

O_hat_mf = reshape(o_hat_mf,M,N);
mesh(x,z,O_hat_mf);
xlabel('x [mm]')
ylabel('z [mm]')
title('Matched filter image')

%
% The optimal linear estimator.
%

figure(fig_n);
clf;
fig_n = fig_n+1;

sigma_o = 1;
% Fopt = inv( P'*P/sigma_e^2 + eye(MN,MN)/sigma_o^2)*P'*inv(sigma_o);
F = P'*P;
F = F + eye(MN,MN) * sigma_o^2/sigma_e^2;
%F = cholinv(F); % Works in Octave.
F = inv(F);
Fopt = F*P';
clear F;

o_hat_opt = Fopt*y;
O_hat_opt = reshape(o_hat_opt,M,N);
mesh(x,z,O_hat_opt);
xlabel('x [mm]')
ylabel('z [mm]')
title('LMMSE (optimal linear) image')

if DAS
  
  fprintf('\n');
  disp('********************************************')
  disp('* Delay-and-sum matrix computation started *')
  disp('********************************************')
  disp(' ')
  
  
  % Allacate space for the delay matrix.
  D = zeros(K*L,MN);
  
  %
  % Transmit SIRs (same for all receive elements)
  %
  
  fprintf(' Computing transmit delays: ');
  
  eval('tic');
  % Transmit "delay-and-sum" matrix. 
  D_t = das_arr(Ro,G,s_par(3:4),delay,m_par(2),foc_met,...
                focal,steer_met,steer_par,apod_met,apod,win_par,'stop');
  
  t = toc;
  if t > 60
    fprintf('elapsed time = %1.1f [min]\n',t/60);
  else
    fprintf('elapsed time = %1.1f [s]\n',t);
  end
  
  fprintf('----------------------------------------------------------------------\n');
  
  % Loop over all receive elements.
  for l=1:L
    
    % Compute the SIRs for a new element.
    fprintf(' Computing delays for receive element %d: ', l);
    eval('tic')
    
    Ro_r = [Ro(:,1)-gx(l) Ro(:,2) Ro(:,3)];
    
    D_r = das(Ro_r, s_par(3:4), delay, m_par(2),'stop');
    %D_r = dreamline(Ro_r, 100, s_par, delay, m_par,'stop');
    
    %
    % Compute pulse-echo responses.
    %
    
    if (l == 1) % Compute the fftw wisdoms only once. 
      [D_tmp,wisdom_1]  = fftconv_p(D_t,D_r,n_cpus); % Double-path delay.
      clear D_r;
      % Delay due to electrical impulse response.
      [D_tmp2,wisdom_2] = fftconv_p(D_tmp,[zeros(1,11) 1 zeros(1,50-12)]',n_cpus);

      clear  D_tmp;
    else
      D_tmp  = fftconv_p(D_t,D_r,n_cpus,wisdom_1); % Double-path delay.
      clear D_r;
      % Delay due to electrical impulse response.
      D_tmp2 = fftconv_p(D_tmp,[zeros(1,11) 1 zeros(1,50-12)]',n_cpus,wisdom_2); 
      clear  D_tmp;
    end
    
    % Add responses for the l:th element.
    %P(K*(l-1)+1:l*K,:) = D_tmp2(1:K,:);
    copy_p(D, [K*(l-1)+1 l*K], [1 size(Ro,1)], D_tmp2, n_cpus); 
    
    clear D_tmp2;
    
    t = toc;
    
    if t > 60
      fprintf('elapsed time = %1.1f [min]\n',t/60);
    else
      fprintf('elapsed time = %1.1f [s]\n',t);
    end
    
  end
  fprintf('----------------------------------------------------------------------\n');

  figure(fig_n);
  clf;
  fig_n = fig_n+1;
  
  o_hat_das = D'*y;
  
  O_hat_das = reshape(o_hat_das,M,N);
  mesh(x,z,O_hat_das);
  xlabel('x [mm]')
  ylabel('z [mm]')
  title('DAS image')
  
end % DAS


