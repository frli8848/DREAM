%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make sound field snapshots and movies using the
% DREAM toolbox.
%
%
% Fredrik Lingvall : 2005-10-25
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if findstr(computer,'linux') && exist('OCTAVE_VERSION')
  disp('***************************************************************');
  disp('*');
  disp('* Your are running Octave on a Linux platform!');
  disp('* You need ffmpeg installed to run this script!');
  disp('*');
  disp('***************************************************************');
  pause(4);
end

% Paths to eps, png, and jpegs.
eps_path = '../eps/';
png_path = 'png/';


eps_file = '50_m20';

% Switches
COMPUTE_SIR = 1;
GENERATE_EPS = 0;
MAKE_MOVIE = 1;
AUTO_SCALE = 0;

mkdir('png');

% 25 MHz
Ts = 0.04;

if COMPUTE_SIR

  % Soundspeed in water (in [m/s]).
  cp = 1500;

  % Length of A-scans.
  K = 1000;

  %
  % Define the array.
  %

  % Spacing between array elements (array pitch).
  arr_samp = 1;                         % [mm]

  % Element size.
  a = 0.9;                              % [mm]
  b = 1.0;                              % [mm]

  % Transmit.
  lt = 8;
  xt =  arr_samp*(-lt+0.5:lt-0.5);
  trans = [xt; a*ones(1,2*lt);];

  fprintf('Array aperture = %1.1f [mm]\n',xt(length(xt))-xt(1));


  % Number of transducer elements.
  L = size(xt,2);

  % Time vector.
  us = (0:Ts:Ts*40000);

  % Max width of ROI.
  x_max = 25;

  idx_z_min = 1;
  idx_z_max = 918;

  % Spatial sampling resolution.
  spat_samp_x = 0.3; % Horizontal resolution [mm].
  spat_samp_z = 0.3; % Vertical resoltion [mm].

  %
  % Horizontal grid.
  %
  x = -x_max:spat_samp_x:x_max;
  N = length(x);

  %
  % Vertical grid.
  %
  z_min = us(idx_z_min)*cp/1e3;
  z_max = us(idx_z_max)*cp/1e3;
  z = (z_min:spat_samp_z:z_max);
  M = length(z);

  %
  % Observation grid matrix.
  %
  [X,Z] = meshgrid(x,z);
  Y = zeros(size(X));
  xo = X(:);
  yo = Y(:);
  zo = Z(:);
  Ro = [xo yo zo];

  % Number of array elements.
  L = size(trans,2);

  % Number of observation points
  MN = size(Ro,1);

  fprintf('z_min = %1.1f [mm]\n',z_min);
  fprintf('z_max = %1.1f [mm]\n',z_max);

  fprintf('KL x MN = %d x %d \n',K*L,M*N);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Electrical impulse response.
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  he_nt=50;
  t = (0:((he_nt-1)))*Ts;

  f0 = 3;                             % Center frequency [MHz].
  fprintf('\nf = %1.2f [MHz]\n',f0);
  lambda = cp/f0/1e3; % [mm].
  fprintf('lambda = %1.2f [mm]\n',lambda);

  t0 = 0.55;                 % Time delay to max amplitude [mus].
  a_n = 20;                        % Envelop parameter.
  h_e = -exp(-a_n.*(t-t0).^2).*cos(2.*pi.*f0.*t);

  % Number of samples to max amplitude of pulse.
  [mx,p_delay] = max(h_e);

  f_tmp = abs(freqz(h_e,1,1024));
  h_e = h_e/max(f_tmp); % Unity gain at center freq.
  f_e = abs(freqz(h_e,1,1024));
  f = (0:1023)/1045/Ts/2;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % DREAM parameters.
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Descretization parameters.
  dx = 0.05;                            % [mm]
  dy = 0.05;                            % [mm]
  dt = Ts;                              % [us]
  nt = 1300;                            % Length of single-path spatial impulse response vector.
  s_par = [dx dy dt nt];


  fprintf('\nDREAM length nt = %d\n\n',nt);

  % Material parameters.
  v     = 1.0;                          % Normal velocity.
  alfa  = 0;                            % Absorbtion (dB/cm Hz).
  m_par = [v cp alfa];

  % Element size for transmit [mm].
  geom_par = [a b];

  % Grid function (position vectors of the transmit elements).
  gx = trans(1,:);
  gx = gx(:);
  gy = zeros(length(gx),1);
  gz = zeros(length(gx),1);
  G = [gx gy gz];

  %
  % User defined focusing.
  %

  foc_dep = 50;
  theta = -20;
  foc_met = 'ud';

  xp = foc_dep*tan(theta*pi/180)*ones(L,1);
  zp = foc_dep*ones(L,1);

  r_max = 0;
  for l=1:L
    tmp_max = norm([xp(l) zp(l)]' - [gx(l) 0]');
    if tmp_max > r_max
      r_max = tmp_max;
    end
  end

  focal = zeros(L,1);
  for n=1:L
    focal(n) = (r_max - norm([xp(n) zp(n)]' - [gx(n) 0]'))*1000/cp;
  end

  fprintf('Array focused at  = %1.1f [mm]\n',foc_dep);
  fprintf('Array steered = %1.1f [deg]\n',theta);

  % Beam steering
  steer_met = 'off'; % Done above by foc_met = 'ud' above.

  phi = 0;                              % Angle in y-direction.
  steer_par = [theta phi];              % Not used here.

  apod_met = 'off';

  % Apodization
  apod_met = 'off';
  apod = ones(length(gx),1);            % Not used here.
  win_par = 1;                          % Not used here.


  z_min = min(Ro(:,3));                 % Depth to closest observation point.
  delay = z_min/cp*1000;                % Single path delay.

  mx_td = max(focal);
  fprintf('----------------------------------------------------------------------\n');
  fprintf(' Max steer delay = %1.2f [us] (= %1.2f [mm], = %d [samples])\n',mx_td, mx_td*1000/cp,...
          round(mx_td/Ts));
  fprintf('----------------------------------------------------------------------\n');

  %
  % Compute transmit SIRs.
  %

  eval('tic');

  %% Serial.
  H_t = dream_arr_rect(Ro,geom_par,G,s_par,delay,m_par,foc_met,...
                       focal,steer_met,steer_par,apod_met,apod,win_par,'stop');

  t = toc;
  if t > 60
    fprintf('elapsed time: dream_arr_rect = %1.1f [min]\n',t/60);
  else
    fprintf('elapsed time: dream_arr_rect = %1.1f [s]\n',t);
  end

  %%
  %% Compute pressure response (i.e., convolve with electro-acoustical impulse response).
  %%

  eval('tic');

  %%P_t = conv_p(H_t,h_e);
  P_t = fftconv_p(H_t,h_e);

  t = toc;
  if t > 60
    fprintf('elapsed time: fftconv_p = %1.1f [min]\n',t/60);
  else
    fprintf('elapsed time: fftconv_p = %1.1f [s]\n',t);
  end
  fprintf('----------------------------------------------------------------------\n');
end

l_z = length(z);

if exist('OCTAVE_VERSION')
  figure(1);
  clf
else
  fig=figure(1);
  set(fig,'DoubleBuffer','on');
  set(gca,'xlim',[min(x) max(x)],'ylim',[z(1) z(l_z)],...
          'NextPlot','replace','Visible','off')

end

kk = 0;
for k=1:3:1100
  kk = kk+1;

  fprintf('Snapshot taken at %1.1f [us]\n', delay+k*Ts)

  ds = Ts*cp/1000;
  u_samp = ceil(spat_samp_z/ds); % n times undersampled.

  if u_samp > 1
    % Take mean of u_samp samples.
    if (round(k-u_samp/2) >=1)
      o = sum(abs(hilbert(P_t(round(k-u_samp/2:k+u_samp/2),:))))'/u_samp;
    else
      o = sum(abs(hilbert(P_t(round(k:k+u_samp/2),:))))'/u_samp;
    end
    O = reshape(o,M,N);
  else
    o = P_t(k,:);
    O = reshape(o,M,N);
  end

  mx = max(abs(o));

  if AUTO_SCALE
    O = 64*abs(O)/mx;
    Im = abs(O);
  else
    g = 0.7;
    Im = g*abs(O);
  end

  if exist('OCTAVE_VERSION')
    image(x,z,Im);
    axis('ij')
    xlabel('x [mm]');
    ylabel('z [mm]');
  else
    h = image(x,z,Im);
    set(gca,'FontSize',16);
    xlabel('x [mm]','FontSize',20);
    ylabel('z [mm]','FontSize',20);
    set(gca,'YTick', 0:10:80)
    drawnow
  end
  grid('off');

  t_str = sprintf('Wavefield Snapshot at %1.1f [us]', delay+k*Ts);
  title(t_str);

  if MAKE_MOVIE

    if exist('OCTAVE_VERSION')
      print(sprintf ('%s/snapshot_%.5d.png', png_path, k));

      if k<10
        eval(['print("' png_path 'snapshot_000' num2str(k) '.png", "-dpng", "-S640,480");']);
      elseif k<100
        eval(['print("' png_path 'snapshot_00' num2str(k) '.png", "-dpng", "-S640,480");']);
      elseif k<1000
        eval(['print("' png_path 'snapshot_0' num2str(k) '.png", "-dpng", "-S640,480");']);
      else
        eval(['print("' png_path 'snapshot_' num2str(k) '.png", "-dpng", "-S640,480");']);
      end
    else % Matlab
      % Add a frame to Matlab movie.
      drawnow;
      snapshots(kk) = getframe;
    end
  end

end % for

if MAKE_MOVIE

  if findstr(computer,'linux') && exist('OCTAVE_VERSION')
    system(['ffmpeg -r 1/5 -start_number 0 -i png/snapshot_%04d.png '...
              '-c:v libx264 -r 30 -pix_fmt yuv420p snapshots.mp4']);
  else
    % Create an AVI movie from MATLAB movie.
    movie2avi(snapshots,'snapshots.avi');
  end

end

t = toc;
if t > 60
  fprintf('elapsed time = %1.1f [min]\n',t/60);
else
  fprintf('elapsed time = %1.1f [s]\n',t);
end

if GENERATE_EPS
  t_str = sprintf('%1.1f', delay+k*Ts);
  if exist('OCTAVE_VERSION')
    eval(['print("' eps_path 'snapshot_' eps_file '_'  t_str '.eps", "-depsc");']);
  else
    eval(['print -deps ' eps_path 'snapshot_' eps_file '_' t_str '.eps']);
  end
end

%
% Electrical impulse response.
%

figure(2);
clf
plot(us(1:length(h_e)),h_e/max(h_e));
grid('on')
axis([0 max(us(length(h_e))) -1.2 1.2])

title('Pulse');

if exist('OCTAVE_VERSION')
  xlabel('t [us]');
  ylabel('Normalized Amplitude');
else
  xlabel('t [\mus]','FontSize',20);
  ylabel('Normalized Amplitude','FontSize',20);
  set(gca,'FontSize',16);
end

if GENERATE_EPS
  if exist('OCTAVE_VERSION')
    eval(['print("' eps_path 'imp_resp.eps", "-depsc");']);
  else
    eval(['print -depsc ' eps_path 'imp_resp.eps']);
  end
end

figure(3);
clf
plot(f,f_e);
grid('on')
axis([0 max(f) 0 1.2])

if exist('OCTAVE_VERSION')
  xlabel('f [MHz]');
  ylabel('Normalized Amplitude');
else
  xlabel('f [MHz]','FontSize',20);
  ylabel('Normalized Amplitude','FontSize',20);
  set(gca,'FontSize',16);
end

title('Pulse Spectrum');

if GENERATE_EPS
  if exist('OCTAVE_VERSION')
    eval(['print("' eps_path 'freq_imp_resp.eps", "-depsc");']);
  else
    eval(['print -depsc ' eps_path 'freq_imp_resp.eps']);
  end
end
