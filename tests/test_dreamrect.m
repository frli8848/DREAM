ONE_POINT = 1;
%ONE_POINT = 0; %Many points.


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
alpha  = 0.0;                              % Absorbtion (dB/cm Hz).
m_par = [v cp alpha];

t_z = z*1e-3/cp * 1e6;          % us

delay = 0;
%delay = t_z;

% ------------- Line Strip Transducer --------------------------

% Geometrical parameters.
a = 10;				% x-size.
b = 15;				% y-size.
geom_par = [a b];

[H,err] = dreamrect(Ro,geom_par,s_par,delay,m_par,'stop');
alpha  = 5.0;                   % Absorbtion (dB/cm Hz).
[Hatt,err] = dreamrect(Ro,geom_par,s_par,delay,[v cp alpha],'stop');
[Ha] = rect_sir(Ro,geom_par,s_par(3:4),delay,m_par(1:2));




figure(1);
clf;

if size(H,2)>1
  mesh(xo,t,H);
  axis('tight');
  view(135,32);
else
  plot(t,H,'b');
  hold on;
  plot(t,Hatt,'r');
  plot(t,Ha,'g');
  ax = axis;
  axis([0 50 ax(3) ax(4)]);
  %xlabel('t [\mus]')
  grid('on');
  legend('Attenuation {\alpha} = 0 [dB/cm MHz]',...
       'Attenuation {\alpha} = 5 [dB/cm MHz]',...
       'Analytical');
end

title('Rectangular transducer')
