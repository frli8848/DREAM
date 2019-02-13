SAVE_EPS = 1;

eps_path = 'eps/';


Fs = 25;   % Sampling freq. in MHz.
Ts = 1/Fs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Observation point(s).
%
z = 20; % [mm]

% One point.
xo =0;
yo = 0;
zo = z;
ro = [xo yo zo];
Ro = ro;
%disp(['Observation point (x,y,z) = ' num2str(Ro)])

% Descretization parameters.
dx = 0.01;                              % [mm].
dy = 0.01;                              % [mm]
dt = Ts;                                % [us].
nt = 1500;                              % Length of spatial impulse response vector.
s_par = [dx dy dt nt];

% Material parameters.
v     = 1.0;                            % Normal velocity.
cp    = 1500;                           % Sound speed.
alfa  = 0;                              % Absorbtion [dB/(cm MHz)].
m_par = [v cp alfa];

% Delay.
t_z = z*1e3/cp;
delay = 0;                              % Start at 0 [us].
%delay = t_z;                           % Start at t_z [us].



% Geometrical parameters.
%r = 0.9
r = 10;
geom_par = [r];


c_dt = 1/2000;                          % [us].
c_nt = 150000;                                  % Length of spatial impulse response vector.
c_par = [c_dt c_nt];

h_c_sir = circ_sir(Ro,geom_par,c_par,delay,m_par(1:2));

t_c = (0:c_dt:c_dt*(c_nt-1)) + delay;


figure(1);
clf

h = plot(t_c,h_c_sir);
set(h,'LineWidth',2);
ax = axis;
axis([0 30 0 1800]);
grid('on');
set(gca,'FontSize',16);

xlabel('t [\mus]','FontSize',20);
ylabel('SIR Amplitude [m/s]','FontSize',20);
%axis('manual');
%axis;

%axis('tight');

if SAVE_EPS;
  print('-depsc',[eps_path 'on_axis_close.eps']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%
% Observation point(s).
%
z = 80; % [mm]

% One point.
xo =0;
yo = 0;
zo = z;
ro = [xo yo zo];
Ro = ro;
%disp(['Observation point (x,y,z) = ' num2str(Ro)])

% Descretization parameters.
dx = 0.01;                              % [mm].
dy = 0.01;                              % [mm]
dt = Ts;                                % [us].
nt = 1500;                              % Length of spatial impulse response vector.
s_par = [dx dy dt nt];

% Material parameters.
v     = 1.0;                            % Normal velocity.
cp    = 1500;                           % Sound speed.
alfa  = 0;                              % Absorbtion [dB/(cm MHz)].
m_par = [v cp alfa];

% Delay.
t_z = z*1e3/cp;
delay = 0;                              % Start at 0 [us].
%delay = t_z;                           % Start at t_z [us].



% Geometrical parameters.
%r = 0.9
r = 10;
geom_par = [r];


c_dt = 1/2000;                          % [us].
c_nt = 150000;                                  % Length of spatial impulse response vector.
c_par = [c_dt c_nt];

h_c_sir = circ_sir(Ro,geom_par,c_par,delay,m_par(1:2));

t_c = (0:c_dt:c_dt*(c_nt-1)) + delay;


figure(2);
clf

h = plot(t_c,h_c_sir);
set(h,'LineWidth',2);
ax = axis;
axis([40 70 0 1800]);
grid('on');
set(gca,'FontSize',16);

xlabel('t [\mus]','FontSize',20);
ylabel('SIR Amplitude [m/s]','FontSize',20);
%axis('manual');
%axis;

%axis('tight');

if SAVE_EPS;
  print('-depsc',[eps_path 'on_axis_far.eps']);
end
