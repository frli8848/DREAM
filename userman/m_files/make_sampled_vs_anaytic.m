
SAVE_EPS = 1;
eps_path = 'eps/';

Fs = 25;   % Sampling freq. in MHz.
Ts = 1/Fs;


z = 50.0; % [mm]

%xo = 25;
xo = 0;
yo = 0;
zo = z;
ro = [xo yo zo];
Ro = ro;
disp(['Observation point (x,y,z) = ' num2str(Ro)])

% Descretization parameters.
dt = Ts; 				% [us]. 
nt = 1500;   				% Length of spatial impulse response vector.
s_par = [dt nt];

% Material parameters.
v     = 1.0; 				% Normal velocity.
cp    = 1500; 				% Sound speed.
alfa  = 0; 				% Absorbtion [dB/(cm MHz)].
m_par = [v cp alfa];

% Delay.
t_z = z*1e3/cp;
delay = 0; 				% Start at 0 [us]. 


figure(1);
clf;

% Geometrical parameters.
r = 1.2;
geom_par = [r];
t = (0:Ts:Ts*(nt-1))+delay;

h_s_sir = scirc_sir(Ro,geom_par,s_par,delay,m_par(1:2),1000);

% Continous response.
c_dt = 1/2000; 				% [us]. 
c_nt = 150000;                          % Length of spatial impulse response vector.
c_par = [c_dt c_nt];

h_c_sir = circ_sir(Ro,geom_par,c_par,delay,m_par(1:2));
t_c = (0:c_dt:c_dt*(c_nt-1)) + delay;
start = 800;
stop = start + 300;

h1 = stem(t(start:stop),h_s_sir(start:stop),'r');
set(h1,'LineWidth',2);

hold('on');

h2 = plot(t_c,h_c_sir,'-');
set(h2,'LineWidth',3);

xlabel('t [\mus]','FontSize',20);
ylabel('SIR Amplitude [m/s]','FontSize',20);
%title('SIR for a Circular Transducer','FontSize',16);
set(gca,'FontSize',16);

ax = axis;
axis([33.25 33.5 0 ax(4)+ax(4)*0.1]) % Large 3 mm.

h_l = legend(['Sampled SIR'],['Analytic SIR']);

for n=830:860
  h_l = line([t(n) t(n)]+Ts/2, [0 ax(4)+ax(4)*0.1]);
  set(h_l,'LineStyle','-.');
end

grid('on');

eps_file = 'sampled_vs_analythic_small';

if SAVE_EPS
  print('-depsc',[eps_path 'circ_sir_' eps_file '.eps']);
end



figure(2);
clf;

% Geometrical parameters.
r = 3;
geom_par = [r];
t = (0:Ts:Ts*(nt-1))+delay;

h_s_sir = scirc_sir(Ro,geom_par,s_par,delay,m_par(1:2),1000);

% Continous response.
c_dt = 1/2000; 				% [us]. 
c_nt = 150000;                          % Length of spatial impulse response vector.
c_par = [c_dt c_nt];

h_c_sir = circ_sir(Ro,geom_par,c_par,delay,m_par(1:2));
t_c = (0:c_dt:c_dt*(c_nt-1)) + delay;
start = 800;
stop = start + 300;

%h1 = plot(t(start:stop),h_s_sir(start:stop));
h1 = stem(t(start:stop),h_s_sir(start:stop),'r');
set(h1,'LineWidth',2);

hold('on');

h2 = plot(t_c,h_c_sir,'-');
set(h2,'LineWidth',3);

%legend(['a = ' num2str(r1) ' [mm]'],['a = ' num2str(r2) ' [mm]'])

xlabel('t [\mus]','FontSize',20);
ylabel('SIR Amplitude [m/s]','FontSize',20);
%title('SIR for a Circular Transducer','FontSize',16);
set(gca,'FontSize',16);

ax = axis;
axis([33.25 33.5 0 ax(4)+ax(4)*0.1]) % Large 3 mm.

h_l = legend(['Sampled SIR'],['Analytic SIR']);

for n=830:860
  h_l = line([t(n) t(n)]+Ts/2, [0 ax(4)+ax(4)*0.1]);
  set(h_l,'LineStyle','-.');
end

grid('on');

eps_file = 'sampled_vs_analythic_large';

if SAVE_EPS
  print('-depsc',[eps_path 'circ_sir_' eps_file '.eps']);
end
