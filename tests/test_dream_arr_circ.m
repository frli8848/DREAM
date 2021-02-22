%%
%% ------------- Array with rectangular elements --------------------------
%%

setup_dream_parameters

% Element size [mm].
R = 0.5;
geom_par = [R];

% Grid function (position vectors of the elements).
x = -10:1:10;
[gx,gy] = meshgrid(x);
gx = gx(:);
gy = gy(:);
gz = zeros(length(gx),1);
G = [gx gy gz];

%% Focusing parameters.
foc_met = 'off';
%foc_met = 'x';
%foc_met = 'y';
%foc_met = 'xy';
%foc_met = 'x+y';
focal = 100;                    % Focus radius

%% Beam steering.
steer_met = 'off';
%steer_met = 'x';
%steer_met = 'y';
%steer_met = 'xy';

theta  = 0;                     % Angle in x-direction.
phi    = 0;                     % Angle in y-direction.
steer_par = [theta phi];

%% Apodization.
apod_met = 'off';
%apod_met = 'ud';                       % User defined.
%apod_met = 'triangle';
%apod_met = 'gauss';
%apod_met = 'raised';                   % Raised cosine.
%apod_met = 'simply';                   % Simply supported.
%apod_met = 'clamped';
apod = ones(length(gx),1);              % Apodization weights for 'ud'.
win_par = 1;                            % Parameter for raised cos and Gaussian apodization functions.

[H,err] = dream_arr_circ(Ro,geom_par,G,s_par,delay,...
                         m_par,foc_met,focal,steer_met,steer_par,apod_met,apod,win_par,'stop');

alpha  = 5.0;                   % Absorbtion (dB/cm Hz).
[Hatt,err] = dream_arr_circ(Ro,geom_par,G,s_par,delay,...
                            [v cp alpha],foc_met,focal,steer_met,steer_par,apod_met,apod,win_par,'stop');

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
  ax = axis;
  axis([0 50 ax(3) ax(4)]);
  xlabel('t [\mus]')
  ylabel('h_{SIR} [m/s]')
  grid('on');
  legend('Attenuation {\alpha} = 0 [dB/cm MHz]',...
         'Attenuation {\alpha} = 5 [dB/cm MHz]');

end
title('Array with cirular elements')
