%%
%% ------------- Array with rectangular elements --------------------------
%%

setup_dream_parameters

% Grid function (position vector of the elements).
gr = [3.5 ...
       3.6 5 ...
       5.1 6 ...
       6.1 7 ...
       7.1 8];

G = gr(:);

% Focusing parameters.
foc_met = 'on';
%foc_met = 'off';
focal = 25;                     % Focus radius

% Apodization.
apod_met = 'off';
%apod_met = 'ud';                       % User defined.
%apod_met = 'triangle';
%apod_met = 'gauss';
%apod_met = 'raised';                   % Raised cosine.
%apod_met = 'simply';                   % Simply supported.
%apod_met = 'clamped';
apod = ones(length(gr),1);              % Apodization weights for 'ud'.
win_par = 1;                            % Parameter for raised cos and Gaussian apodization functions.

[H,err] = dream_arr_annu(Ro, G, s_par, delay,...
                         m_par,...
                         foc_met, focal,...
                         apod_met, apod, win_par,'stop');

alpha  = 5.0;                   % Absorbtion (dB/cm Hz).
[Hatt,err]= dream_arr_annu(Ro, G, s_par, delay,...
                           [v cp alpha],...
                           foc_met, focal,...
                           apod_met, apod, win_par,'stop');

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
title('Annular Array')
