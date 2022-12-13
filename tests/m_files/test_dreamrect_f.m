%%
%% ------------- Focused Rectangular Transducer --------------------------
%%

setup_dream_parameters

% Geometrical parameters.
a = 10;				% x-size.
b = 15;				% y-size.
geom_par = [a b];

% Focusing parameters.
%foc_met = 'off';
%foc_met = 'x';
%foc_met = 'y';
foc_met = 'xy';
%foc_met = 'x+y';
focal = 10;                     % Focus radius

[H,err] = dreamrect_f(Ro,geom_par,s_par,delay,m_par,foc_met,focal,'stop');

alpha  = 5.0;                   % Absorbtion (dB/cm Hz).
[Hatt,err] = dreamrect_f(Ro,geom_par,s_par,delay,[v cp alpha],foc_met,focal,'stop');

if (exist('DO_PLOTTING'))

  figure(1);
  clf;

  if (size(H,2)>1)
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
  title_str = sprintf('Focused Rectangular transducer at z=%1.1f [mm]',focal);
  title(title_str)

end
