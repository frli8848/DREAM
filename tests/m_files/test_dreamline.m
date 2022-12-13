%%
%% ------------- Line Strip Transducer --------------------------
%%

setup_dream_parameters

% Geometrical parameters.
a = 10;				% x-size [mm].
                                %width = dy;				% Strip size (= dy).
geom_par = [a];

[H,err] = dreamline(Ro,geom_par,s_par,delay,m_par,'stop');

alpha  = 5.0;                   % Absorbtion (dB/cm Hz).
[Hatt,err] = dreamline(Ro,geom_par,s_par,delay,[v cp alpha],'stop');

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
           'Attenuation {\alpha} = 5 [dB/cm MHz]')
  end

  title('Line strip transducer')
  fprintf('dreamline\n');

end
