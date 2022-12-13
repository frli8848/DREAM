%%
%% ------------- Focused Circular Transducer --------------------------
%%

setup_dream_parameters

%% Focused/concave

% Geometrical parameters.
a = 10;				% x-size of the transducer.
b = 20;				% y-size of the transducer.
Rcurv = 10;                    % Radius of the curvature.
geom_par = [a b Rcurv];

[Hf,err] = dreamcylind(Ro,geom_par,s_par,delay,m_par,'stop');

alpha  = 5.0;                   % Absorbtion (dB/cm Hz).
[Hf_att,err] = dreamcylind(Ro,geom_par,s_par,delay,[v cp alpha],'stop');

if (exist('DO_PLOTTING'))

  figure(1);
  clf;

  if size(Hf,2)>1
    mesh(xo,t,Hf);
    axis('tight');
    view(135,32);
  else
    subplot(211)
    plot(t,Hf,'b');
    hold on;
    plot(t,Hf_att,'r');
    ax = axis;
    axis([0 50 ax(3) ax(4)]);
    xlabel('t [\mus]')
    ylabel('h_{SIR} [m/s]')
    grid('on');
    legend('Attenuation {\alpha} = 0 [dB/cm MHz]',...
           'Attenuation {\alpha} = 5 [dB/cm MHz]');
  end
  title_str = sprintf('Focused cylindrical transducer at Rcurv=%1.1f [mm]',Rcurv);
  title(title_str)

end

%% Defocused/convex

% Geometrical parameters.
a = 10;				% x-size of the transducer.
b = 20;				% y-size of the transducer.
Rcurv = -10;                   % Radius of the curvature (negative)..
geom_par = [a b Rcurv];

[Hd,err] = dreamcylind(Ro,geom_par,s_par,delay,m_par,'stop');

alpha  = 5.0;                   % Absorbtion (dB/cm Hz).
[Hd_att,err] = dreamcylind(Ro,geom_par,s_par,delay,[v cp alpha],'stop');


if (exist('DO_PLOTTING'))

  figure(2);
  clf;

  if size(Hd,2)>1
    figure(2)
    clf
    mesh(xo,t,Hd);
    axis('tight');
    view(135,32);
  else
    subplot(212)
    plot(t,Hd,'b');
    hold on;
    plot(t,Hd_att,'r');
    ax = axis;
    axis([0 50 ax(3) ax(4)]);
    xlabel('t [\mus]')
    ylabel('h_{SIR} [m/s]')
    grid('on');
    legend('Attenuation {\alpha} = 0 [dB/cm MHz]',...
           'Attenuation {\alpha} = 5 [dB/cm MHz]');
  end
  title_str = sprintf('Defocused cylindrical transducer at Rcurv=%1.1f [mm]',Rcurv);
  title(title_str)

end
