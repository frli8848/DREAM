%%
%% ------------- Rectangular Transducer --------------------------
%%

setup_dream_parameters

% Geometrical parameters.
a = 10;				% x-size.
b = 15;				% y-size.
geom_par = [a b];

[H,err] = dreamrect(Ro,geom_par,s_par,delay,m_par,'stop');

alpha  = 5.0;                   % Absorbtion (dB/cm Hz).
[Hatt,err] = dreamrect(Ro,geom_par,s_par,delay,[v cp alpha],'stop');

%% Sampled analytical
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
  xlabel('t [\mus]')
  ylabel('h_{SIR} [m/s]')
  grid('on');
  legend('Attenuation {\alpha} = 0 [dB/cm MHz]',...
       'Attenuation {\alpha} = 5 [dB/cm MHz]',...
       'Analytical');
end

title('Rectangular transducer')
