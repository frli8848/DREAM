
% Generate the figures for the Lossy Media Section
% 
% Fredrik Lingvall 2007-03-30.

Fs = 25;   % Sampling freq. in MHz.
Ts = 1/Fs;

%
% Observation points.
%

%  Points along z-axis.
d  = 15;
zo = (10:d:50);
yo = zeros(length(zo),1);
xo = zeros(length(zo),1);
Ro = [xo(:) yo(:) zo(:)];


% Descretization parameters.
dt = Ts; 				% us. 
nt = 1600;   				% Length of spatial impulse response vector.
s_par = [dt nt];

t = 0:Ts:Ts*(nt-1);

% Material parameters.
cp    = 1540; 			        % Sound speed.
alfa  = 1.0; 				% Absorbtion (dB/cm MHz).
m_par = [cp alfa];

delay = 0;
%delay = Ro(:,3)/cp*1000;                % This removes the
                                        % prapagation delay in the plots.

H = dream_att(Ro,s_par,delay,m_par);


figure(1);
clf;

for n=1:size(H,2)
  switch (n)
    
   case 1
    plot(t,H(:,n),'LineWidth',2);
    
   case 2
    plot(t,H(:,n),'--','LineWidth',2);
    
   case 3
    plot(t,H(:,n),'-.','LineWidth',2);
    
   case 4
    plot(t,H(:,n),':','LineWidth',2);
    
   case 5
    plot(t,H(:,n),'r--','LineWidth',2);
    
  end
  str{n} = sprintf(['z = ' num2str(zo(n)) ' [mm]']);
  hold('on');
end
grid('on');
hold('off');
legend(str);
axis([0 40 0 0.7])
title(['Casual Greens Functions for Lossy Media (\alpha = ' num2str(alfa) ' [dB /cm MHz])']);
xlabel('Time [\mus]');

print -depsc ../eps/lossy_time.eps


figure(2);
clf;

Hf = fft(H);
f = (0:size(H,1)-1)/size(H,1)*Fs;

for n=1:size(H,2)
  switch (n)
   
   case 1
    plot(f,20*log10(abs(Hf(:,n))),'LineWidth',2);
    
   case 2
    plot(f,20*log10(abs(Hf(:,n))),'--','LineWidth',2);
    
   case 3
    plot(f,20*log10(abs(Hf(:,n))),'-.','LineWidth',2);
    
   case 4
    plot(f,20*log10(abs(Hf(:,n))),':','LineWidth',2);
    
   case 5
    plot(f,20*log10(abs(Hf(:,n))),'r--','LineWidth',2);
    
   otherwise
    plot(f,20*log10(abs(Hf(:,n))));
    
  end
  str{n} = sprintf(['z = ' num2str(zo(n)) ' [mm]']);
  hold('on');
end
hold('off');
grid('on');
legend(str);
ax = axis;
%axis([0 Fs/2 ax(3) ax(4)])
axis([0 3 -10 0])
xlabel('Frequency [MHz]');
ylabel('Amplitude [dB]');
title(['Attenuation coefficient \alpha = ' num2str(alfa) ' [dB /cm MHz]']);

print -depsc ../eps/lossy_freq.eps
