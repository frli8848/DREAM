function [h,t,f_e,f,p_delay] = el_imp_resp(Ts)
% [h,f_e,t,p_delay] = el_imp_resp(Ts)
%
% Get el. imp. resp. from hydrophone data.    
% 
%
% Fredrik Lingvall 2005-10-18.

%dp = '/disk2/users/staff/fl/projects/phd/array/data';
%[D,xd,yd,us] = ndt2mat([dp '/el' num2str(32) '.ndt']);
%[Mh,Nh,Kh] = size(D);
%B = reshape(D,Nh,Kh);
%B = B - mean(mean(B));
%[Bh,xh] = hist(B(:),10000);
%[mx,ind] = max(hist(B(:),10000));
%bias = xh(ind);
%B = B - bias;
%h = B(177,64:130);
%h = h(:);
%h = conv(h,h); % Forward + backward el. imp. reponse (assumed to be identical). 


%
% Simulated response.
%
nt=50;
t = (0:((nt-1)))*Ts;
%t = t(1:200);

f0 = 3;                             % Center frequency [MHz].
fprintf('\nf = %1.2f [MHz]\n',f0);
%%lambda = cp/f0/1e3; % [mm].
%%fprintf('lambda = %1.2f [mm]\n',lambda);

t0 = 0.55;                 % Time delay to max amplitude [mus].
a_n = 20;                        % Envelop parameter.
h = -exp(-a_n.*(t-t0).^2).*cos(2.*pi.*f0.*t);


% Number of samples to max amplitude of pulse.
%[mx,p_delay] = max(abs(hilbert(h)));
[mx,p_delay] = max(h);

%h = diff(h);
%h = diff(diff(h));

% Don't need to deconv the SIR since its a dirac at focus!
%h_sir = dreamcylind_f([0 0 190],[0.9 33 190],[0.05 0.05 Ts ...
%	  100], 190/cp*1000, [1.0 cp 0]); 

%h_e = h(7:100);
h_e = h;
f_tmp = abs(freqz(h_e,1,1024));
h_e = h_e/max(f_tmp); % Unity gain at center freq.

f_e = abs(freqz(h_e,1,1024));

f = (0:1023)/1045/Ts/2;
%f_e = f_e(1:512);
