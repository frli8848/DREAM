%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1D convolution algorithms
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L1 = 1000;
L2 = 400;
N = 50;
X = randn(L1,N);
Y = randn(L2,N);

%%
%% Serial Time-domain
%%

Z0 = zeros(L1+L2-1,N);
eval('tic')
for n=1:N
  Z0(:,n) = conv( X(:,n), Y(:,n));
end
t0= toc;
fprintf('conv: %f [s]\n\n', t0);

if (exist('fftconv'))
  Z0_fft = zeros(L1+L2-1,N);
  eval('tic')
  for n=1:N
    Z0_fft(:,n) = fftconv( X(:,n), Y(:,n));
  end
  t0_fft= toc;
  fprintf('fftconv: %f [s]\n\n',t0_fft);
end


%%
%% Threaded time-domain
%%

eval('tic')
Z1 = conv_p(X, Y);
t1= toc;
fprintf('conv_p: %f [s]\n', t1);
fprintf('Numerical error: ||conv-conv_p|| = %e\n\n', sum(sum((Z0-Z1).^2)));

%%
%% Threaded frequency-domain
%%

eval('tic')
Z2 = fftconv_p(X, Y);
t2= toc;
fprintf('fftconv_p: %f [s]\n', t2);
fprintf('Numerical error: ||conv-fftconv_p|| = %e\n\n', sum(sum((Z0-Z2).^2)));

%%
%% Threaded frequency-domain using the overlap-and-add method
%%

if (exist('fftconv_ola'))
  block_len = 2*256+1 - L2; % Make the FFT length a power of 2 (radix-2).
  eval('tic')
  Z3 = fftconv_ola(X, Y, block_len);
  t3= toc;
  fprintf('fftconv_ola: %f [s]\n', t3);
  fprintf('Numerical error: ||conv-fftconv_ola|| = %e\n\n', sum(sum((Z0-Z3).^2)));
end

%%
%% Test FFTW wisdoms
%%

[y,wisdom] = fftconv_p(X(:,1), Y(:,1));
eval('tic')
Z2_2 = fftconv_p(X, Y, wisdom);
t2 = toc;
fprintf('fftconv_p with wisdom: %f [s]\n',t2);
fprintf('Numerical error: ||conv-fftconv_p (w. wisdom)|| = %e\n\n',...
        sum(sum((Z0-Z2_2).^2)));

%%
%% In-place modes
%%

Z2 = zeros(size(Z2));
eval('tic')
fftconv_p(X, Y, Z2);
t2 = toc;
fprintf('fftconv_p in-place: %f [s]\n', t2);
fprintf('Numerical error: ||conv-fftconv_p (in-place)|| = %e\n\n',...
        sum(sum((Z0-Z2).^2)));

eval('tic')
fftconv_p(X, Y, Z2, '=');
t2 = toc;
fprintf('fftconv_p in-place = mode: %f [s]\n\n', t2);

eval('tic')
fftconv_p(X, Y, Z2, '+=');
t2 = toc;
fprintf('fftconv_p in-place += mode: %f [s]\n\n', t2);

eval('tic')
fftconv_p(X, Y, Z2, '-=');
t2 = toc;
fprintf('fftconv_p in-place -= mode: %f [s]\n\n', t2);

%%
%% sum_ffftconv
%%

h_len = 500;                    % Typically the SIR length.
u_len = 100;                    % Typically input signal length.
fft_len = h_len + u_len - 1;

L = 64;                      % Typically the number of array elements.
N = 500;                     % Typically number of observation points.
H = randn(h_len, N, L);      % Typically the Tx SRIs.
U = randn(u_len, L); % Typically the input signals for each array element.

eval('tic')
YF = zeros(fft_len,N);
for l=1:L
  for n=1:N
    YF(:,n) = YF(:,n) + fft(H(:,n,l), fft_len).* fft(U(:,l), fft_len);
  end
end
Y = real(ifft(YF))/fft_len; % Normalize w.r.t fft_len.
t2 = toc;
fprintf('Naive implementation of sum_fftconv: %f [s]\n\n', t2);

eval('tic')
Y2 = sum_fftconv(H, U);
t2 = toc;
fprintf('sum_fftconv: %f [s]\n', t2);
fprintf('Numerical error: ||Naive sum_fftconv - sum_fftconv|| = %e\n\n',...
        sum(sum((Y-Y2).^2)));

Y2 = zeros(fft_len,N);
eval('tic')
sum_fftconv(H, U, Y2);
t2 = toc;
fprintf('sum_fftconv in-place: %f [s]\n', t2);
fprintf('Numerical error: ||Naive sum_fftconv - sum_fftconv (in-place)|| = %e\n\n',...
        sum(sum((Y-Y2).^2)));
