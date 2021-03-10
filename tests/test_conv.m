%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1D convolution algorithms
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L1 = 10000;
L2 = 4000;
N = 1000;
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
fprintf('conv\n');

%%
%% Threaded time-domain
%%

eval('tic')
Z1 = conv_p(X, Y);
t1= toc;
fprintf('conv_p\n');

%%
%% Threaded frequency-domain
%%

eval('tic')
Z2 = fftconv_p(X, Y);
t2= toc;
fprintf('fftconv_p\n');

%%
%% Threaded frequency-domain using the overlap-and-add method
%%

block_len = 2*4096+1 - L2; % Make the FFT length a power of 2 (radix-2).
eval('tic')
Z3 = fftconv_ola(X, Y, block_len);
t3= toc;
fprintf('fftconv_ola\n\n');

fprintf('Elapsed time conv: %f, conv_p: %f, fftconv_p: %f, fftconv_ola: %f\n\n',t0,t1,t2,t3);
fprintf('Numerical error: ||conv-conv_p|| = %e, ||conv-fftconv_p|| = %e ||conv-fftconv_ola|| = %e\n\n',...
        sum(sum((Z0-Z1).^2)),...
        sum(sum((Z0-Z2).^2)),...
        sum(sum((Z0-Z3).^2)));
