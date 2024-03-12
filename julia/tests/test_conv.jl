#################################################################
#
# 1D convolution algorithms
#
#################################################################

using LinearAlgebra
#using Plots
using TickTock
using Printf
using DSP

using conv_p_m
using fftconv_p_m

L1 = 1000;
L2 = 400;
N = 50;
X = randn(L1,N);
Y = randn(L2,N);

#
# Serial Time-domain
#

Z0 = zeros(L1+L2-1,N);
tick()
for n=1:N
  Z0[:,n] = conv( X[:,n], Y[:,n]);
end
t0= tok();
@printf("conv: %f [s]\n\n", t0);

#
# Threaded time-domain
#

tick()
Z1 = conv_p_m.conv_p(X, Y)
t1 = tok()
@printf("conv_p: %f [s]\n", t1)
@printf("Numerical error: ||conv-conv_p|| = %e\n\n", sum(sum((Z0-Z1).^2)))

#
# Threaded frequency-domain
#

tick()
Z2 = fftconv_p_m.fftconv_p(X, Y)
t2 = tok()
@printf("fftconv_p: %f [s]\n", t2)
@printf("Numerical error: ||conv-fftconv_p|| = %e\n\n", sum(sum((Z0-Z2).^2)))
