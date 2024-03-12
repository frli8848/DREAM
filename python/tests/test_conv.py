#!/usr/bin/env python

#################################################################
#
# 1D convolution algorithms
#
#################################################################

import math
import numpy as np
import scipy as sp
#import matplotlib.pyplot as plt
import time

#import sys
#sys.path.insert(0, '../')

import fftconv_p as ft_p

L1 = 1000
L2 = 400
N = 50
X = np.random.randn(L1, N)
Y = np.random.randn(L2, N)

#
# Serial Time-domain
#

Z0 = np.zeros((L1+L2-1, N))
t = time.time()
for n in range(N):
    Z0[:,n] = sp.signal.convolve(X[:,n], Y[:,n])
t0 = time.time() - t
print("conv: %f [s]" % (t0))

#
# Threaded frequency-domain
#

# NB. The first call to ft_p.fftconv_p(X, Y) is slow.
t = time.time()
Z2 = ft_p.fftconv_p(X, Y)
t2 = time.time() - t
print("fftconv_p: %f [s]\n" %  (t2))
sse = np.sum( (Z0-Z2)**2)
print("Numerical error: ||conv-fftconv_p|| = %e\n" % (sse))
