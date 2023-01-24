#!/usr/bin/python

import time
import math
import numpy as np
import matplotlib.pyplot as plt

#import sys
#sys.path.insert(0, '../')

import dreamrect as dr

#plt.rcParams['text.usetex'] = True

#setup_dream_parameters
Fs = 4.0;                       # Sampling freq. [MHz].
Ts = 1/Fs;                      # [us].

# Descretization parameters.
dx = 0.05;                # [mm].
dy = 0.05;                # [mm]
dt = Ts;                  # [us].
nt = 300.0;               # Length of spatial impulse response vector.
s_par = np.asmatrix([dx,dy,dt,nt]);

t = np.linspace(0.0, Ts*float(nt-1), num=int(nt), endpoint=True)

# Material parameters.
v     = 1.0;                    # Normal velocity.
cp    = 1000.0;                 # Sound speed.
alpha = 0.0;                    # Absorbtion [dB/cm Hz].
m_par = np.asmatrix([v,cp,alpha]);

#t_z = z*1.0e-3/cp * 1.0e6;      # us
delay = np.asmatrix(0.0);

# Element size [mm].
a = 10.0;
b = 15.0;
geom_par = np.asmatrix([a,b]);

# Make sure the number of obs points are a 64 (the OpenCL workgroup size)
x = np.linspace(-15, 15, 61)
y = np.linspace(-15, 15, 61)
z = np.linspace(1.25, 65, 256)
X,Y,Z = np.meshgrid(x, y, z, indexing='ij');

M = np.shape(X)[0]
N = np.shape(X)[1]
K = np.shape(X)[2]
Ro_len = M*N*K
print("Num obs points: %d" % Ro_len)

X = np.reshape(X, (M*N*K, 1), order='F');
Y = np.reshape(Y, (M*N*K, 1), order='F');
Z = np.reshape(Z, (M*N*K, 1), order='F');

Ro = X
Ro = np.hstack((Ro, Y))
Ro = np.hstack((Ro, Z))

t = time.time()
H = dr.dreamrect(Ro,geom_par,s_par,delay,m_par,"stop")
elapsed = time.time() - t
print("dreamrect elapesed time: %f" % elapsed)

#tic,[H,err] = dreamrect(Ro,geom_par,s_par,delay,m_par,'stop','gpu');toc
