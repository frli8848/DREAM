#!/usr/bin/python

import math
import numpy as np
import matplotlib.pyplot as plt

#import sys
#sys.path.insert(0, '../')

import dreamcirc as dc

#plt.rcParams['text.usetex'] = True

#setup_dream_parameters
Fs = 4.0;                       # Sampling freq. [MHz].
Ts = 1/Fs;                      # [us].

z = 10.0;

# One point.
xo = 0;
yo = 0;
zo = z;
ro = np.asmatrix([xo,yo,zo]);
Ro = ro;
print("Observation point (x,y,z) = (%f, %f, %f)" % (xo,yo,zo))

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
R = 10.0;
geom_par = np.asmatrix([R]);

H = dc.dreamcirc(Ro,geom_par,s_par,delay,m_par,"stop")
alpha  = 5.0;                   # Absorbtion (dB/cm Hz).
m_par = np.asmatrix([v,cp,alpha]);
Hatt = dc.dreamcirc(Ro,geom_par,s_par,delay,m_par,"stop");

if 'DO_PLOTTING' in locals():

    fig = plt.figure(1)
    plt.clf()
    plt.plot(t, H, 'b-', label='Attenuation alpha = 0 [dB/cm MHz]', linewidth=2, markersize=12)
    plt.plot(t, Hatt, 'r--',label='Attenuation alpha = 5 [dB/cm MHz]', linewidth=2, markersize=12)
    plt.xlabel('t [us]')
    plt.ylabel('h_{SIR} [m/s]')
    plt.title('Circular transducer')
    plt.legend()
    plt.show()
