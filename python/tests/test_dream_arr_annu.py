#!/usr/bin/python

import math
import numpy as np
import matplotlib.pyplot as plt

#import sys
#sys.path.insert(0, '../')

import dream_arr_annu as daa

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
alpha = 0.0;                    # Absorbtion (dB/cm Hz).
m_par = np.asmatrix([v,cp,alpha]);

#t_z = z*1.0e-3/cp * 1.0e6;      # us
delay = np.asmatrix(0.0);

# Grid function (position vectors of the elements).
gr = np.asmatrix([3.5, 3.6, 5.0, 5.1, 6.0, 6.1, 7.0, 7.1, 8.0])
gr = gr.reshape(9,1)
[num_elements,n] = np.shape(gr)
G = gr

# Focusing parameters.
#foc_met = "off";
foc_met = "x";
#foc_met = "y";
#foc_met = "xy";
#foc_met = "x+y";
focal = np.asmatrix(25.0);      # Focus depth

## Apodization.
apod_met = "off";
#apod_met = "ud";                       # User defined.
#apod_met = "triangle";
#apod_met = "gauss";
#apod_met = "raised";                   # Raised cosine.
#apod_met = "simply";                   # Simply supported.
#apod_met = "clamped";
apod = np.ones((num_elements,1));     # Apodization weights for "ud".
win_par = 1.0; # Parameter for raised cos and Gaussian apodization functions.

H = daa.dream_arr_annu(Ro,G,s_par,delay,m_par,foc_met,focal,apod_met,apod,win_par,"stop")

alpha  = 5.0;                   # Absorbtion (dB/cm Hz).
m_par = np.asmatrix([v,cp,alpha]);
Hatt = daa.dream_arr_annu(Ro,G,s_par,delay,m_par,foc_met,focal,apod_met,apod,win_par,"stop");

if 'DO_PLOTTING' in locals():

    fig = plt.figure(1)
    plt.clf()
    plt.plot(t, H, 'b-', label='Attenuation alpha = 0 [dB/cm MHz]', linewidth=2, markersize=12)

    plt.plot(t, Hatt, 'r--',label='Attenuation alpha = 5 [dB/cm MHz]', linewidth=2, markersize=12)
    plt.title('Annular Array')
    plt.xlabel('t [us]')
    plt.ylabel('h_{SIR} [m/s]')
    plt.legend()
    plt.show()
