#!/usr/bin/python

import math
import numpy as np
import matplotlib.pyplot as plt

#import sys
#sys.path.insert(0, '../')

import dreamsphere as dsp

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
nt = 300.0;                 # Length of spatial impulse response vector.
s_par = np.asmatrix([dx,dy,dt,nt]);

t = np.linspace(0.0, Ts*float(nt-1), num=int(nt), endpoint=True)

# Material parameters.
v     = 1.0;                    # Normal velocity.
cp    = 1000.0;                 # Sound speed.
alpha = 0.0;                    # Absorbtion (dB/cm Hz).
m_par = np.asmatrix([v,cp,alpha]);

#t_z = z*1.0e-3/cp * 1.0e6;      # us
delay = np.asmatrix(0.0);

# Focused

# Geometrical parameters.
R = 10;
Rcurv = 25.0;
geom_par = np.asmatrix([R,Rcurv]);

H = dsp.dreamsphere(Ro,geom_par,s_par,delay,m_par,"stop")

fig, axs = plt.subplots(2)

axs[0].plot(t, H, 'b-', label='Attenuation alpha = 0 [dB/cm MHz]', linewidth=2, markersize=12)

alpha  = 5.0;                   # Absorbtion (dB/cm Hz).
m_par = np.asmatrix([v,cp,alpha]);
Hatt = dsp.dreamsphere(Ro,geom_par,s_par,delay,m_par,"stop");
axs[0].plot(t, Hatt, 'r--',label='Attenuation alpha = 5 [dB/cm MHz]', linewidth=2, markersize=12)

axs[0].set_title('Focused spherical transducer at Rcurv=25 [mm]')
axs[0].legend()

axs[0].set_xlabel('t [us]')
axs[0].set_ylabel('h_{SIR} [m/s]')

# Defocused

Rcurv = -25.0;
geom_par = np.asmatrix([R,Rcurv]);

H = dsp.dreamsphere(Ro,geom_par,s_par,delay,m_par,"stop")

alpha  = 0.0;                   # Absorbtion (dB/cm Hz).
m_par = np.asmatrix([v,cp,alpha]);

H = dsp.dreamsphere(Ro,geom_par,s_par,delay,m_par,"stop")
axs[1].plot(t, H, 'b-', label='Attenuation alpha = 0 [dB/cm MHz]', linewidth=2, markersize=12)

alpha  = 5.0;                   # Absorbtion (dB/cm Hz).
m_par = np.asmatrix([v,cp,alpha]);
Hatt = dsp.dreamsphere(Ro,geom_par,s_par,delay,m_par,"stop");
axs[1].plot(t, Hatt, 'r--',label='Attenuation alpha = 5 [dB/cm MHz]', linewidth=2, markersize=12)

axs[1].set_title('Defocused spherical transducer at Rcurv=-25 [mm]')
axs[1].legend()

axs[1].set_xlabel('t [us]')
axs[1].set_ylabel('h_{SIR} [m/s]')


plt.show()
