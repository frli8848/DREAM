#!/usr/bin/python

import math
import numpy as np
import matplotlib.pyplot as plt

#import sys
#sys.path.insert(0, '../')

import dream_arr_circ as dac

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

# Element size [mm].
R = 0.4;
geom_par = np.asmatrix([R]);

# Grid function (position vectors of the elements).
gx = np.linspace(-10.0, 10.0, num=21, endpoint=True)
gx = gx.reshape(21,1)
[num_elements,n] = np.shape(gx)
gy = np.zeros((num_elements,1), dtype=float);
gz = np.zeros((num_elements,1), dtype=float);
G = gx
G = np.hstack((G, gy))
G = np.hstack((G, gz))

# Focusing parameters.
#foc_met = "off";
foc_met = "x";
#foc_met = "y";
#foc_met = "xy";
#foc_met = "x+y";
focal = np.asmatrix(10.0);      # Focus depth

## Beam steering.
steer_met = "off";
#steer_met = "x";
#steer_met = "y";
#steer_met = "xy";

theta  = 0.0;                   # Angle in x-direction.
phi    = 0.0;                   # Angle in y-direction.
steer_par = np.asmatrix([theta,phi]);

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

H = dac.dream_arr_circ(Ro,geom_par,G,s_par,delay,m_par,foc_met,focal,steer_met,steer_par,apod_met,apod,win_par,"stop")

fig = plt.figure(1)
plt.clf()

plt.plot(t, H, 'b-', label='Attenuation alpha = 0 [dB/cm MHz]', linewidth=2, markersize=12)

alpha  = 5.0;                   # Absorbtion (dB/cm Hz).
m_par = np.asmatrix([v,cp,alpha]);
Hatt = dac.dream_arr_circ(Ro,geom_par,G,s_par,delay,m_par,foc_met,focal,steer_met,steer_par,apod_met,apod,win_par,"stop");
plt.plot(t, Hatt, 'r--',label='Attenuation alpha = 5 [dB/cm MHz]', linewidth=2, markersize=12)

plt.title('Array with circular elements')

plt.xlabel('t [us]')
plt.ylabel('h_{SIR} [m/s]')

plt.legend()
plt.show()
