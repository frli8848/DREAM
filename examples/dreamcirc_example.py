import math
import numpy as np
import matplotlib.pyplot as plt

import dreamcirc as dc

Fs = 4.0                        # Sampling freq. [MHz].
Ts = 1/Fs                       # [us].

# Depth
z = 10.0                        # [mm]

xo = np.linspace(0.0, 50.0, num=50, endpoint=True)
xo = gx.reshape(51,1)
[m,n] = np.shape(xo)
gy = np.zeros((m,1), dtype=float)
gz = z*np.ones((m,1), dtype=float)
Ro = gx
Ro = np.hstack((G, gy))
Ro = np.hstack((G, gz))

# Descretization parameters.
dx = 0.03                # [mm].
dy = 0.03                # [mm].
dt = Ts                  # [us].
nt = 400.0               # Length of spatial impulse response vector.
s_par = np.asmatrix([dx,dy,dt,nt])

# Material parameters.
v     = 1.0                    # Normal velocity.
cp    = 1500.0                 # Sound speed.
alpha = 0.0                    # Absorbtion (dB/cm Hz).
m_par = np.asmatrix([v,cp,alpha])

# Element size [mm].
R = 10.0
geom_par = np.asmatrix([R]);
H = dc.dreamcirc(Ro,geom_par,s_par,delay,m_par,"stop")

fig = plt.figure(1)
plt.clf()
plt.plot(t, H, 'b-', label='Attenuation alpha = 0 [dB/cm MHz]', linewidth=2, markersize=12)
plt.xlabel('t [us]')
plt.ylabel('h_{SIR} [m/s]')
plt.title('Circular transducer')
plt.legend()
plt.show()
