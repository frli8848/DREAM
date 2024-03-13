#!/usr/bin/env python

import math
import numpy as np
import matplotlib.pyplot as plt

#import sys
#sys.path.insert(0, '../')

import dreamrect as dr
import das as das
import das_f as das_f
import fftconv_p as ft_p

#
# ------------- Delay-and-sum --------------------------
#

Fs = 50.0                      # Sampling freq. [MHz].
Ts = 1/Fs                      # [us].

# Descretization parameters.
dx = 0.05                # [mm].
dy = 0.05                # [mm]
dt = Ts                  # [us].
nt = 2000                # Length of spatial impulse response vector.
s_par = np.asmatrix([dx,dy,dt,nt])

#t = np.linspace(0.0, Ts*float(nt-1), num=int(nt), endpoint=True)

# Material parameters.
v     = 1.0                    # Normal velocity.
cp    = 1500                   # Sound speed.
alpha  = 0.0                   # Absorbtion (dB/cm Hz).
m_par = np.asmatrix([v,cp,alpha])

# Simulate a single point scatterer at (x=0,z=10)
z_pt = 10

#
# Generate simulated data
#

# Simulated electrical impulse response.

nt_he = 150
t = np.linspace(0.0, Ts*float(nt_he-1), num=int(nt_he), endpoint=True)

f0 = 2.5                             # Center frequency [MHz].
t0 = 0.55                            # Time delay to max amplitude [us].
a_n = 10                             # Envelop parameter.

system_delay = t0 + 0.21 # Delay to the max of the pulse.

h_e = - np.exp(-a_n * (t - t0)** 2) * np.cos(2*np.pi*f0 * t)

print("f = %1.2f [MHz]" % (f0))
lmbd = cp/f0/1.0e3 # [mm].
print("lambda = %1.2f [mm]" % (lmbd))

# f_e = abs(freqz(h_e,1,1024))
#h_e = h_e/max(f_e) # Unity gain at center freq.

if 'DO_PLOTTING' in locals():

    fig = plt.figure(1)
    plt.clf()

    #subplot(211)
    plt.plot(t, h_e)
    plt.xlabel("t [\\mu s]")
    plt.title('System impulse response')

    #subplot(212)
    #f = (0:1023)/1024/Ts/2
    #plot(f,20*log10(f_e/max(f_e)))
    #xlabel!("f [MHz]")
    #)

    #plt.legend()
    plt.show(block=False)

#
# RCA TFM data (full matrix capture - a.k.a FMC)
#

d  = 0.5                        # Array pitch
#xo = (-25:d:25)
xo = np.linspace(-25.0,25.0, num=int((50.0/d)+1), endpoint=True)
yo = np.zeros((xo.shape[0], 1))
zo = z_pt * np.ones((xo.shape[0], 1))

# We cannot append column-wise in python as in MATLAB/Octave so
# we have to do it row-wise and transpose.
Ro_t = np.asmatrix([xo.flatten('F'), yo.flatten('F'), zo.flatten('F')])
Ro_t = Ro_t.T

# Crossed transmit and receve elemets.

# Geometrical parameters.
a = 0.4                         # element x-size [mm].
b = 50                          # element y-size [mm].
geom_par_t = np.asmatrix([a,b])

delay = np.asmatrix([0.0])
Ht = dr.dreamrect(Ro_t,geom_par_t,s_par,delay,m_par,"stop")

geom_par_r = np.asmatrix([b,a])

Ro_r = np.asmatrix([yo.flatten('F'), xo.flatten('F'), zo.flatten('F')])
Ro_r = Ro_r.T

Hr = dr.dreamrect(Ro_r,geom_par_r,s_par,delay,m_par,"stop")

L = xo.shape[0]
Yfmc = np.zeros((nt+nt-1+nt_he-1, L*L))

# Loop over all transmit elements
n_t = 0
for n in range(0, L*L, L) :
    Hdp = ft_p.fftconv_p(Hr, Ht[:,n_t]) # Double-path SIRs for the n_t:th transmit
    #print("n = %d n_t = %d M = %d x M = %d" % (n, n_t, Hdp.shape[0], Hdp.shape[1]))
    Yfmc[:,n:(n+L)] = ft_p.fftconv_p(Hdp, h_e)
    n_t += 1

Yfmc = Yfmc / np.max(np.max(np.abs(Yfmc),axis=0),axis=0) # Normalize amplitudes

if 'DO_PLOTTING' in locals():

    t_dp = np.linspace(0.0, Ts*Yfmc.shape[0], num=int(Yfmc.shape[0]), endpoint=True)

    if False : # NB. This is plots very slow
        fig = plt.figure(2)
        plt.clf()
        a_scan_idx = np.linspace(1.0, Yfmc.shape[1], num=int(Yfmc.shape[1]), endpoint=True)
        plt.pcolor(a_scan_idx, t_dp, Yfmc)
        plt.gca().invert_yaxis()
        plt.title("FMC RCA B-scan");
        plt.xlabel("A-scan index");
        plt.ylabel("t [\\mu s]");
        plt.show(block=False)

    fig = plt.figure(3)
    plt.clf()
    plt.pcolor(xo, t_dp, Yfmc[:,1:L*L:L])
    plt.gca().invert_yaxis()
    plt.title("FMC RCA B-scan");
    plt.xlabel("A-scan index");
    plt.ylabel("t [\\mu s]");
    plt.show(block=False)

num_elements = xo.shape[0]

# Transmit element positions
xt = xo.flatten('F')
yt = np.zeros((num_elements,1)).flatten('F')
zt = np.zeros((num_elements,1)).flatten('F')
Gt = np.asmatrix([xt,yt,zt])
Gt = Gt.T

# Receive element positions
xr = np.zeros((num_elements,1)).flatten('F')
yr = xo.flatten('F')
zr = np.zeros((num_elements,1)).flatten('F')
Gr = np.asmatrix([xr,yr,zr])
Gr = Gr.T

# Observation points for DAS
x = np.linspace(-25.0,25.0, num=int((50.0/d)+1), endpoint=True)
z = np.linspace(0.0,20.0, num=int(64), endpoint=True) # Make sure its a multiple of 64 (the OpenCL work group size).

X, Z = np.meshgrid(x,z)
Y = np.zeros((X.shape[0], X.shape[1]))
Ro_rca = np.asmatrix([X.flatten('F'), Y.flatten('F'), Z.flatten('F')])
Ro_rca = Ro_rca.T

delay = np.matrix(system_delay) # Compensate for the pulse/system (transducer) delay.

#
# CPU - double precision
#

Im_rca = das.das(Yfmc, Gt, Gr, Ro_rca, dt, delay, cp, "rca", "ignore", "cpu")

if 'DO_PLOTTING' in locals():

    fig = plt.figure(4)
    plt.clf()
    py_Im = np.reshape(Im_rca, (x.shape[0], z.shape[0]))
    plt.pcolor(x, z, py_Im.T)
    plt.gca().invert_yaxis()
    plt.title("RCA CPU Reconstruction")
    plt.xlabel("x [mm]")
    plt.ylabel("z [mm]")
    plt.show(block=False)

#
# GPU - double precision
#

Im_rca_gpu = das.das(Yfmc, Gt, Gr, Ro_rca, dt, delay, cp, "rca", "ignore", "gpu")

if 'DO_PLOTTING' in locals():

    fig = plt.figure(5)
    plt.clf()
    py_Im = np.reshape(Im_rca_gpu, (x.shape[0], z.shape[0]))
    plt.pcolor(x, z, py_Im.T)
    plt.gca().invert_yaxis()
    plt.title("RCA GPU Reconstruction")
    plt.xlabel("x [mm]")
    plt.ylabel("z [mm]")
    plt.show(block=False)


#
# CPU - single precision
#

# Convert args to single precision.
Yfmc_f   = np.float32(Yfmc)
Gt_f     = np.float32(Gt)
Gr_f     = np.float32(Gr)
Ro_rca_f = np.float32(Ro_rca)
dt_f     = np.float32(dt)
delay_f  = np.float32(delay)
cp_f     = np.float32(cp)

Im_rca_f  = das_f.das(Yfmc_f, Gt_f, Gr_f, Ro_rca_f, dt_f, delay_f, cp_f, "rca", "ignore", "cpu")

if 'DO_PLOTTING' in locals():

    fig = plt.figure(6)
    plt.clf()
    py_Im = np.reshape(Im_rca_f, (x.shape[0], z.shape[0]))
    plt.pcolor(x, z, py_Im.T)
    plt.gca().invert_yaxis()
    plt.title("RCA GPU  Single Precision Reconstruction")
    plt.xlabel("x [mm]")
    plt.ylabel("z [mm]")
    plt.show(block=False)

#
# GPU - single precision
#

Im_rca_gpu_f  = das_f.das(Yfmc_f, Gt_f, Gr_f, Ro_rca_f, dt_f, delay_f, cp_f, "rca", "ignore", "gpu")

if 'DO_PLOTTING' in locals():

    fig = plt.figure(7)
    plt.clf()
    py_Im = np.reshape(Im_rca_gpu_f, (x.shape[0], z.shape[0]))
    plt.pcolor(x, z, py_Im.T)
    plt.gca().invert_yaxis()
    plt.title("RCA GPU  Single Precision Reconstruction")
    plt.xlabel("x [mm]")
    plt.ylabel("z [mm]")
    plt.show(block=False)
