#
# dreamrect
#

using LinearAlgebra
using Plots
using LazyGrids

using dreamrect_m

Fs = 4.0;                       # Sampling freq. [MHz].
Ts = 1/Fs;                      # [us].

z = 10.0;

# One point.
xo = 0;
yo = 0;
zo = z;
ro = [xo yo zo];
Ro = ro;
#print("Observation point (x,y,z) = (%f, %f, %f)" % (xo,yo,zo))

# Descretization parameters.
dx = 0.05;                # [mm].
dy = 0.05;                # [mm]
dt = Ts;                  # [us].
nt = 300.0;               # Length of spatial impulse response vector.
s_par = [dx,dy,dt,nt];

# Material parameters.
v     = 1.0;                    # Normal velocity.
cp    = 1000.0;                 # Sound speed.
alpha = 0.0;                    # Absorbtion [dB/cm Hz].
m_par = [v,cp,alpha];

#t_z = z*1.0e-3/cp * 1.0e6;      # us
delay = [0.0];

# Element size [mm].
a = 10.0;
b = 15.0;
geom_par = [a,b];


# Make sure the number of obs points are a 64 (the OpenCL workgroup size)
x = range(-15, 15, 61)
y = range(-15, 15, 61)
z = range(1.25, 65, 256)
(X,Y,Z) = (xg, yg, zg) = ndgrid(x, y, z);
Ro = [vec(X) vec(Y) vec(Z)];
size(Ro)

# This disables Julia's (REPL) signal handler for SIGINT
# so that we can CTRL-C and interrupt our code.
ccall(:jl_exit_on_sigint, Cvoid, (Cint,), 0)

@time H = dreamrect_m.dreamrect(Ro,geom_par,s_par,delay,m_par,"stop");
