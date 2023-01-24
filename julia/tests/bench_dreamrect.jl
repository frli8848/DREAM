#
# dreamrect
#

using LinearAlgebra
using Plots
using Printf
using LazyGrids

using dreamrect_m

Fs = 4.0;                       # Sampling freq. [MHz].
Ts = 1/Fs;                      # [us].

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
x = range(-15, 15, 61);
y = range(-15, 15, 61);
z = range(1.25, 65, 256);
(X,Y,Z) = ndgrid(x, y, z); # Same as meshgrid
Ro = [vec(X) vec(Y) vec(Z)];
(num_obs_points, n) = size(Ro);
@printf("Num obs points: %d\n", num_obs_points);

# This disables Julia's (REPL) signal handler for SIGINT
# so that we can CTRL-C and interrupt our code if we want.
ccall(:jl_exit_on_sigint, Cvoid, (Cint,), 0)

@time H = dreamrect_m.dreamrect(Ro,geom_par,s_par,delay,m_par,"stop");
