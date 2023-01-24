#
# dreamsphere
#

using LinearAlgebra
using Plots

using dreamsphere_m

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


t = range(0.0, Ts*float(nt-1), length=Int(nt));

# FIXME: Make subplots instead.

# Focused

# Geometrical parameters.
R = 10;
Rcurv = 25.0;
geom_par = [R,Rcurv];

H = dreamsphere_m.dreamsphere(Ro,geom_par,s_par,delay,m_par,"stop");

plot(t, H, label="Focused: Attenuation alpha = 0 [dB/cm MHz]", linewidth=2)

alpha  = 5.0;                   # Absorbtion (dB/cm Hz).
m_par = [v,cp,alpha];
Hatt = dreamsphere_m.dreamsphere(Ro,geom_par,s_par,delay,m_par,"stop");

plot!(t, Hatt, label="Focused: Attenuation alpha = 5 [dB/cm MHz]", linewidth=2)

#title!("Focused spherical transducer at Rcurv=25 [mm]")
#xlabel!("t [us]")
#ylabel!("h_{SIR} [m/s]")

# Defocused

Rcurv = -25.0;
geom_par = [R,Rcurv];

alpha  = 0.0;                   # Absorbtion (dB/cm Hz).
m_par = [v,cp,alpha];

Hd = dreamsphere_m.dreamsphere(Ro,geom_par,s_par,delay,m_par,"stop");

plot!(t, Hd, label="Defocused: Attenuation alpha = 0 [dB/cm MHz]", linewidth=2)

alpha  = 5.0;                   # Absorbtion (dB/cm Hz).
m_par = [v,cp,alpha];
Hattd = dreamsphere_m.dreamsphere(Ro,geom_par,s_par,delay,m_par,"stop");

plot!(t, Hattd, label="Defocused: Attenuation alpha = 5 [dB/cm MHz]", linewidth=2)

title!("Spherical transducer with Rcurv = 25 or -25 [mm]")
xlabel!("t [us]")
ylabel!("h_{SIR} [m/s]")

gui()
