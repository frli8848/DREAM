#
# dreamcirc_f
#

using LinearAlgebra
using Plots

using dreamcirc_f_m

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
R = 10.0;
geom_par = [R];

# Focusing parameters.
#foc_met = "off";
#foc_met = "x";
#foc_met = "y";
foc_met = "xy";
#foc_met = "x+y";
focal = [10.0];

H = dreamcirc_f_m.dreamcirc_f(Ro,geom_par,s_par,delay,m_par,foc_met,focal,"stop");

t = range(0.0, Ts*float(nt-1), length=Int(nt));

plot(t, H, label="Attenuation alpha = 0 [dB/cm MHz]", linewidth=2)

alpha  = 5.0;                   # Absorbtion (dB/cm Hz).
m_par = [v,cp,alpha];
Hatt = dreamcirc_f_m.dreamcirc_f(Ro,geom_par,s_par,delay,m_par,foc_met,focal,"stop");

plot!(t, Hatt, label="Attenuation alpha = 5 [dB/cm MHz]", linewidth=2)

title!("Focused circular transducer")
xlabel!("t [us]")
ylabel!("h_{SIR} [m/s]")

gui()
