#
# dream_arr_annu
#

using LinearAlgebra
using Plots

using dream_arr_annu_m

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

# Grid function (position vectors of the elements).
gr = [3.5, 3.6, 5.0, 5.1, 6.0, 6.1, 7.0, 7.1, 8.0];
gr = vec(gr);
num_elements = 9;
G = gr

# Focusing parameters.
#foc_met = "off";
#foc_met = "x";
#foc_met = "y";
foc_met = "xy";
#foc_met = "x+y";
focal = [25.0];

# Apodization.
apod_met = "off";
#apod_met = "ud";                       # User defined.
#apod_met = "triangle";
#apod_met = "gauss";
#apod_met = "raised";                   # Raised cosine.
#apod_met = "simply";                   # Simply supported.
#apod_met = "clamped";
apod = vec(ones(num_elements,1));     # Apodization weights for "ud".
win_par = 1.0; # Parameter for raised cos and Gaussian apodization functions.

H = dream_arr_annu_m.dream_arr_annu(Ro,G,s_par,delay,m_par,foc_met,focal,apod_met,apod,win_par,"stop");

t = range(0.0, Ts*float(nt-1), length=Int(nt));

plot(t, H, label="Attenuation alpha = 0 [dB/cm MHz]", linewidth=2)

alpha  = 5.0;                   # Absorbtion (dB/cm Hz).
m_par = [v,cp,alpha];
Hatt = dream_arr_annu_m.dream_arr_annu(Ro,G,s_par,delay,m_par,foc_met,focal,apod_met,apod,win_par,"stop");

plot!(t, Hatt, label="Attenuation alpha = 5 [dB/cm MHz]", linewidth=2)

title!("Annular Array")
xlabel!("t [us]")
ylabel!("h_{SIR} [m/s]")

gui()
