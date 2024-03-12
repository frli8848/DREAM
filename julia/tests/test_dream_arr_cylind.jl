#
# dream_arr_cylind
#

using LinearAlgebra
using Plots

using dream_arr_cylind_m

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
a = 0.8;
b = 20.0;
Rcurv = 10.0;
geom_par = [a,b,Rcurv];

# Grid function (position vectors of the elements).
num_elements = 21;
gx = range(-10.0, 10, length=num_elements);
gx = vec(gx)
gy = zeros(num_elements,1);
gz = zeros(num_elements,1);
G = [gx gy gz];

# Focusing parameters.
#foc_met = "off";
#foc_met = "x";
#foc_met = "y";
foc_met = "xy";
#foc_met = "x+y";
focal = [10.0];

## Beam steering.
steer_met = "off";
#steer_met = "x";
#steer_met = "y";
#steer_met = "xy";

theta  = 0.0;                   # Angle in x-direction.
phi    = 0.0;                   # Angle in y-direction.
steer_par = [theta,phi];

## Apodization.
apod_met = "off";
#apod_met = "ud";                       # User defined.
#apod_met = "triangle";
#apod_met = "gauss";
#apod_met = "raised";                   # Raised cosine.
#apod_met = "simply";                   # Simply supported.
#apod_met = "clamped";
apod = vec(ones(num_elements,1)); # Apodization weights for "ud". NB Must be a "1D" array (vector).
win_par = 1.0; # Parameter for raised cos and Gaussian apodization functions.

# Focused

H = dream_arr_cylind_m.dream_arr_cylind(Ro,geom_par,G,s_par,delay,m_par,foc_met,focal,steer_met,steer_par,apod_met,apod,win_par,"stop")

t = range(0.0, Ts*float(nt-1), length=Int(nt));

alpha  = 5.0;                   # Absorbtion (dB/cm Hz).
m_par = [v,cp,alpha];
Hatt = dream_arr_cylind_m.dream_arr_cylind(Ro,geom_par,G,s_par,delay,m_par,foc_met,focal,steer_met,steer_par,apod_met,apod,win_par,"stop");

if (@isdefined(DO_PLOTTING))

    plot(t, H, label="Focused: Attenuation alpha = 0 [dB/cm MHz]", linewidth=2)
    plot!(t, Hatt, label="Focused: Attenuation alpha = 5 [dB/cm MHz]", linewidth=2)

    #title!("Array with rectangular elements")
    #xlabel!("t [us]")
    #ylabel!("h_{SIR} [m/s]")
end

# Defocused

Rcurv = -10.0;
geom_par = [a,b,Rcurv];

alpha  = 0.0;                   # Absorbtion (dB/cm Hz).
m_par = [v,cp,alpha];

Hd = dream_arr_cylind_m.dream_arr_cylind(Ro,geom_par,G,s_par,delay,m_par,foc_met,focal,steer_met,steer_par,apod_met,apod,win_par,"stop")

alpha  = 5.0;                   # Absorbtion (dB/cm Hz).
m_par = [v,cp,alpha];

Hattd = dream_arr_cylind_m.dream_arr_cylind(Ro,geom_par,G,s_par,delay,m_par,foc_met,focal,steer_met,steer_par,apod_met,apod,win_par,"stop");

if (@isdefined(DO_PLOTTING))

    plot!(t, Hd, label="Defocused: Attenuation alpha = 0 [dB/cm MHz]", linewidth=2)
    plot!(t, Hattd, label="Defocused: Attenuation alpha = 5 [dB/cm MHz]", linewidth=2)

    title!("Array with cylindrical transducer elements with Rcurv = 10.0 and -10.0")
    xlabel!("t [us]")
    ylabel!("h_{SIR} [m/s]")

    gui()
end
