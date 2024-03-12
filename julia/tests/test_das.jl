#################################################################
#
# 1D convolution algorithms
#
#################################################################

using LinearAlgebra
using Plots
using TickTock
using Printf
using DSP

# DREAM functions
using fftconv_p_m
using dreamrect_m
using das_m

#
# ------------- Delay-and-sum --------------------------
#

Fs = 50.0;                      # Sampling freq. [MHz].
Ts = 1/Fs;                      # [us].

# Descretization parameters.
dx = 0.05;                # [mm].
dy = 0.05;                # [mm]
dt = Ts;                  # [us].
nt = 2000;                # Length of spatial impulse response vector.
s_par = [dx,dy,dt,nt];

t = 0:Ts:Ts*(nt-1);

# Material parameters.
v     = 1.0;                    # Normal velocity.
cp    = 1500;                   # Sound speed.
alpha  = 0.0;                   # Absorbtion (dB/cm Hz).
m_par = [v,cp,alpha];

# Simulate a single point scatterer at (x=0,z=10)
z_pt = 10;

#
# Generate simulated data
#

# Simulated electrical impulse response.

nt_he = 150;
t = (0:((nt_he-1)))*Ts;

f0 = 2.5;                             # Center frequency [MHz].
t0 = 0.55;                            # Time delay to max amplitude [us].
a_n = 10;                             # Envelop parameter.

system_delay = t0 + 0.21; # Delay to the max of the pulse.

h_e = .- exp.( .- a_n .* (t .- t0) .^ 2) .* cos.(2 .* pi .* f0 .* t);

@printf("\nf = %1.2f [MHz]\n", f0);
lambda = cp/f0/1e3; # [mm].
@printf("lambda = %1.2f [mm]\n",lambda);

# f_e = abs(freqz(h_e,1,1024));
#z = PolynomialRatio(1, h_e)
#f_e, w = freqresp(z);
#h_e = h_e/max(f_e); # Unity gain at center freq.

if (@isdefined(DO_PLOTTING))
  figure(1);
  clf;
  subplot(211)
  plot(t, h_e);
  xlabel("t [\\mu s]")
  title("System impulse response")

  subplot(212)
  f = (0:1023)/1024/Ts/2;
  plot(f,20*log10(f_e/max(f_e)));
  xlabel("f [MHz]")
end

# Geometrical parameters.
a = 0.4;                        # x-size.
b = 15;				# y-size.
geom_par = [a,b];

#
# TFM data (full matrix capture - a.k.a FMC)
#

d  = 0.5;                       # Array pitch
xo = (-25:d:25);
yo = zeros(length(xo),1);
zo = z_pt*ones(length(xo),1);
Ro = [xo[:] yo[:] zo[:]];

delay = [0.0];

#H, err = dreamrect_m.dreamrect(Ro,geom_par,s_par,delay,m_par,"stop");
H = dreamrect_m.dreamrect(Ro,geom_par,s_par,delay,m_par,"stop");

L = length(xo);
Yfmc = zeros(nt+nt-1+nt_he-1,L^2);

# Loop over all transmit elements
n_t = 1;
for n=1:L:L^2

    # Julia is very picky about data types so we need
    # to make the args of fftconv_p a ::Matrix{Float64}!
    H_n_t = reshape(H[:,n_t],:,1);
    H_e = reshape(h_e,:,1);

    Hdp = fftconv_p_m.fftconv_p(H, H_n_t); # Double-path SIRs for the n_t:th transmit
    Yfmc[:,n:(n+L-1)] = fftconv_p_m.fftconv_p(Hdp, H_e);
    global n_t = n_t + 1;
end

Yfmc = Yfmc ./ maximum(maximum(abs.(Yfmc),dims=1),dims=1); # Normalize amplitudes

if (@isdefined(DO_PLOTTING))
  figure(2);
  clf;
  t_dp = 0:Ts:Ts .* (size(Yfmc,1)-1);
  imagesc(1:L^2,t_dp,Yfmc)
  title("FMC B-scan")
  xlabel("A-scan index")
  ylabel("t [\\mu s]")

  figure(3);
  clf;
  t_dp = 0:Ts:Ts*(size(Yfmc,1)-1);
  imagesc(xo,t_dp,Yfmc(:,1:L:L^2))
  title("FMC B-scan when transmit element = receive element")
  xlabel("x [mm]")
  ylabel("t [\\mu s]")
end

num_elements = size(xo,1);
Gt = [xo[:] zeros(num_elements,1) zeros(num_elements,1)];
Gr = Gt;

# Observation points for DAS
x = -25:0.5:25;
z = (0:63)/64*20; # Make sure its a factor of 64 (the OpenCL work group size).
X, Z = meshgrid(x,z);
Y = zeros(size(X));
Ro_tfm = [X[:] Y[:] Z[:]];

delay = system_delay; # Compensate for the pulse/system (transducer) delay.
Im_tfm = das(Yfmc, Gt, Gr, Ro_tfm, dt, delay, cp,"tfm");

if (@isdefined(DO_PLOTTING))
  figure(4);
  clf;
  imagesc(x,z,reshape(Im_tfm,length(z),length(x)))
  title("TFM Reconstruction")
  xlabel("x [mm]")
  ylabel("z [mm]")
end

Im_tfm_gpu = das(Yfmc, Gt, Gr, Ro_tfm, dt, delay, cp,"tfm", "ignore","gpu");

if (@isdefined(DO_PLOTTING))
  figure(5);
  clf;
  imagesc(x,z,reshape(Im_tfm_gpu,length(z),length(x)))
  title("TFM Reconstruction")
  xlabel("x [mm]")
  ylabel("z [mm]")
end
