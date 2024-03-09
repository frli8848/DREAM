#################################################################
#
# 1D convolution algorithms
#
#################################################################

using LinearAlgebra
#using Plots
using PyPlot
using TickTock
using Printf
using DSP
using VectorizedRoutines

# DREAM functions
using fftconv_p_m
using dreamrect_m
using das_m

#DO_PLOTTING=1;

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
    #clf;
    #subplot(211)
    #plot(
    plot(t, h_e)
    xlabel("t [\\mu s]")
    title("System impulse response")

    #subplot(212)
    #f = (0:1023)/1024/Ts/2;
    #plot(f,20*log10(f_e/max(f_e)));
    #xlabel!("f [MHz]");
    #)
end

#
# RCA TFM data (full matrix capture - a.k.a FMC)
#

d  = 0.5;                       # Array pitch
xo = (-25:d:25);
yo = zeros(length(xo),1);
zo = z_pt*ones(length(xo),1);
Ro_t = [xo[:] yo[:] zo[:]];

# Crossed transmit and receve elemets.

# Geometrical parameters.
a = 0.4;                        # x-size.
b = 50;				# y-size.
geom_par_t = [a,b];

delay = [0.0];
Ht  = dreamrect_m.dreamrect(Ro_t,geom_par_t,s_par,delay,m_par,"stop");

geom_par_r = [b,a];
Ro_r = [yo[:] xo[:] zo[:]];
Hr = dreamrect_m.dreamrect(Ro_r,geom_par_r,s_par,delay,m_par,"stop");

L = length(xo);
Yfmc = zeros(nt+nt-1+nt_he-1,L^2);

# Loop over all transmit elements
n_t = 1;
for n=1:L:L^2

    # Julia is very picky about data types so we need
    # to make the args of fftconv_p a ::Matrix{Float64}!
    H_n_t = reshape(Ht[:,n_t],:,1);
    H_e = reshape(h_e,:,1);

    Hdp = fftconv_p_m.fftconv_p(Hr, H_n_t); # Double-path SIRs for the n_t:th transmit
    Yfmc[:,n:(n+L-1)] = fftconv_p_m.fftconv_p(Hdp, H_e);
    global n_t = n_t + 1;
end

Yfmc = Yfmc ./ maximum(maximum(abs.(Yfmc),dims=1),dims=1); # Normalize amplitudes

if (@isdefined(DO_PLOTTING))

    # NB. This plots so slow in Julia that we have switched it off!
    if (false)
        figure(2);
        #clf;
        t_dp = 0:Ts:(Ts*(size(Yfmc,1)-1));
        #imagesc(1:L^2, t_dp, Yfmc)
        pcolor(1:L^2, t_dp, Yfmc)
        title("FMC RCA B-scan");
        xlabel("A-scan index");
        ylabel("t [\\mu s]");
    end

    figure(3);
    #clf;
    t_dp = 0:Ts:(Ts*(size(Yfmc,1)-1));
    #imagesc(xo,t_dp,Yfmc[:,1:L:L^2])
    pcolor(xo, t_dp, Yfmc[:,1:L:L^2])
    title("FMC RCA B-scan when transmit element = receive element");
    xlabel("x [mm]");
    ylabel("t [\\mu s]");
end

num_elements = size(xo,1);
Gt = [xo[:] zeros(num_elements,1) zeros(num_elements,1)];
Gr = [zeros(num_elements,1) xo[:] zeros(num_elements,1)];

# Observation points for DAS
x = -25:0.5:25;
z = (0:63)/64*20; # Make sure its a factor of 64 (the OpenCL work group size).
X, Z = Matlab.meshgrid(x,z); # From VectorizedRoutines
Y = zeros(size(X));
Ro_rca = [X[:] Y[:] Z[:]];

delay = [system_delay]; # Compensate for the pulse/system (transducer) delay.
Im_rca = das_m.das(Yfmc, Gt, Gr, Ro_rca, dt, delay, cp,"rca","ignore","cpu");

if (@isdefined(DO_PLOTTING))
    figure(4);
    #clf;
    #imagesc(x,z,reshape(Im_rca,length(z),length(x)));
    pcolor(x, z, reshape(Im_rca,length(z),length(x)));
    title("RCA CPU Reconstruction");
    xlabel("x [mm]");
    ylabel("z [mm]");
end

Im_rca_gpu = das_m.das(Yfmc, Gt, Gr, Ro_rca, dt, delay, cp,"rca", "ignore","gpu");

if (@isdefined(DO_PLOTTING))
    figure(5);
    #clf;
    #images(x,z,reshape(Im_rca_gpu,length(z),length(x)));
    pcolor(x, z, reshape(Im_rca_gpu,length(z),length(x)));
    title("RCA GPU Reconstruction");
    xlabel("x [mm]");
    ylabel("z [mm]");
end
