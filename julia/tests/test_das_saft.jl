#################################################################
#
# 1D convolution algorithms
#
#################################################################

using LinearAlgebra
using TickTock
using Printf
using DSP
using VectorizedRoutines

if (@isdefined(DO_PLOTTING))
    #using Plots
    using PyPlot
end

# DREAM functions
using fftconv_p_m
using dreamrect_m
using das_m
using das_f_m

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

# Geometrical parameters.
a = 0.4;                        # x-size.
b = 15;				# y-size.
geom_par = [a,b];

#
# SAFT data
#

d  = 0.5;                       # Array pitch
xo = (-25:d:25);
yo = zeros(length(xo),1);
zo = z_pt*ones(length(xo),1);
Ro = [xo[:] yo[:] zo[:]];

delay = [0.0];

#H, err = dreamrect_m.dreamrect(Ro,geom_par,s_par,delay,m_par,"stop");
H = dreamrect_m.dreamrect(Ro,geom_par,s_par,delay,m_par,"stop");

# Julia is very picky about data types so we need
# to make the args of fftconv_p a ::Matrix{Float64}!
H_e = reshape(h_e,:,1);

Hdp = fftconv_p_m.fftconv_p(H, H); # Double-path SIRs
Ysaft = fftconv_p_m.fftconv_p(Hdp, H_e);
Ysaft = Ysaft ./ maximum(maximum(abs.(Ysaft),dims=1),dims=1); # Normalize amplitudes

if (@isdefined(DO_PLOTTING))
    figure(2);
    #clf;
    t_dp = 0:Ts:(Ts*(size(Ysaft,1)-1));
    #imagesc(xo,t_dp,Ysaft)
    pcolor(xo, t_dp, Ysaft)
    title("SAFT B-scan");
    xlabel("x [mm]");
    ylabel("t [\\mu s]");
end

num_elements = size(xo,1);
Gt = [xo[:] zeros(num_elements,1) zeros(num_elements,1)];
Gr = zeros(0,0); # Not used in SAFT mode.

# Observation points for DAS
x = -25:0.5:25;
z = (0:63)/64*20; # Make sure its a factor of 64 (the OpenCL work group size).
X, Z = Matlab.meshgrid(x,z); # From VectorizedRoutines
Y = zeros(size(X));
Ro_saft = [X[:] Y[:] Z[:]];

delay = [system_delay]; # Compensate for the pulse/system (transducer) delay.

#
# CPU - double precision
#

Im_saft = das_m.das(Ysaft, Gt, Gr, Ro_saft, dt, delay, cp,"saft", "ignore", "cpu");

if (@isdefined(DO_PLOTTING))
    figure(3);
    #clf;
    #imagesc(x,z,reshape(Im_saft,length(z),length(x)))
    pcolor(x, z, reshape(Im_saft,length(z),length(x)))
    title("SAFT CPU Reconstruction");
    xlabel("x [mm]");
    ylabel("z [mm]");
end

#
# GPU - double precision
#

Im_saft_gpu = das_m.das(Ysaft, Gt, Gr, Ro_saft, dt, delay, cp, "saft", "ignore", "gpu");

if (@isdefined(DO_PLOTTING))
    figure(4);
    #clf;
    #images(x,z,reshape(Im_saft_gpu,length(z),length(x)))
    pcolor(x, z, reshape(Im_saft_gpu,length(z),length(x)))
    title("SAFT GPU Reconstruction")
    xlabel("x [mm]")
    ylabel("z [mm]")
end

#
# CPU - single precision
#

# Convert args to single precision.
Ysaft_f   = convert(Matrix{Float32}, Ysaft)
Gt_f     = convert(Matrix{Float32}, Gt)
Gr_f     = convert(Matrix{Float32}, Gr)
Ro_saft_f = convert(Matrix{Float32}, Ro_saft)
dt_f     = convert(Float32, dt)
delay_f  = convert(Vector{Float32}, delay)
cp_f     = convert(Float32, cp)

Im_saft_f = das_f_m.das(Ysaft_f, Gt_f, Gr_f, Ro_saft_f, dt_f, delay_f, cp_f, "saft", "ignore", "cpu");

if (@isdefined(DO_PLOTTING))
    figure(3);
    #clf;
    #imagesc(x,z,reshape(Im_saft,length(z),length(x)))
    pcolor(x, z, reshape(Im_saft_f, length(z), length(x)))
    title("SAFT CPU Single Precision Reconstruction");
    xlabel("x [mm]");
    ylabel("z [mm]");
end

#
# GPU - single precision
#

Im_saft_gpu_f = das_f_m.das(Ysaft_f, Gt_f, Gr_f, Ro_saft_f, dt_f, delay_f, cp_f, "saft", "ignore", "gpu");

if (@isdefined(DO_PLOTTING))
    figure(4);
    #clf;
    #images(x,z,reshape(Im_saft_gpu,length(z),length(x)))
    pcolor(x, z, reshape(Im_saft_gpu_f, length(z), length(x)))
    title("SAFT GPU Single Precision Reconstruction")
    xlabel("x [mm]")
    ylabel("z [mm]")
end
