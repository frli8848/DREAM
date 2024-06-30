#!/usr/bin/env python

import math
import numpy as np
import matplotlib.pyplot as plt
import time

#import sys
#sys.path.insert(0, '../')

# Check if we have vtk installed
try:
    import vtk
    have_vtk =  True

except:
    have_vtk =  False

import dreamrect as dr
import das as das
import das_f as das_f
import fftconv_p as ft_p

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
s_par = np.asmatrix([dx,dy,dt,nt])

# Material parameters.
v     = 1.0;                    # Normal velocity.
cp    = 1500;                   # Sound speed.
alpha  = 0.0;                   # Absorbtion (dB/cm Hz).
m_par = np.asmatrix([v,cp,alpha])

# Simulate a single point scatterer at (x=0,z=10)
z_pt = 10;

#
# Generate simulated data
#

# Simulated electrical impulse response.

nt_he = 150;
t = np.linspace(0.0, Ts*float(nt_he-1), num=int(nt_he), endpoint=True)

f0 = 2.5;                             # Center frequency [MHz].
t0 = 0.55;                            # Time delay to max amplitude [us].
a_n = 10;                             # Envelop parameter.

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

#
# 3D Observation points for DAS RCA
#

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
y = np.linspace(-25.0,25.0, num=int((50.0/d)+1), endpoint=True)
z = np.linspace(0.0,20.0, num=int(64), endpoint=True) # Make sure its a multiple of 64 (the OpenCL work group size).

X,Y,Z = np.meshgrid(x,y,z)
Ro_rca = np.asmatrix([X.flatten('F'), Y.flatten('F'), Z.flatten('F')])
Ro_rca = Ro_rca.T

Nx = x.shape[0]
Ny = y.shape[0]
Nz = z.shape[0]

delay = np.matrix(system_delay) # Compensate for the pulse/system (transducer) delay.

t = time.time()
Im_rca = das.das(Yfmc, Gt, Gr, Ro_rca, dt, delay, cp, "rca", "ignore", "cpu")
t1 = time.time() - t
print("das RCA CPU: %f [s]\n" %  (t1))

if 'DO_PLOTTING' in locals():

    fig = plt.figure(2)
    plt.clf()

    O_cpu = np.reshape(Im_rca, (Nz, Nx*Ny))
    O_cpu = O_cpu.T
    c_scan_cpu = np.reshape(np.max(np.abs(O_cpu),axis=1), (Nx, Ny))
    mx = np.max(np.max(c_scan_cpu,axis=0),axis=0);

    #imagesc(x, y, 20.0*log10(c_scan_cpu/mx));
    plt.pcolor(x, y, 20.0*np.log10(c_scan_cpu/mx))
    #h_cb = colorbar;
    #h_cb_title = get(h_cb,'Title');
    #set(h_cb_title,'String','Normalized Amplitude [dB]')
    #axis square;
    plt.xlabel('x [mm]');
    plt.ylabel('y [mm]');
    plt.title('C-scan RCA DAS beamformed data');
    plt.show(block=False)


t = time.time()
Im_rca_gpu = das.das(Yfmc, Gt, Gr, Ro_rca, dt, delay, cp, "rca", "ignore", "gpu")
t2 = time.time() - t
print("\ndas RCA GPU: %f [s]\n" %  (t2))

if 'DO_PLOTTING' in locals():

    fig = plt.figure(3);
    plt.clf()

    O_gpu = np.reshape(Im_rca_gpu, (Nz, Nx*Ny));
    O_gpu = O_gpu.T
    c_scan_gpu = np.reshape(np.max(np.abs(O_gpu),axis=1), (Nx, Ny));
    mx = np.max(np.max(c_scan_gpu,axis=0),axis=0);

    #imagesc(x, y, 20.0*log10(c_scan_gpu/mx));
    plt.pcolor(x, y, 20.0*np.log10(c_scan_gpu/mx))

    #h_cb = colorbar;
    #h_cb_title = get(h_cb,'Title');
    #set(h_cb_title,'String','Normalized Amplitude [dB]')
    #axis square;
    plt.xlabel('x [mm]');
    plt.ylabel('y [mm]');
    plt.title('C-scan RCA GPU DAS beamformed data');
    plt.show(block=False)

#
# Single precision
#

print('*** Single precision ***');

# Convert args to single precision.
Yfmc_f   = np.float32(Yfmc)
Gt_f     = np.float32(Gt)
Gr_f     = np.float32(Gr)
Ro_rca_f = np.float32(Ro_rca)
dt_f     = np.float32(dt)
delay_f  = np.float32(delay)
cp_f     = np.float32(cp)

t = time.time()
Im_rca_f  = das_f.das(Yfmc_f, Gt_f, Gr_f, Ro_rca_f, dt_f, delay_f, cp_f, "rca", "ignore", "cpu")
t3 = time.time() - t
print("das single precision RCA CPU: %f [s]\n" %  (t3))

if 'DO_PLOTTING' in locals():

    fig = plt.figure(4)
    plt.clf()

    O_cpu = np.reshape(Im_rca_f, (Nz, Nx*Ny))
    O_cpu = O_cpu.T
    c_scan_cpu = np.reshape(np.max(np.abs(O_cpu),axis=1), (Nx, Ny))
    mx = np.max(np.max(c_scan_cpu,axis=0),axis=0);

    #imagesc(x, y, 20.0*log10(c_scan_cpu/mx));
    plt.pcolor(x, y, 20.0*np.log10(c_scan_cpu/mx))
    #h_cb = colorbar;
    #h_cb_title = get(h_cb,'Title');
    #set(h_cb_title,'String','Normalized Amplitude [dB]')
    #axis square;
    plt.xlabel('x [mm]');
    plt.ylabel('y [mm]');
    plt.title('C-scan RCA CPU DAS Single Precision beamformed data');
    plt.show(block=False)

t = time.time()
Im_rca_gpu_f  = das_f.das(Yfmc_f, Gt_f, Gr_f, Ro_rca_f, dt_f, delay_f, cp_f, "rca", "ignore", "gpu")
t4 = time.time() - t
print("\ndas single precision RCA GPU: %f [s]\n" %  (t4))

if 'DO_PLOTTING' in locals():

    fig = plt.figure(5)
    plt.clf()

    O_cpu = np.reshape(Im_rca_gpu_f, (Nz, Nx*Ny))
    O_cpu = O_cpu.T
    c_scan_cpu = np.reshape(np.max(np.abs(O_cpu),axis=1), (Nx, Ny))
    mx = np.max(np.max(c_scan_cpu,axis=0),axis=0);

    #imagesc(x, y, 20.0*log10(c_scan_cpu/mx));
    plt.pcolor(x, y, 20.0*np.log10(c_scan_cpu/mx))
    #h_cb = colorbar;
    #h_cb_title = get(h_cb,'Title');
    #set(h_cb_title,'String','Normalized Amplitude [dB]')
    #axis square;
    plt.xlabel('x [mm]');
    plt.ylabel('y [mm]');
    plt.title('C-scan RCA GPU DAS Single Precision beamformed data');
    plt.show(block=False)

Ntot = Nx*Ny*Nz*num_elements*num_elements
print("* Total number of DAS operations: %d x %d x %d x %d x %d = %d" % (Nx,Ny,Nz,num_elements,num_elements, Ntot))
print("* CPU double precision %1.2f [MDAS/s] " % (Ntot/t1/1e6))
print("* GPU double precision %1.2f [MDAS/s] " % (Ntot/t2/1e6))
print("* CPU single precision %1.2f [MDAS/s] " % (Ntot/t3/1e6))
print("* GPU single precision %1.2f [MDAS/s] " % (Ntot/t4/1e6))

#
# 3D VTK plot
#

if have_vtk:

    #maxNumPoints = 1e6
    zMin = -0.05
    zMax = 0.05
    vtkPolyData = vtk.vtkPolyData()

    #clearPoints()
    vtkPoints = vtk.vtkPoints()
    vtkCells = vtk.vtkCellArray()
    vtkDepth = vtk.vtkDoubleArray()
    vtkDepth.SetName('DepthArray')
    vtkPolyData.SetPoints(vtkPoints)
    vtkPolyData.SetVerts(vtkCells)
    vtkPolyData.GetPointData().SetScalars(vtkDepth)
    vtkPolyData.GetPointData().SetActiveScalars('DepthArray')

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(vtkPolyData)
    mapper.SetColorModeToDefault()
    mapper.SetScalarRange(zMin, zMax)
    mapper.SetScalarVisibility(1)
    vtkActor = vtk.vtkActor()
    vtkActor.SetMapper(mapper)

    # Add data
    for k in range(np.size(Im_rca, 0)):
        point = [Ro_rca[k,0],Ro_rca[k,1],Ro_rca[k,2]]
        #point = 20.0*(np.random.rand(3)-0.5)
        pointId = vtkPoints.InsertNextPoint(point[:])
        vtkDepth.InsertNextValue(Im_rca[k]/mx)
        vtkCells.InsertNextCell(1)
        vtkCells.InsertCellPoint(pointId)
        vtkCells.Modified()
        vtkPoints.Modified()
        vtkDepth.Modified()

    # Renderer
    renderer = vtk.vtkRenderer()
    renderer.AddActor(vtkActor)
    #renderer.SetBackground(.2, .3, .4)
    renderer.SetBackground(0.0, 0.0, 0.0)
    renderer.ResetCamera()

    # Render Window
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)

    # Interactor
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    # Begin Interaction
    renderWindow.Render()
    renderWindow.SetWindowName("XYZ Data Viewer")
    renderWindowInteractor.Start()
