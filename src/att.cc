/***
 *
 * Copyright (C) 2008,2009,2015 Fredrik Lingvall
 *
 * This file is part of the DREAM Toolbox.
 *
 * The DREAM Toolbox is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * The DREAM Toolbox is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with the DREAM Toolbox; see the file COPYING.  If not, write to the 
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 *
 ***/

// $Revision: 855 $ $Date: 2015-05-05 13:57:35 +0200 (Tue, 05 May 2015) $ $LastChangedBy: frli8848 $

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "fft.h"
typedef std::complex<double> Complex;
#include "att.h"

// This file is just a C++ version of att.c for Octave. It's for
// enabling the toolbox to be build with MSVC (which isn't iso C99 
// so one cannot use <complex.h> and one therefore needs to use the 
// <complex> C++ class instead). 

#ifdef USE_FFTW

#ifdef PARALLEL
#include <thread>
#include <mutex>
#include <signal.h>
#endif

/***
 *
 * 1) Allocate nthreads working vectors xc and buf.
 *
 * 2) At each att call lock a semaphore and let the next thread use the next unlocked xc and buf, respectively.
 *
 * 3) unlock the semaphore when att(.) is done so that the buffers can be reused.
 *
 * At each call do:
 *
 * for n = 1 to nthreads
 *
 *    if XC[n] == locked
 *      do nothing (break)
 *    else
 *      xc = XC[n];
 *      buf = BUF[n];
 *      lock XC[n] and BUF[n].
 *   end
 * end 
 *
 ***/

//
// Globals
//

std::complex<double> **XC = NULL;
double       **BUF = NULL;
dream_idx_type nthreads = 1;

#ifdef PARALLEL
std::mutex **buffer_locks = NULL;
#endif

/***
 *
 * Allocate buffers and setup mutexes.
 *
 ***/

void att_init(dream_idx_type nt, dream_idx_type n_threads)
{
  dream_idx_type n;
  
  nthreads = n_threads;
  
  // Allocate memory.
#ifdef PARALLEL
  if (nthreads > 1)
    buffer_locks = (std::mutex**) malloc(nthreads*sizeof(std::mutex*));
#endif
  
  //XC = (fftw_complex**) fftw_malloc(nthreads*sizeof(fftw_complex*)); 
  XC = reinterpret_cast<std::complex<double>**>(fftw_malloc(nthreads*sizeof(fftw_complex*))); 
  BUF = (double**) fftw_malloc(nthreads*sizeof(double*)); 
  
  for (n=0; n<nthreads; n++) {
    XC[n] = reinterpret_cast<std::complex<double>*>(fftw_malloc((nt+1)*sizeof(fftw_complex))); 
    BUF[n] = (double*) fftw_malloc((nt+1)*sizeof(double));
    
#ifdef PARALLEL
    if (nthreads > 1)
      buffer_locks[n] = new std::mutex();
#endif
  }
  
  // Init FFTW plan etc.
  fft_init(nt, reinterpret_cast<fftw_complex*>(XC[0]),BUF[0]);
  
  return;
}

/***
 *
 * Free allocated buffer memory and mutexes.
 *
 ***/

void att_close()
{
  dream_idx_type n;
  
  for (n=0; n<nthreads; n++) {
    fftw_free(XC[n]);
    fftw_free(BUF[n]);
    
#ifdef PARALLEL
    if (nthreads > 1)
      delete &(buffer_locks[n]);
#endif
  }
  
  fftw_free(XC);
  fftw_free(BUF);

#ifdef PARALLEL
  if (nthreads > 1)
    free(buffer_locks);
#endif
  
  fft_close();
  
  return;
}

#endif

/***
 *
 * Impulse response for absorption.
 *
 ***/

void att(double alfa, double rj, dream_idx_type it, double dt, double cp, double *h, dream_idx_type nt, double ai)
{
  const double mille = 1000.0;
  double pi,pi2;
  double a0,a1,b;
  dream_idx_type i, k;
  double b1, b2, x1;
  double dw;
  double w_n, t;
  double w;
  double Fs;

  pi = 4.0 * atan(1.0);
  pi2 = pi*pi;

#ifdef USE_FFTW
  
  //fftw_complex *xc = NULL;
  std::complex<double> *xc = NULL;
  double          *buf = NULL;
  int             this_buffer = 0;
  
  // pthread_mutex_trylock behaves identically to pthread_mutex_lock, except
  // that it does not block the calling  thread  if  the  mutex  is  already
  // locked  by  another  thread  (or by the calling thread in the case of a
  // ``fast'' mutex).  Instead,  pthread_mutex_trylock  returns  immediately
  // with the error code EBUSY.
  
  if (nthreads > 1) {
    // Use the first set of buffers that are not locked.
    for(k=0; k<nthreads; k++) {
      
#ifdef PARALLEL
      if (buffer_locks[k]->try_lock())
	buffer_locks[k]->lock();
      
#endif
      this_buffer = k;
      xc = XC[k];
      buf = BUF[k];
#ifdef PARALLEL
    }
#endif
  }
  
  else { // nthreads = 1.
    xc = XC[0];
    buf = BUF[0];
  }
  
#else
  
  double *xr,*xi,*buf;
  
  // Allocate memory.
  xr  = (double*) malloc((nt+1)*sizeof(double)); 
  xi  = (double*) malloc((nt+1)*sizeof(double)); 
  buf = (double*) malloc((nt+1)*sizeof(double));
  
  // Clear the vectors (really not needed).
  for (i=0; i<=nt; i++) {
    xr[i]  = (double) 0.0;
    xi[i]  = (double) 0.0;
    buf[i] = (double) 0.0;
  }

#endif
  
  //
  // Calculate frequency-domain Green's function.
  //

  /* absor = alfa en 1/m Hz: */
  //alfa /= (double) 8.86*10000.0; // dB per cm MHz to Neper per m Hz conversion? (8.686 = 20/log(10)).
  alfa /= (double) 8.686*10000.0; 

  /*  	dt s, a1 s, */
  dt /= (mille * mille);	// Sampling period [s].
  Fs  = 1/dt;			// Sampling frequecy [Hz].
  dw  = (double) 2.0 * pi / (double) nt; // Angular freq. sampling step [rad].
  t   = it * dt;

  /* Change  units rj [m] */
  rj /= mille;
  a0  = alfa;
  // The 0.95 constant controls the phase only (causality). See:
  // K. Aki and P. G. Richards, "Quantative Seismology: Theory and Methods",
  // San Francisco, CA, Freeman, 1980.
  a1  = dt * (double) 0.95 / pi; 
  //a1  = dt / pi; 
  
#ifdef USE_FFTW
  xc[0] = Complex(1.0,0.0);	// w = 0 (f = 0 Hz).
#else
  xr[0] = (double) 1.0;		// w = 0 (f = 0 Hz).
  xi[0] = (double) 0.0;		// w = 0 (f = 0 Hz).
#endif

  for (k=1; k<(nt/2+1); k++) { // Loop over all freqs.
    
    w_n = ((double) k) * dw;	// Normalized angular freq. 
    w   = w_n*Fs;		// Angular freq.
    x1  = exp(-rj*a0 * w/(2*pi)); // Amplitude of the Green's function.
    b   = log ( ((double) 1.0) / (a1 * w));
    
    // (Roughly) Eq.(A6) in Piwakowski and Sbai, IEEE UFFC, vol 46,
    // No 2, March 1999, p. 422--440. 
    b1 = cos( -w*rj*a0/(pi2)*b - t*w ); // Real part.
    b2 = sin( -w*rj*a0/(pi2)*b - t*w ); // Imag part.
    
#ifdef USE_FFTW
    xc[k]    = Complex(x1*b1,x1*b2);
    xc[nt-k] = Complex(x1*b1,-x1*b2);
    //xc[nt-k] = conj(xc[k]);

    // Slower than above.
    //xc[k] = cexp(-I*(w*rj*a0/pi2*b - t*w));
    //xc[nt-k] = conj(xc[k]);
    
#else
    // Real part.    
    xr[k]    = x1 * b1;
    xr[nt-k] = x1 * b1;
    
    // Imag part.
    xi[k]    =  x1 * b2;
    xi[nt-k] = -x1 * b2;
#endif
  }

  //
  // Do inverse Fourier transform to get time-domain Green's function.
  //

#ifdef USE_FFTW
  cr_ifft(reinterpret_cast<fftw_complex*>(xc),buf,nt); // Complex input, real output.
#else
  cr_ifft(xr,xi,buf,nt); // Complex input, real output.
#endif
  
  for (i=0; i<nt; i++) {
    h[i] += ai * buf[i];
  }

#ifdef USE_FFTW

#ifdef PARALLEL
  // We are done. Unlock the buffers
  if (nthreads > 1)
    buffer_locks[this_buffer]->unlock();
#endif

#else
  free(xr);
  free(xi);
  free(buf);  
#endif

  return;
} /* att */



/***
 *
 * Impulse response for absorption - annular array.
 *
 ***/

void att_annu(double alfa, double rj, dream_idx_type it, double  dt, double cp, double *h, dream_idx_type nt, double ai, int ns, int isize)
{
  /* Initialized data */
  const double mille = 1000.0;
  double pi;
  double a,b;
  dream_idx_type i, k;
  double a1;
  double b1, b2, b3, b4, x1;
  double dw;
  double w_n, tq, t;
  double w;
  
  pi = 4.0 * atan(1.0);


#ifdef USE_FFTW

  //fftw_complex *xc = NULL;
  std::complex<double> *xc = NULL;
  double          *buf = NULL;
  int             this_buffer = 0;
  
  if (nthreads > 1) {
    // Use the first set of buffers that are not locked.
    for(k=0; k<nthreads; k++) {

#ifdef PARALLEL
      if (buffer_locks[k]->try_lock()) {
	buffer_locks[k]->lock();
#endif
	this_buffer = k;
	xc = XC[k];
	buf = BUF[k];
#ifdef PARALLEL
      }
#endif
    }
  }
  else { // nthreads = 1.
    xc = XC[0];
    buf = BUF[0];
  }

#else

  double *xr,*xi,*buf;

  // Allocate memory.
  xr  = (double*) malloc((nt+1)*sizeof(double)); 
  xi  = (double*) malloc((nt+1)*sizeof(double)); 
  buf = (double*) malloc((nt+1)*sizeof(double));

  // Clear the vectors (really not needed).
  for (i=0; i<=nt; i++) {
    xr[i]  = (double) 0.0;
    xi[i]  = (double) 0.0;
    buf[i] = (double) 0.0;
  }

#endif
  
  //
  // Calculate frequency-domain Green's function.
  //

  // dB per cm MHz to Neper per m Hz conversion? (8.686 = 20/log(10)).
  alfa /= (double) 8.686*10000.0; 

  /*  	dt s, a1 s, */
  dt /= (mille * mille);     // Sampling period [s].
  // The 0.95 constant controls the phase only (causality). See:
  // K. Aki and P. G. Richards, "Quantative Seismology: Theory and Methods",
  // San Francisco, CA, Freeman, 1980.
  a1  = dt * (double) 0.95 / pi; 
  //a1  = dt / pi; 
  dw  = (double) 2.0 * pi / (double) nt; // Freq. sampling step.

  t = it * dt;

  /* Change  units rj [m] */
  rj /= mille;
  a = -(rj * alfa) / ( (double) 2.0 * pi);
  tq = rj * alfa * cp / (pi * pi);

#ifdef USE_FFTW
  xc[0] = Complex(1.0,0.0);	// w = 0 (f = 0 Hz).
#else
  xr[0] = (double) 1.0;		// w = 0 (f = 0 Hz).
  xi[0] = (double) 0.0;		// w = 0 (f = 0 Hz).
#endif
  for (k=1; k<(nt/2+1); k++) {
    
    w_n = (double) k * dw;
    w = w_n / dt;
    x1 = exp(a * w);    
    b = log(1.0 / (a1 * w));
    
    // t temp arrive en sec, ftb1 retard de focal/baleyage en sec 
    b4 = -t * w_n / dt;
    b3 = -(b * tq) * w_n / (dt * cp);
    
    b = b4 + b3;
    b1 = cos(b);
    b2 = sin(b);

#ifdef USE_FFTW
    xc[k]    = Complex(x1*b1,x1*b2);
    xc[nt-k] = Complex(x1*b1, -x1*b2);
#else
    // Real part.    
    xr[k]    = x1 * b1;
    xr[nt-k] = x1 * b1;
    
    // Imag part.
    xi[k]    =  x1 * b2;
    xi[nt-k] = -x1 * b2;
#endif
  }

  // Should these be set to zero ?
  //xr[nt/2] = 0.0;
  //xi[nt/2] = 0.0;
  
  //
  // Do inverse Fourier transform to get time-domain Green's function.
  //

#ifdef USE_FFTW
  cr_ifft(reinterpret_cast<fftw_complex*>(xc),buf,nt); // Complex input, real output.
#else
  cr_ifft(xr,xi,buf,nt); // Complex input, real output.
#endif
  
  for (i=0; i<nt; i++) {
    h[i + ns*nt] += ai * buf[i];
  }
#ifdef USE_FFTW

#ifdef PARALLEL
  // We are done. Unlock the buffers
  if (nthreads > 1) 
    buffer_locks[this_buffer]->unlock();
#endif

#else
  free(xr);
  free(xi);
  free(buf);  
#endif

  return;
}
