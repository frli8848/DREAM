/***
 *
 * Copyright (C) 2008,2009,2015,2021 Fredrik Lingvall
 *
 * Ths file is part of the DREAM Toolbox.
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


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <complex>
#include "attenuation.h"

// This file is just a C++ version of att.c for Octave. It's for
// enabling the toolbox to be build with MSVC (which isn't iso C99
// so one cannot use <complex.h> and one therefore needs to use the
// <complex> C++ class instead).

#if defined DREAM_OCTAVE || defined HAVE_FFTW

#include <thread>
#include <mutex>
#include <signal.h>

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
std::mutex **buffer_locks = NULL;

#endif

void Attenuation::att(FFTCVec &xc_vec, FFTVec &x_vec, double r, dream_idx_type it, double *h, double ai, dream_idx_type element)
{
  double a0,a1,b;
  dream_idx_type i, k;
  double b1, b2, x1;
  double dw;
  double w_n, t;
  double w;
  double Fs;

  std::complex<double> *xc = xc_vec.get();
  double *x = x_vec.get();

  //
  // Calculate frequency-domain Green's function.
  //

  // Absorbtion : alpha [dB/m Hz]
  // dB per cm MHz to Neper per m Hz conversion? (8.686 = 20/log(10)).
  double alpha = m_alpha / (8.686*10000.0);

  double dt = m_dt*1.0e-6;          // Sampling period [s].
  Fs  = 1/dt;                       // Sampling frequecy [Hz].
  dw  = 2.0 * M_PI / double(m_len); // Angular freq. sampling step [rad].
  t   = double(it) * dt;            // [s]
  r *= 1.0e-3;                      // [m]

  // The 0.95 constant controls the phase only (causality). See:
  // K. Aki and P. G. Richards, "Quantative Seismology: Theory and Methods",
  // San Francisco, CA, Freeman, 1980.
  a1  = dt * 0.95 / M_PI;

  xc[0] = std::complex<double>(1.0, 0.0); // w = 0 (f = 0 Hz).

  for (k=1; k<(m_len/2+1); k++) { // Loop over all freqs.

    w_n = double(k) * dw;	// Normalized angular freq.
    w   = w_n*Fs;		// Angular freq.
    x1  = exp(-r*alpha * w/(2.0*M_PI)); // Amplitude of the Green's function.
    b   = log (1.0/(a1 * w));

    // (Roughly) Eq.(A6) in M_Piwakowski and Sbai, IEEE UFFC, vol 46,
    // No 2, March 1999, p. 422--440.
    b1 = cos( -w*r*alpha/(M_PI*M_PI)*b - t*w ); // Real part.
    b2 = sin( -w*r*alpha/(M_PI*M_PI)*b - t*w ); // Imag part.

    xc[k]       = std::complex<double>(x1*b1, x1*b2);
    xc[m_len-k] = std::complex<double>(x1*b1,-x1*b2);

    // Slower than above.
    //xc[k] = cexp(-I*(w*r*a0*b/(M_PI*M_PI) - t*w));
    //xc[nt-k] = conj(xc[k]);
  }

  // Do inverse Fourier transform to get time-domain Green's function.
  m_fft->ifft(xc_vec, x_vec);

  if (element>0) {
    for (i=0; i<m_len; i++) {
      h[i+element*m_len] += ai * x[i];
    }
  } else {
    for (i=0; i<m_len; i++) {
      h[i] += ai * x[i];
    }
  }

  return;
}
