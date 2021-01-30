/***
*
* Copyright (C) 2003,2006,2007,2008,2009,2012,2014,2021 Fredrik Lingvall
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


#include <math.h>
#include <stdlib.h>

#include <octave/oct.h>
#include "fft.h"

#ifdef HAVE_FFTW

#include <fftw3.h>

// Global plans.
fftw_plan    p_forward;
fftw_plan    p_backward;

void fft_init(dream_idx_type n, fftw_complex *xc, double *y)
{
  p_backward = fftw_plan_dft_c2r_1d(n, xc, y, FFTW_MEASURE);

  return;
}

void fft_close()
{
  fftw_destroy_plan(p_backward);

  return;
}

#endif

/***
 *
 * Inverse DFT of a complex vector (returns the real part).
 *
 ***/

#ifdef HAVE_FFTW
void cr_ifft(fftw_complex *xc, double *y, dream_idx_type n)
#else
void cr_ifft(double *xir, double *xii, double *y, dream_idx_type n)
#endif
{
  dream_idx_type k;

#ifdef HAVE_FFTW

  fftw_execute_dft_c2r(p_backward,xc,y);

  // Normalize.
  for (k=0; k<n; k++) {
    y[k] /= (double) n;
  }

#else

  ComplexMatrix xc(n,1);
  ComplexMatrix yc;

  // Copy input vector to an Octave Matrix.
  for (k=0; k<n; k++) {
    xc(k) = Complex(xir[k],xii[k]);
  }

  // Let Octave do the job!
  yc = xc.ifourier();

  // Use the real part only.
  for (k=0; k<n; k++)
    y[k] = yc(k).real();

#endif

  return;
}

/***
 *
 * Inverse DFT of a complex vector. (Not used in DREAM)
 *
 ***/


#ifdef HAVE_FFTW
void cc_ifft(fftw_complex *xc, fftw_complex *yc, dream_idx_type n)
#else
void cc_ifft(double *xir, double *xii, double *yor, double *yoi, dream_idx_type n)
#endif
{

#ifdef HAVE_FFTW
  ; // Nothing here yet.
#else
  dream_idx_type k;

  ComplexMatrix xc(n,1);
  ComplexMatrix yc;

  // Copy input vector to a Ocatve Matrix.
  for (k=0; k<n; k++) {
    xc(k) = Complex(xir[k],xii[k]);
  }

  // Let Octave do the job!
  yc = xc.ifourier();

  // Copy data.
  for (k=0; k<n; k++) {
    yor[k] = yc(k).real();
    yoi[k] = yc(k).imag();
  }

#endif

  return;
}
