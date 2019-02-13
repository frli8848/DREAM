/***
*
* Copyright (C) 2003,2006,2007,2008,2009,2012,2014 Fredrik Lingvall
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

// $Revision: 770 $ $Date: 2014-05-16 09:06:42 +0200 (Fri, 16 May 2014) $ $LastChangedBy: frli8848 $

#include <math.h>
#include <stdlib.h>

#include <octave/oct.h>
#include "fft.h"

#ifdef USE_FFTW
#include <fftw3.h>
#endif

#ifdef USE_FFTW

// Global plans.
fftw_plan    p_forward;
fftw_plan    p_backward;

extern "C"
void fft_init(dream_idx_type n, fftw_complex *xc, double *RESTRICT y)
{
  p_backward = fftw_plan_dft_c2r_1d(n, xc, y, FFTW_MEASURE);

  return;
}

extern "C"
void fft_close()
{
  fftw_destroy_plan(p_backward); 

  return;
}

#endif

/***
 *
 * Amplitude spectrum of a real vector. (Not used in DREAM)
 *
 ***/

extern "C"
void dream_fft(double *RESTRICT x, double *RESTRICT y, dream_idx_type n)
{
  dim_vector dims(n,1);
  NDArray xi(dims);
  ComplexNDArray yc;
  NDArray yr;
  dream_idx_type i;
  
  // Copy input vector to a Matlab Matrix.
  for (i=0; i<n; i++)
    xi(i) = x[i];
  
  // Let Octave do the job!
  yc = xi.fourier();
  yr = real(yc);

  for (i=0; i<n; i++)
    y[i] = yr(i);
  
  return;
}

/***
 *
 * Inverse DFT of a real vector (returns the real part). (Not used in DREAM)
 *
 ***/

extern "C"
void dream_ifft(double *RESTRICT x, double *RESTRICT y, dream_idx_type n)
{
  dim_vector dims(n,1);
  NDArray xi(dims);
  ComplexNDArray yc;
  NDArray yr;
  dream_idx_type i;
  
  // Copy input vector to a Octave Matrix.
  for (i=0; i<n; i++)
    xi(i) = x[i];
  
  // Let Octave do the job!
  yc = xi.ifourier();
  yr = real(yc);

  for (i=0; i<n; i++)
    y[i] = yr(i);
  
  return;
}

/***
 *
 * Inverse DFT of a complex vector (returns the real part).
 *
 ***/

#ifdef USE_FFTW
extern "C"
void cr_ifft(fftw_complex *xc, double *RESTRICT y, dream_idx_type n)
#else
extern "C"
void cr_ifft(double *RESTRICT xir, double *RESTRICT xii, double *RESTRICT y, dream_idx_type n)
#endif
{
  dream_idx_type k;

#ifdef USE_FFTW

  fftw_execute_dft_c2r(p_backward,xc,y);
  
  // Normalize.
  for (k=0; k<n; k++)
    y[k] /= (double) n;

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


#ifdef USE_FFTW
extern "C"
void cc_ifft(fftw_complex *xc, fftw_complex *yc, dream_idx_type n)
#else
extern "C"
void cc_ifft(double *RESTRICT xir, double *RESTRICT xii, double *RESTRICT yor, double *RESTRICT yoi, dream_idx_type n)
#endif
{

#ifdef USE_FFTW
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

