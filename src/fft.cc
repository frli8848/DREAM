/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2013,2014,2016 Fredrik Lingvall
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

#include <uchar.h>
#include "mex.h"
#include "fft.h"

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

#ifdef HAVE_FFTW

// Global plans.
fftw_plan    p_forward;		// Double precision.
fftw_plan    p_backward;

// Double precision.
void fft_init(dream_idx_type n, fftw_complex *xc, double *y)
{
  p_backward = fftw_plan_dft_c2r_1d(n, xc, y, FFTW_MEASURE);

  return;
}

// Double precision.
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

// Double precision.

#ifdef HAVE_FFTW
void cr_ifft(fftw_complex *xc, double *y, dream_idx_type n)
#else
void cr_ifft(double *xir, double *xii, double *y, dream_idx_type n)
#endif
{

#ifdef HAVE_FFTW
  dream_idx_type k;

  fftw_execute_dft_c2r(p_backward,xc,y);

  // Normalize.
  for (k=0; k<n; k++)
    y[k] /= (double) n;

#else

  mxArray *X, *Y;
  double *yr, *yi, *xr, *xi;
  dream_idx_type k;

  X = mxCreateDoubleMatrix(n,1,mxCOMPLEX);
  mxComplexDouble *x = mxGetComplexDoubles(X);

  // Copy input vector to a Matlab Matrix.
  for (k=0; k<n; k++) {
    xr[k] = x[k].real;
    xi[k] = x[k].imag;
  }

  // Let Matlab do the job!
  mexCallMATLAB(1, &Y, 1, &X, "ifft");
  mxComplexDouble *yo = mxGetComplexDoubles(Y);

  // Use the real part only.
  for (k=0; k<n; k++) {
    y[k] = yo[k].real;
  }

  mxDestroyArray(X);
  mxDestroyArray(Y);

#endif
  return;
}

/***
 *
 * Inverse DFT of a complex vector (not used in DREAM).
 *
 ***/

#ifdef HAVE_FFTW
void cc_ifft(fftw_complex *xc, fftw_complex *yc, dream_idx_type n)
#else
void cc_ifft(double *xir, double *xii, double *yor, double *yoi, dream_idx_type n)
#endif
{

#ifdef HAVE_FFTW

  dream_idx_type k;

  fftw_execute_dft(p_backward,xc,yc);

  // Normalize.
  std::complex<double> *y = reinterpret_cast<std::complex<double>*>(yc);
  std::complex<double> len_inv(1.0/double(n), 0.0);
  for (k=0; k<n; k++) {
    y[k] *= len_inv;
  }

#else

  dream_idx_type i;
  mxArray *X,*Y;
  double *yr, *yi, *xr, *xi;

  X = mxCreateDoubleMatrix(n,1,mxCOMPLEX);
  mxComplexDouble *x = mxGetComplexDoubles(X);
  //xr = mxGetPr(X);
  //xi = mxGetPi(X);

  // Copy input vector to a Matlab Matrix.
  for (i=0; i<n; i++) {
    //xr[i] = xir[i];
    //xi[i] = xii[i];
    x[i].real = xir[i];
    x[i].imag = xii[i];
  }

  // Let Matlab do the job!
  mexCallMATLAB(1, &Y, 1, &X, "ifft");
  mxComplexDouble *y = mxGetComplexDoubles(Y);
  //yr = mxGetPr(Y);
  //yi = mxGetPi(Y);

  // Copy data.
  for (i=0; i<n; i++) {
    //yor[i] = yr[i];
    //yoi[i] = yi[i];
    yor[i] = y[i].real;
    yoi[i] = y[i].imag;
  }

  mxDestroyArray(X);
  mxDestroyArray(Y);

#endif

  return;
}
