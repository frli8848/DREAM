/***
*
* Copyright (C) 2021 Fredrik Lingvall
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

/***
 *
 * Frequency domain convolution of two vectors.
 *
 ***/

#include "fftconv.h"

void fftconv(FFT &fft,
             double *xr, octave_idx_type nx, double *yr, octave_idx_type ny, double *zr,
             double *a,  double *b, double *c,
             std::complex<double> *af, std::complex<double> *bf, std::complex<double> *cf,
             ConvMode conv_mode)
{
  octave_idx_type n, fft_len;

  fft_len = nx+ny-1;

  //
  // Copy and zero-pad.
  //

  for (n=0; n < nx; n++) {
    a[n] = xr[n];
  }

  for (n=nx; n < fft_len; n++)
    a[n] = 0.0; // Zero-pad.

  for (n=0; n < ny; n++)
    b[n] = yr[n];
  for (n=ny; n < fft_len; n++)
    b[n] = 0.0; // Zero-pad.

  // Fourier transform xr.
  FFTVec a_v(fft_len, a);
  FFTCVec af_v(fft_len, af);
  fft.fft(a_v, af_v);

  // Fourier transform yr.
  FFTVec b_v(fft_len, b);
  FFTCVec bf_v(fft_len, bf);
  fft.fft(b_v, bf_v);

  // Do the filtering.
  for (n = 0; n < fft_len; n++) {
    cf[n] = (af[n] * bf[n]);
  }

  //
  // Compute the inverse DFT of the filtered data.
  //

  FFTCVec cf_v(fft_len, cf);
  FFTVec c_v(fft_len, c);
  fft.ifft(cf_v, c_v);

  switch (conv_mode) {

  case ConvMode::sum:
    {
      // in-place '+=' operation.
      for (n = 0; n < fft_len; n++) {
        zr[n] += c[n];
      }
    }
    break;

  case ConvMode::neg:
    {
      // in-place '-=' operation.
      for (n = 0; n < fft_len; n++) {
        zr[n] -= c[n];
      }
    }
    break;

  case ConvMode::equ:
  default:
    {
      // '=' operation.
      memcpy(zr,c,fft_len*sizeof(double));
      break;
    }
  }
}
