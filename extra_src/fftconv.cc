/***
*
* Copyright (C) 2021,2022 Fredrik Lingvall
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
             const double *xr, dream_idx_type nx, const double *yr, dream_idx_type ny, double *zr,
             double *a,  double *b, double *c,
             std::complex<double> *af, std::complex<double> *bf, std::complex<double> *cf,
             ConvMode conv_mode)
{
  dream_idx_type n, fft_len;

  fft_len = nx+ny-1;

  //
  // Copy and zero-pad.
  //

  for (dream_idx_type n=0; n < nx; n++) {
    a[n] = xr[n];
  }

  for (dream_idx_type n=nx; n < fft_len; n++)
    a[n] = 0.0; // Zero-pad.

  for (dream_idx_type n=0; n < ny; n++) {
    b[n] = yr[n];
  }

  for (dream_idx_type n=ny; n < fft_len; n++) {
    b[n] = 0.0; // Zero-pad.
  }

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

/***
 *
 * Convolution of two vectors for each input (array element) and a summation (accumulation) of the result.
 *
 * H    : Handle to to all L impulse response (transmit) matrices (3D matrix) with size L x (h_len x no)
 *        where no is the number of observation points
 * n    : Observation point index [0, no-1]
 * U    : Input/signal matrix with size u_len x L (each column is the input signal/pulse to the corresponding
 *        array element
 * z    : Output vector with length h_len+u_len-1
 *
 ***/

void add_fftconv(FFT &fft,
                 double **H, dream_idx_type L, dream_idx_type h_len, // 3D impulse response matrix.
                 dream_idx_type n, // observation point index.
                 double *U, dream_idx_type u_len, // Input signal matrix.
                 double *z,                        // OUtput vector.
                 double *a, double *b, double *c,
                 std::complex<double> *af, std::complex<double> *bf, std::complex<double> *cf,
                 ConvMode conv_mode)
{
  dream_idx_type fft_len = h_len+u_len-1;

  // Clear output Fourier coefficeints.
  for (dream_idx_type it=0; it < fft_len; it++) {
      cf[it] = 0.0;
  }

  // Loop aver all L inputs (array elements).
  for (dream_idx_type l=0; l<L; l++) {

    // Copy and zero-pad l:th response vector for the n:th observation point.
    for (dream_idx_type it=0; it<h_len; it++) {
      a[it] = (H[l])[it + n*h_len];
    }

    for (dream_idx_type it=h_len; it<fft_len; it++) {
      a[it] = 0.0; // Zero-pad.
    }

    // Copy l:th column (input signal) and zero-pad.
    for (dream_idx_type it=0; it<u_len; it++) {
      b[it] = U[it + l*u_len];
    }

    for (dream_idx_type it=u_len; it<fft_len; it++) {
      b[it] = 0.0; // Zero-pad.
    }

    // Fourier transform a.
    FFTVec a_v(fft_len, a);
    FFTCVec af_v(fft_len, af);
    fft.fft(a_v, af_v);

    // Fourier transform b.
    FFTVec b_v(fft_len, b);
    FFTCVec bf_v(fft_len, bf);
    fft.fft(b_v, bf_v);

    // Do the filtering and add (accumulate) the results
    for (dream_idx_type i_f=0; i_f < fft_len; i_f++) {
      cf[i_f] += (af[i_f] * bf[i_f])  / double(fft_len);
    }
  }

  //
  // Compute the inverse DFT of the filtered and summed data.
  //

  FFTCVec cf_v(fft_len, cf);
  FFTVec c_v(fft_len, c);
  fft.ifft(cf_v, c_v);

  //fftw_execute_dft_c2r(p_backward,reinterpret_cast<fftw_complex*>(cf),c);

  switch (conv_mode) {

  case ConvMode::sum:
    {
      // in-place '+=' operation.
      for (dream_idx_type it = 0; it < fft_len; it++) {
        z[it] += c[it];
      }
    }
    break;

  case ConvMode::equ:
  default:
    {
      // '=' operation.
      memcpy(z, c, fft_len*sizeof(double));
      break;
    }
  }
 }
