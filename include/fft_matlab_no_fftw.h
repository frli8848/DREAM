/***
*
* Copyright (C) 2022 Fredrik Lingvall
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

#pragma once

#include <complex> // Must be included before fftw3.h
#include <mutex>
#include <cstring>
#include <iostream>

#include "mex.h"
#include "dream.h"

class FFTVec
{
 public:

 FFTVec(dream_idx_type len) :
  m_is_allocated(false)
    {
      m_len = len;
      m_v = new double[m_len];
      m_is_allocated = true;
    };

 FFTVec(dream_idx_type len, double *v) :
  m_is_allocated(false)
    {
      m_len = len;
      m_v = v;
    };

  ~FFTVec() {
    if (m_is_allocated) {
      delete[] m_v;
    }
  };

  double* get() {return m_v;};
  dream_idx_type len() {return m_len;};

 private:

  double *m_v;
  dream_idx_type m_len;
  bool m_is_allocated;
};

class FFTCVec
{
 public:

 FFTCVec(dream_idx_type len) :
  m_is_allocated(false)
    {
      m_len = len;
      m_vc = new std::complex<double>[len];
    };

 FFTCVec(dream_idx_type len, std::complex<double> *vc) :
  m_is_allocated(false)
    {
      m_len = len;
      m_vc = vc;
    }

  ~FFTCVec() {
    if (m_is_allocated) {
      delete[] m_vc;
    }
  };

  std::complex<double>* get() {return m_vc;};
  dream_idx_type len() {return m_len;};

 private:

  std::complex<double>* m_vc;
  dream_idx_type m_len;
  bool m_is_allocated;
};

class FFT
{
 public:

  FFT(dream_idx_type fft_len, std::mutex *fft_mutex = nullptr, int plan_method = 4) {
    m_fft_len = fft_len;
    m_fft_mutex = fft_mutex;
  };

  ~FFT() = default;

  void fft(FFTVec &x, FFTCVec &yc) {
    mxArray *X, *Y;
    dream_idx_type len=x.len();

    // FIXME: Matlab thread safe.
    //
    // Here we try to lock the mexCallMATLAB code
    // with a mutex but it does not seem to be enough
    // to avoid segfaults when using mutilple threads!?
    // The FFT must run in the main Matlab thread only!?

    //if (m_fft_mutex) {
    //  const std::lock_guard<std::mutex> lock(m_fft_mutex[0]);

    X = mxCreateDoubleMatrix(len,1,mxREAL);

    //X = mxCreateDoubleMatrix(5,1,mxREAL);
    mxDouble *xp = mxGetDoubles(X);

    // Copy input vector to a Matlab Matrix.
    std::memcpy(xp, x.get(), len*sizeof(double));

    // Let Matlab do the job!
    mexCallMATLAB(1, &Y, 1, &X, "fft");

    // Copy Matlab Matrix to output  vector

#if MX_HAS_INTERLEAVED_COMPLEX
    mxComplexDouble *yp = mxGetComplexDoubles(Y);
    std::memcpy(yc.get(), reinterpret_cast< std::complex<double>*>(yp), len*sizeof(std::complex<double>));
#else
    double *yr= mxGetPr(Y), *yi= mxGetPr(Y);
    std::complex<double>* yc_p = yc.get();
    for (dream_idx_type n=0; n<fft_len; n++) {
      yc_p[n] = std::complex<double>(yr[n], yi[n]);
    }
#endif
    mxDestroyArray(Y);
    mxDestroyArray(X);
    //} // mutex
  }

  void ifft(FFTCVec &xc, FFTVec &y) {
    mxArray *X, *Y;
    dream_idx_type k, len=xc.len();

    // FIXME: Matlab thread safe (see fft above).
    //if (m_fft_mutex) {
    //  const std::lock_guard<std::mutex> lock(m_fft_mutex[0]);

    X = mxCreateDoubleMatrix(len,1,mxCOMPLEX);

    // Copy input vector to a Matlab Matrix.

#if MX_HAS_INTERLEAVED_COMPLEX
    mxComplexDouble *xp = mxGetComplexDoubles(X);
    std::memcpy(xp, xc.get(), len*sizeof(std::complex<double>));
#else
    double *xr= mxGetPr(X), *xi= mxGetPr(X);
    std::complex<double>* xc_p = xc.get();
    for (dream_idx_type n=0; n<fft_len; n++) {
      xr[n] = real(xc_p[n]);
      xi[n] = imag(xc_p[n]);
    }
#endif

    // Let Matlab do the job!
    mexCallMATLAB(1, &Y, 1, &X, "ifft");

    // Use the real part only.
    double *y_out = y.get();

    // We need to check if the result of the FFT is complex
    // otherwise it will crash when calling mxGetComplexDoubles
    // and mxGetPr cannot be used when Y is complex if we have
    // MX_HAS_INTERLEAVED_COMPLEX defined (i.e., for R2018a and above).

    if (mxIsComplex(Y)) {

#if MX_HAS_INTERLEAVED_COMPLEX
      std::complex<double> *yp = reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(Y));
      for (k=0; k<len; k++) {
        y_out[k] = std::real(yp[k]);
      }
#else
      double *yr= mxGetPr(Y);
      for (k=0; k<len; k++) {
        y_out[k] = yr[k];
      }
#endif
    } else {

      double *yr= mxGetPr(Y);
      for (k=0; k<len; k++) {
        y_out[k] = yr[k];
      }
    }

    mxDestroyArray(Y);
    mxDestroyArray(X);
    //} // mutex
  }

  void cc_ifft(const FFTCVec &xc, FFTCVec &yc);

  bool import_wisdom(char *str) {return false;}
  bool import_wisdom(std::string str) {return false;}
  void forget_wisdom() {;}
  bool is_wisdom(char *str) {return false;}
  bool is_wisdom(std::string str) {return false;}
  std::string get_wisdom() {
    char *str = "";
    std::string w_str(str);
    return w_str;
  }

 private:

  dream_idx_type m_fft_len;
  std::mutex *m_fft_mutex;
};
