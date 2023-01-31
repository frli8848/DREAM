/***
*
* Copyright (C) 2002,2003,2004,2006,2007,2008,2009,2012,2014,2021 Fredrik Lingvall
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

#include "dream.h"

#if defined DREAM_OCTAVE || defined HAVE_FFTW
#include <fftw3.h>

class FFTVec
{
 public:

  FFTVec(dream_idx_type len) :
    m_is_allocated(false)
  {
    m_len = len;
    m_v = (double*) fftw_malloc(sizeof(double)*m_len);
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
      fftw_free(m_v);
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
    m_vc = reinterpret_cast<std::complex<double>*>(fftw_malloc(sizeof(fftw_complex)*len));
    m_is_allocated = true;
  };

  FFTCVec(dream_idx_type len, std::complex<double> *vc) :
    m_is_allocated(false)
  {
    m_len = len;
    m_vc = vc;
  }

  ~FFTCVec() {
    if (m_is_allocated) {
      fftw_free(m_vc);
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

  FFT(dream_idx_type fft_len, std::mutex *fft_mutex = nullptr, unsigned plan_method = 4) {

    m_fft_len = fft_len;
    m_fft_mutex = fft_mutex;

    // FFTW's plan functions is not thread safe so (optionally) protect them with a mutex

    FFTVec x(fft_len);
    FFTCVec xc(fft_len);

    unsigned int pm = FFTW_ESTIMATE;

    switch (plan_method) {

    case 1:
      pm = FFTW_EXHAUSTIVE; // Very slow
      break;

    case 2:
      pm = FFTW_PATIENT; // Slow
      break;

    case 3:
      pm = FFTW_MEASURE; // Too slow on long FFTs.
      break;

    default:
    case 4:
      pm = FFTW_ESTIMATE;
      break;
    }

    if (fft_mutex) { //  Guard the FFTW planers with a mutex.

      const std::lock_guard<std::mutex> lock(m_fft_mutex[0]);

      m_p_forward = fftw_plan_dft_r2c_1d( (int) fft_len, x.get(), reinterpret_cast<fftw_complex*>(xc.get()), pm);
      m_p_backward = fftw_plan_dft_c2r_1d( (int) fft_len, reinterpret_cast<fftw_complex*>(xc.get()), x.get(), pm);

    } else { // Don't use a mutex

      m_p_forward = fftw_plan_dft_r2c_1d( (int) fft_len, x.get(), reinterpret_cast<fftw_complex*>(xc.get()), pm);
      m_p_backward = fftw_plan_dft_c2r_1d( (int)fft_len, reinterpret_cast<fftw_complex*>(xc.get()), x.get(), pm);
    }
  };

  ~FFT() = default;

  void fft(FFTVec &x, FFTCVec &yc) {
    fftw_execute_dft_r2c(m_p_forward, x.get(), reinterpret_cast<fftw_complex*>(yc.get()));
  }

  void ifft(FFTCVec &xc, FFTVec &y) {
    fftw_execute_dft_c2r(m_p_backward, reinterpret_cast<fftw_complex*>(xc.get()), y.get());

    // FFTW's FFT's are not normalized so normalize it:
    for (dream_idx_type k=0; k<m_fft_len; k++) {
      y.get()[k] /= double(m_fft_len);
    }
  }

  void cc_ifft(const FFTCVec &xc, FFTCVec &yc);

  bool import_wisdom(char *str) {return fftw_import_wisdom_from_string(str);}
  bool import_wisdom(std::string str) {return fftw_import_wisdom_from_string(str.c_str());}
  void forget_wisdom() {fftw_forget_wisdom();}
  bool is_wisdom(char *str) {return std::strcmp("fftw_wisdom",str);}
  bool is_wisdom(std::string str) {return std::strcmp("fftw_wisdom",str.c_str());}
  std::string get_wisdom() {
    char *str = fftw_export_wisdom_to_string();
    std::string w_str(str);
    return w_str;
  }

 private:

  dream_idx_type m_fft_len;
  std::mutex *m_fft_mutex;
  fftw_plan m_p_forward;
  fftw_plan m_p_backward;
};

#endif

#if defined(DREAM_MATLAB) && not defined(HAVE_FFTW)
#include "fft_matlab_no_fftw.h"
#endif
