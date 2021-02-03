/***
*
* Copyright (C) 2002,2003,2004,2006,2007,2008,2009,2012,2014,2015,2021 Fredrik Lingvall
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
* Free
* Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
* 02110-1301, USA.
*
***/

#include <memory>

#include "dream.h"
#include "fft.h"

/***
 *
 * Header file for attenuation functions (att.cc)
 *
 ***/

#if defined DREAM_OCTAVE || defined HAVE_FFTW
void att_init(dream_idx_type nt, dream_idx_type n_threads);
void att_close();
#endif

void att(double alpha, double r, dream_idx_type it, double dt, double cp, double *h, dream_idx_type nt, double ai);
void att_annu(double alpha, double r, dream_idx_type it, double  dt, double cp, double *h, dream_idx_type nt, double ai, int ns, int num_elements);

class Attenuation
{
 public:

  Attenuation(dream_idx_type len, double dt, double cp, double alpha, size_t n_threads=1) {
    m_fft = std::make_unique<FFT>(len);
    m_len = len;
    m_dt = dt;
    m_cp = cp;
    m_alpha = alpha;
  };
  ~Attenuation() = default;

  void att(FFTCVec &xc_vec, FFTVec &x_vec, double r, dream_idx_type it, double *h, double ai);
  void att_annu(FFTCVec &xc_vec, FFTVec &x_vec, double r, dream_idx_type it, double *h, double ai, int num_elements);

 private:

  std::unique_ptr<FFT> m_fft;
  dream_idx_type m_len;
  double m_dt;
  double m_cp;
  double m_alpha;
};
