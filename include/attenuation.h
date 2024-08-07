/***
*
* Copyright (C) 2002,2003,2004,2006,2007,2008,2009,2012,2014,2015,2021,2024 Fredrik Lingvall
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

#pragma once

#include <memory>

#include "dream.h"
#include "fft.h"

class Attenuation
{
 public:

  Attenuation(dream_idx_type len, double dt, double alpha) {
    m_fft = std::make_unique<FFT>(len);
    m_len = len;
    m_dt = dt;
    m_alpha = alpha;
  };
  ~Attenuation() = default;

  void att(FFTCVec &xc_vec, FFTVec &x_vec, double r, dream_idx_type it, double *h, double ai, dream_idx_type element=0);
  void print_parameters() {
    std::cout << "len: " << m_len << std::endl;
    std::cout << "dt: " << m_dt << std::endl;
    std::cout << "alpha: " << m_alpha << std::endl;
  }
 private:

  std::unique_ptr<FFT> m_fft;
  dream_idx_type m_len;
  double m_dt;
  double m_alpha;
};
