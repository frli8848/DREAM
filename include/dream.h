/***
*
* Copyright (C) 2007-2009,2012,2014-2016,2019,2021,2023 Fredrik Lingvall
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


#ifndef __DREAM__
#define __DREAM__

#include <memory>

enum class ErrorLevel {
  none,
  stop,
  warn,
  ignore,
  parallel_stop
};

enum class DelayType {
  single,
  multiple
};

enum class FocusMet {
  none,
  x,
  y,
  xy,
  x_y,
  ud
};

enum class SteerMet {
  none,
  x,
  y,
  xy
};

enum class ApodMet {
  triangle,
  gauss,
  raised_cosine,
  simply_supported,
  clamped,
  ud
};

enum class ConvMode {
  equ,
  sum,
  neg
};

enum class DASType {
  saft,
  tfm,
  rca
};

#ifdef DREAM_OCTAVE

// Use the index type that Octave was build with.
#include <octave-config.h>
typedef octave_idx_type dream_idx_type;

// Macros
#define mxGetM(N)   args(N).matrix_value().rows()
#define mxGetN(N)   args(N).matrix_value().cols()
#define mxIsChar(N) args(N).is_string()

#endif

#ifdef DREAM_MATLAB
#include "mex.h"
typedef mwSize dream_idx_type;
#endif

#ifdef DREAM_PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
typedef pybind11::ssize_t dream_idx_type;
#endif

// FIXME: Check what Julia uses as matrix/vector index type!
#ifdef DREAM_JULIA
typedef size_t dream_idx_type;
#endif

#if !defined(DREAM_OCTAVE) && !defined(DREAM_MATLAB) && !defined(DREAM_PYTHON) && !defined(DREAM_JULIA)
typedef size_t dream_idx_type;
#define HAVE_FFTW
#endif



/***
 *
 * Class for spatial impulse response (SIR) data
 *
 ***/

class SIRData
{
 public:

  SIRData(dream_idx_type len, dream_idx_type no) {
    m_size = len;
    m_len = len;
    m_no = no;
    m_H = std::make_unique<double[]>(len*no);
    m_h = m_H.get();
  }

  SIRData(double *h, dream_idx_type len) {
    m_h = h;
    m_size = len;
    m_len = len;
    m_no = 1;
  }

  SIRData(double *h, dream_idx_type len, dream_idx_type no) {
    m_h = h;
    m_size = len*no;
    m_len = len;
    m_no = no;
  }

  ~SIRData() = default;

  void clear() {
    for (dream_idx_type it=0; it<m_size; it++) {
      m_h[it] = 0.0;
    };
  };

  double* get() {
    return m_h;
  };

  dream_idx_type get_len() {
    return m_len;
  };

  dream_idx_type get_no() {
    return m_no;
  };

 private:

  std::unique_ptr<double[]> m_H;
  double *m_h;
  dream_idx_type m_size;
  dream_idx_type m_len;
  dream_idx_type m_no;

};

#endif // __DREAM__
