/***
*
* Copyright (C) 2002,2003,2005,2006,2007,2008,2009,2021,2023,2024 Fredrik Lingvall
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

#include <thread>

#include "dream.h"
#include "attenuation.h"
#include "dream_error.h"

class Cylind
{
 public:

 Cylind() = default;
  ~Cylind() = default;

  SIRError dreamcylind(double alpha,
                       double *Ro, dream_idx_type No,
                       double a, double b, double Rcurv,
                       double dx, double dy, double dt, dream_idx_type nt,
                       DelayType delay_type, double *delay,
                       double v, double cp,
                       double *h, ErrorLevel err_level);

  static void abort(int signum);
  bool is_running();

  // public since dream_arr_cylind uses these too.
  SIRError dreamcylind_serial(double xo, double yo, double zo,
                              double a, double b, double Rcurv,
                              double dx, double dy, double dt,
                              dream_idx_type nt, double delay, double v, double cp,
                              double *h,
                              ErrorLevel err_level,
                              double weight=1.0);

  SIRError dreamcylind_serial(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                              double xo, double yo, double zo,
                              double a, double b, double Rcurv,
                              double dx, double dy, double dt,
                              dream_idx_type nt, double delay, double v, double cp,
                              double *h,
                              ErrorLevel err_level,
                              double weight=1.0);
private:

  void* smp_dream_cylind(void *arg);
  std::thread cylind_thread(void *arg) {
    return std::thread(&Cylind::smp_dream_cylind, this, arg);
  }

  double cylind_f(double xs, double ys,
                  double Rcurv, double z_Rcurv,
                  double xo, double yo, double zo,
                  double &du);

  double cylind_d(double xs, double ys,
                  double Rcurv, double z_Rcurv,
                  double xo, double yo, double zo,
                  double &du);
};
