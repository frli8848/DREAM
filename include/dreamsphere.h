/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2021,2023 Fredrik Lingvall
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

#include <thread>

#include "dream.h"
#include "attenuation.h"
#include "dream_error.h"

class Sphere
{
 public:

 Sphere()
   : m_out_err(ErrorLevel::none)
    {;}

  ~Sphere()  = default;

  ErrorLevel dreamsphere(double alpha,
                         double *Ro, dream_idx_type No,
                         double R, double Rcurv,
                         double dx, double dy, double dt, dream_idx_type nt,
                         DelayType delay_type, double *delay,
                         double v, double cp,
                         double *h, ErrorLevel err_level);

  static void abort(int signum);
  bool is_running();

 private:

  void* smp_dream_sphere(void *arg);
  std::thread sphere_thread(void *arg) {
    return std::thread(&Sphere::smp_dream_sphere, this, arg);
  }

  ErrorLevel dreamsphere_serial(double xo, double yo, double zo,
                                double R, double Rcurv,
                                double dx, double dy, double dt,
                                dream_idx_type nt, double delay, double v, double cp,
                                double  *h, ErrorLevel err_level);

  ErrorLevel dreamsphere_serial(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                                double xo, double yo, double zo,
                                double R, double Rcurv,
                                double dx, double dy, double dt,
                                dream_idx_type nt, double delay, double v, double cp,
                                double  *h, ErrorLevel err_level);

  double sphere_f(double xs, double ys,
                  double Rcurv,
                  double xo, double yo, double zo);

  double sphere_d(double xs, double ys,
                  double Rcurv,
                  double xo, double yo, double zo);

  ErrorLevel m_out_err;
};
