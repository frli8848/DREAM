/***
*
* Copyright (C) 2003,2006,2007,2008,2009,2021,2023 Fredrik Lingvall
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

class Rect_f
{
 public:

 Rect_f()
   : m_out_err(ErrorLevel::none)
    {;}

  ~Rect_f()  = default;

ErrorLevel dreamrect_f(double alpha,
                       double *Ro, dream_idx_type No,
                       double a, double b,
                       FocusMet foc_met, double focal,
                       double dx, double dy, double dt, dream_idx_type nt,
                       DelayType delay_type, double *delay,
                       double v, double cp,
                       double *h, ErrorLevel err_level);

 static void abort(int signum);
 bool is_running();

 private:

 void* smp_dream_rect_f(void *arg);
 std::thread rect_f_thread(void *arg) {
   return std::thread(&Rect_f::smp_dream_rect_f, this, arg);
 }

 double focus_delay_rect(FocusMet foc_met, double focal,
                         double xs, double ys,
                         double x_max, double y_max, double r_max,
                         double cp);

 ErrorLevel dreamrect_f_serial(double xo, double yo, double zo,
                               double a, double b,
                               FocusMet foc_met, double focal,
                               double dx, double dy, double dt,
                               dream_idx_type nt,
                               double delay,
                               double v, double cp,
                               double *h,
                               ErrorLevel err_level);

 ErrorLevel dreamrect_f_serial(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                               double xo, double yo, double zo,
                               double a, double b,
                               FocusMet foc_met, double focal,
                               double dx, double dy, double dt,
                               dream_idx_type nt,
                               double delay,
                               double v, double cp,
                               double *h,
                               ErrorLevel err_level);

 ErrorLevel m_out_err;
};
