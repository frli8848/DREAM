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

class Rect
{
 public:

 Rect()
   : m_out_err(ErrorLevel::none)
    {;}

  ~Rect()  = default;

ErrorLevel dreamrect(double alpha,
                     double *ro, dream_idx_type no,
                     double a, double b,
                     double dx, double dy, double dt, dream_idx_type nt,
                     DelayType delay_type, double *delay,
                     double v, double cp,
                     double *h, ErrorLevel err_level);

#ifdef USE_OPENCL
  int cl_dreamrect(const double *Ro,
                   int No,
                   double a,
                   double b,
                   double dx,
                   double dy,
                   double dt,
                   int nt,
                   double delay,
                   double v,
                   double cp,
                   double *h);
#endif


 static void abort(int signum);
 bool is_running();

 private:

 void* smp_dream_rect(void *arg);
 std::thread rect_thread(void *arg) {
   return std::thread(&Rect::smp_dream_rect, this, arg);
 }

 ErrorLevel m_out_err;
};

// FIXME. We keep the serial code out of the class until
// we have benchmarked how much overhead it is to have them as
// member functions in the array class(es) that also uses these
// functions.

ErrorLevel dreamrect_serial(double xo, double yo, double zo,
                            double a, double b,
                            double dx, double dy, double dt,
                            dream_idx_type nt,
                            double delay,
                            double v, double cp,
                            double *h,
                            ErrorLevel err_level,
                            double weight=1.0);

ErrorLevel dreamrect_serial(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                            double xo, double yo, double zo,
                            double a, double b,
                            double dx, double dy, double dt,
                            dream_idx_type nt,
                            double delay,
                            double v, double cp,
                            double *h,
                            ErrorLevel err_level,
                            double weight=1.0);
