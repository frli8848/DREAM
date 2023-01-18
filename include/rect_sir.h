/***
*
* Copyright (C) 2004,2006,2007,2008,2009,2023 Fredrik Lingvall
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
#include "dream_error.h"

class RectSir
{
 public:

 RectSir()
   : m_out_err(ErrorLevel::none)
    {;}

  ~RectSir()  = default;

ErrorLevel rect_sir(double *ro, dream_idx_type no,
                    double a, double b,
                    double dt, dream_idx_type nt,
                    DelayType delay_type, double *delay,
                    double v, double cp,
                    double *h);

#ifdef USE_OPENCL
 int cl_rect_sir(const double *ro,
                 int no,
                 double a,
                 double b,
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

 void* smp_rect_sir(void *arg);
 std::thread rect_sir_thread(void *arg) {
   return std::thread(&RectSir::smp_rect_sir, this, arg);
 }

 void rect_sir_serial(double xo,
                      double yo,
                      double zo,
                      double a,
                      double b,
                      double dt,
                      dream_idx_type nt,
                      double delay,
                      double v,
                      double cp,
                      double *h);

  ErrorLevel m_out_err;
};
