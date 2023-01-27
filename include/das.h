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

enum class DASType {
  saft,
  tfm
};

class DAS
{
public:

  DAS()
    : m_out_err(ErrorLevel::none)
  {;}

  ~DAS()  = default;

  ErrorLevel das(double *Y, dream_idx_type a_scan_len,
                 double *ro, dream_idx_type no,
                 double *gt, dream_idx_type num_t_elements,
                 double *gr, dream_idx_type num_r_elements, // SAFT if num_r_elements = 0;
                 double dt,
                 DelayType delay_type, double *delay,
                 double cp,
                 double *Im,
                 ErrorLevel err_level);

  static void abort(int signum);
  bool is_running();

#ifdef USE_OPENCL
  int cl_das_tfm(const double *Y, // Data
                 dream_idx_type a_scan_len,
                 const double *Ro, dream_idx_type No,
                 const double *gt, dream_idx_type num_t_elements,
                 const double *gr, dream_idx_type num_r_elements,
                 double ddt,
                 double delay,
                 double cp,
                 double *Im);
#endif

private:

  void* smp_das(void *arg);
  std::thread das_thread(void *arg) {
    return std::thread(&DAS::smp_das, this, arg);
  }

  ErrorLevel das_saft_serial(double *Y, // Size: nt x num_elements
                             dream_idx_type a_scan_len,
                             double *g, dream_idx_type num_elements,
                             double xo, double yo, double zo,
                             double dt, double delay,
                             double cp, double &im, ErrorLevel err_level);

  ErrorLevel das_tfm_serial(double *Y, // Size: nt x num_t_elements*num_r_elements (=FMC)
                            dream_idx_type a_scan_len,
                            double *gt, dream_idx_type num_t_elements,
                            double *gr, dream_idx_type num_r_elements,
                            double xo, double yo, double zo,
                            double dt, double delay,
                            double cp, double &im, ErrorLevel err_level);

  ErrorLevel m_out_err;
};
