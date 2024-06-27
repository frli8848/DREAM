/***
*
* Copyright (C) 2002,2003,2005,2006,2007,2008,2009,2019,2021,2023,2024 Fredrik Lingvall
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

class ArrAnnu
{
public:

  ArrAnnu() = default;
  ~ArrAnnu() = default;

  SIRError dream_arr_annu(double alpha,
                          double *Ro, dream_idx_type No,
                          double dx, double dy, double dt, dream_idx_type nt,
                          DelayType delay_type, double *delay,
                          double v, double cp,
                          dream_idx_type num_radii, double *Gr,
                          FocusMet foc_met, double *focal,
                          double *apod, bool do_apod, ApodMet apod_met, double apod_par,
                          double *h, ErrorLevel err_level);

  static void abort(int signum);
  bool is_running();

 private:

  void* smp_dream_arr_annu(void *arg);
  std::thread arr_annu_thread(void *arg) {
    return std::thread(&ArrAnnu::smp_dream_arr_annu, this, arg);
  }

  void annular_disc_radii(double *ring_r, double *Gr, dream_idx_type num_radii);
  double focusing_annular(FocusMet foc_met, double focal, double ring_r, double ring_r_max, double cp);
  double apodization_annular(ApodMet apod_met, dream_idx_type n, double *apod, double ring_r,
                             double ring_r_max, double apod_par);
  void resp_annular(double *h_disc, double *h_ring, dream_idx_type nt, dream_idx_type n);
  void superpos_annular(double *h_ring, double *h, dream_idx_type nt,
                        double weight, double foc_delay, dream_idx_type n, double dt);

  SIRError dream_arr_annu_serial(double xo, double yo, double zo,
                                 double dx, double dy, double dt, dream_idx_type nt,
                                 double delay,
                                 double v, double cp,
                                 dream_idx_type num_radii, double *Gr,
                                 FocusMet foc_met, double *focal,
                                 double *apod, bool do_apod, ApodMet apod_met, double param,
                                 double *h, ErrorLevel err_level);

  SIRError dream_arr_annu_serial(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                                 double xo, double yo, double zo,
                                 double dx, double dy, double dt, dream_idx_type nt,
                                 double delay,
                                 double v, double cp,
                                 dream_idx_type num_radii, double *Gr,
                                 FocusMet foc_met, double *focal,
                                 double *apod, bool do_apod, ApodMet apod_met, double param,
                                 double *h, ErrorLevel err_level);
};
