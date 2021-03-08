/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2014,2019,2021 Fredrik Lingvall
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

#include <cmath>

#include "dream_arr_rect.h"
#include "dreamrect.h"
#include "arr_functions.h"

ErrorLevel dream_arr_rect(double xo, double yo, double zo,
                          double a, double b,
                          double dx, double dy, double dt,
                          dream_idx_type nt,
                          double delay, double v, double cp,
                          dream_idx_type num_elements, double *gx, double *gy, double *gz,
                          FocusMet foc_met, double *focal,
                          SteerMet steer_met, double theta, double phi,
                          double *apod, bool do_apod, ApodMet apod_met, double apod_par,
                          double *h, ErrorLevel err_level)
{
  ErrorLevel err = ErrorLevel::none, out_err = ErrorLevel::none;

  double r_max, x_max, y_max;
  max_dim_arr(&x_max, &y_max, &r_max, gx, gy, gz, num_elements);

  for (dream_idx_type n=0; n<num_elements; n++) {

    double foc_delay = 0.0;
    if (foc_met != FocusMet::ud) {
      focusing(foc_met, focal[0], gx[n], gy[n], x_max, y_max, r_max, cp, &foc_delay);
    } else {
      focusing(foc_met, focal[n], gx[n], gy[n], x_max, y_max, r_max, cp, &foc_delay);
    }

    double steer_delay = 0.0;
    beamsteering(steer_met, theta, phi, gx[n], gy[n], x_max, y_max, r_max, cp, &steer_delay);

    double weight = 1.0;
    if (do_apod) {
      apodization(apod_met, n, apod, &weight, gx[n], gy[n], r_max, apod_par);
    }

    // Compute the response for the n:th elemen and add it to the impulse response vector h.
    err = dreamrect(xo - gx[n], yo - gy[n], zo - gz[n],
                    a, b, dx, dy, dt, nt,
                    delay - foc_delay - steer_delay,
                    v, cp,
                    h,
                    err_level,
                    weight);

    if (err != ErrorLevel::none) {
      out_err = err;
    }
  }

  return out_err;
}

ErrorLevel dream_arr_rect(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                          double xo, double yo, double zo,
                          double a, double b,
                          double dx, double dy, double dt,
                          dream_idx_type nt, double delay, double v, double cp,
                          dream_idx_type num_elements, double *gx, double *gy, double *gz,
                          FocusMet foc_met, double *focal,
                          SteerMet steer_met, double theta, double phi,
                          double *apod, bool do_apod, ApodMet apod_met, double apod_par,
                          double *h, ErrorLevel err_level)
{
  ErrorLevel err = ErrorLevel::none, out_err = ErrorLevel::none;

  for (dream_idx_type i=0; i<nt; i++) {
    h[i] = 0.0;
  }


  double r_max, x_max, y_max;
  max_dim_arr(&x_max, &y_max, &r_max, gx, gy, gz, num_elements);

  for (dream_idx_type n=0; n<num_elements; n++) {

    double foc_delay = 0.0;
    if (foc_met != FocusMet::ud) {
      focusing(foc_met, focal[0], gx[n], gy[n], x_max, y_max, r_max, cp, &foc_delay);
    } else {
      focusing(foc_met, focal[n], gx[n], gy[n], x_max, y_max, r_max, cp, &foc_delay);
    }

    double steer_delay = 0.0;
    beamsteering(steer_met, theta, phi, gx[n], gy[n], x_max, y_max, r_max, cp, &steer_delay);

    double weight = 1.0;
    if (do_apod) {
      apodization(apod_met, n, apod, &weight, gx[n], gy[n], r_max, apod_par);
    }

    // Compute the response for the n:th elemen and add it to the impulse response vector h.
    err = dreamrect(att, xc_vec, x_vec,
                    xo - gx[n], yo - gy[n], zo - gz[n],
                    a, b, dx, dy, dt, nt,
                    delay - foc_delay - steer_delay,
                    v, cp,
                    h,
                    err_level,
                    weight);
    if (err != ErrorLevel::none) {
      out_err = err;
    }
  }

  return out_err;
}
