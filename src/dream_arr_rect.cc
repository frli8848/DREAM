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


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "dream_arr_rect.h"
#include "arr_functions.h"
#include "dream_error.h"

//
// Function prototypes
//

int rect_ab(double xo, double yo, double zo,
            double x, double y, double z,
            double a, double b,
            double dx, double dy, double dt,
            dream_idx_type nt,
            double delay, double foc_delay, double steer_delay,
            double v, double cp, double weight,
            double *h, int err_level);

int rect_ab(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
            double xo, double yo, double zo,
            double x, double y, double z,
            double a, double b,
            double dx, double dy, double dt,
            dream_idx_type nt,
            double delay, double foc_delay, double steer_delay,
            double v, double cp, double weight,
            double *h, int err_level);

/***
 *
 * dream_arr_rect
 *
 ***/

int dream_arr_rect(double xo, double yo, double zo,
                   double a, double b,
                   double dx, double dy, double dt,
                   dream_idx_type nt, double delay, double v, double cp,
                   dream_idx_type num_elements, double *gx, double *gy, double *gz,
                   int foc_type, double *focal,
                   int steer_type, double theta, double phi,
                   double *apod, bool do_apod, int apod_type, double param,
                   double *h, int err_level)
{
  int err = NONE, out_err = NONE;

  for (dream_idx_type i=0; i<nt; i++) {
    h[i] = 0.0;
  }

  double foc_delay=0.0, steer_delay=0.0, weight=1.0;
  double ramax, xamax, yamax;
  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, num_elements);

  for (dream_idx_type n=0; n<num_elements; n++) {

    if (foc_type != FOCUS_UD) {
      focusing(foc_type, focal[0], gx[n], gy[n], xamax, yamax, ramax, cp, &foc_delay);
    } else {
      focusing(foc_type, focal[n], gx[n], gy[n], xamax, yamax, ramax, cp, &foc_delay);
    }

    beamsteering(steer_type, theta, phi, gx[n], gy[n], xamax, yamax, ramax, cp, &steer_delay);

    if (do_apod) {
      apodization(apod_type, n, apod, &weight, gx[n], gy[n], ramax, param);
    }

    // Compute the response for the n:th elemen and add it to the impulse response vector ha.
    err = rect_ab(xo, yo, zo, gx[n], gy[n], gz[n], a, b, dx, dy, dt, nt,
                  delay, foc_delay, steer_delay, v, cp, weight, h, err_level);

    if (err != NONE) {
      out_err = err;
    }
  }

  return out_err;
}

int dream_arr_rect(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                   double xo, double yo, double zo,
                   double a, double b,
                   double dx, double dy, double dt,
                   dream_idx_type nt, double delay, double v, double cp,
                   dream_idx_type num_elements, double *gx, double *gy, double *gz,
                   int foc_type, double *focal,
                   int steer_type, double theta, double phi,
                   double *apod, bool do_apod, int apod_type, double param,
                   double *h, int err_level)
{
  int err = NONE, out_err = NONE;

  for (dream_idx_type i=0; i<nt; i++) {
    h[i] = 0.0;
  }


  double foc_delay=0.0, steer_delay=0.0, weight=1.0;
  double ramax, xamax, yamax;
  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, num_elements);

  for (dream_idx_type n=0; n<num_elements; n++) {

    if (foc_type != FOCUS_UD) {
      focusing(foc_type, focal[0], gx[n], gy[n], xamax, yamax, ramax, cp, &foc_delay);
    } else {
      focusing(foc_type, focal[n], gx[n], gy[n], xamax, yamax, ramax, cp, &foc_delay);
    }

    beamsteering(steer_type, theta, phi, gx[n], gy[n], xamax, yamax, ramax, cp, &steer_delay);

    if (do_apod) {
      apodization(apod_type, n, apod, &weight, gx[n], gy[n], ramax, param);
    }

    // Compute the response for the n:th elemen and add it to the impulse response vector ha.
    err = rect_ab(att, xc_vec, x_vec,
                  xo, yo, zo, gx[n], gy[n], gz[n], a, b, dx, dy, dt, nt,
                  delay, foc_delay, steer_delay, v, cp, weight, h, err_level);

    if (err != NONE) {
      out_err = err;
    }
  }

  return out_err;
}

/***
 *
 * rect_ab
 *
 * Computes the impulse respone of one rectangular element for use within an array.
 *
 * NB. We add (superimpose) the response to impulse response vector h!
 *
 ***/

int rect_ab(double xo, double yo, double zo,
            double x, double y, double z, // Element center position
            double a, double b,
            double dx, double dy, double dt,
            dream_idx_type nt,
            double delay, double foc_delay, double steer_delay,
            double v, double cp, double weight,
            double *h, int err_level)
{
   int err = NONE;

  double ds = dx * dy;

  double xsmin = x - a/2;
  double xsmax = x + a/2;

  double ysmin = y - b/2;
  double ysmax = y + b/2;

  double zs = 0.0;

  // Loop over all surface elements (xs, ys)

  double ys = ysmin + dy/2.0;
  while (ys <= ysmax) {

    double xs = xsmin + dx/2.0;
    while (xs <= xsmax) {

      double r;
      distance(xo, yo, zo, xs, ys, zs, &r);
      double t = r * 1.0e3/cp; // Propagation delay in micro seconds.
      dream_idx_type it = (dream_idx_type) rint((t - delay + foc_delay + steer_delay)/dt);

      if ((it < nt) && (it >= 0)) {

        double ai = weight * v * ds / (2*M_PI*r);
        ai /= dt;
        ai *= 1.0e3;            // Convert to SI units.
        h[it] += ai;

      } else {

        if  (it >= 0)
          err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
        else
          err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

        if ( (err_level == PARALLEL_STOP) || (err_level == STOP) )
          return err; // Bail out.
      }

      xs += dx;
    }
    ys += dy;
  }

  return err;
}


int rect_ab(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
            double xo, double yo, double zo,
            double x, double y, double z,
            double a, double b,
            double dx, double dy, double dt,
            dream_idx_type nt,
            double delay, double foc_delay, double steer_delay,
            double v, double cp, double weight,
            double *h, int err_level)
{
  int err = NONE;

  double ds = dx * dy;

  double xsmin = x - a/2;
  double xsmax = x + a/2;

  double ysmin = y - b/2;
  double ysmax = y + b/2;

  double zs = 0.0;

  double ys = ysmin + dy/2.0;
  while (ys <= ysmax) {

    double xs = xsmin + dx/2.0;
    while (xs <= xsmax) {

      double r;
      distance(xo, yo, zo, xs, ys, zs, &r);
      double t = r * 1.0e3/cp; // Propagation delay in micro seconds.
      dream_idx_type it = (dream_idx_type) rint((t - delay + foc_delay + steer_delay)/dt);

      if ((it < nt) && (it >= 0)) {

        double ai = weight * v * ds / (2*M_PI*r);
        ai /= dt;
        ai *= 1.0e3;            // Convert to SI units.

        att.att(xc_vec, x_vec, r, it, h, ai);

      } else {

        if  (it >= 0)
          err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
        else
          err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

        if ( (err_level == PARALLEL_STOP) || (err_level == STOP) ) {
          return err; // Bail out.
        }
      }

      xs += dx;
    }
    ys += dy;
  }

  return err;
}
