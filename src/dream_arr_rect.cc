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
            double xs, double ys, double zs,
            double a, double b,
            double dx, double dy, double dt,
            dream_idx_type nt,
            double delay, double foc_delay, double steer_delay,
            double v, double cp, double weight,
            double *h, int err_level);

int rect_ab(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
            double xo, double yo, double zo,
            double xs, double ys, double zs,
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
  double steer_delay;
  dream_idx_type i;
  double ramax, xamax, yamax;
  double foc_delay, weight;
  int err = NONE, out_err = NONE;

  for (i=0; i<nt; i++) {
    h[i] = 0.0;
  }

  foc_delay   = 0.0;
  steer_delay = 0.0;
  weight = 1.0;

  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, num_elements);

  for (i=0; i<num_elements; i++) {

    if (foc_type != FOCUS_UD) {
      focusing(foc_type, focal[0], gx[i], gx[i], xamax, yamax, ramax, cp, &foc_delay);
    } else {
      focusing(foc_type, focal[i], gx[i], gx[i], xamax, yamax, ramax, cp, &foc_delay);
    }

    beamsteering(steer_type, theta, phi, gx[i], gy[i], xamax, yamax, ramax, cp, &steer_delay);

    if (do_apod) {
      apodization(apod_type, i, apod, &weight, gx[i], gy[i], ramax, param);
    }

    // Compute the response for the i:th elemen and add it to the impulse response vector ha.
    err = rect_ab(xo, yo, zo, gx[i], gy[i], gz[i], a, b, dx, dy, dt, nt,
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
  double steer_delay;
  dream_idx_type i;
  double ramax, xamax, yamax;
  double foc_delay, weight;
  int err = NONE, out_err = NONE;

  for (i=0; i<nt; i++) {
    h[i] = 0.0;
  }

  foc_delay   = 0.0;
  steer_delay = 0.0;
  weight = 1.0;

  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, num_elements);

  for (i=0; i<num_elements; i++) {

    if (foc_type != FOCUS_UD) {
      focusing(foc_type, focal[0], gx[i], gx[i], xamax, yamax, ramax, cp, &foc_delay);
    } else {
      focusing(foc_type, focal[i], gx[i], gx[i], xamax, yamax, ramax, cp, &foc_delay);
    }

    beamsteering(steer_type, theta, phi, gx[i], gy[i], xamax, yamax, ramax, cp, &steer_delay);

    if (do_apod) {
      apodization(apod_type, i, apod, &weight, gx[i], gy[i], ramax, param);
    }

    // Compute the response for the i:th elemen and add it to the impulse response vector ha.
    err = rect_ab(att, xc_vec, x_vec,
                  xo, yo, zo, gx[i], gy[i], gz[i], a, b, dx, dy, dt, nt,
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
            double xs, double ys, double zs,
            double a, double b,
            double dx, double dy, double dt,
            dream_idx_type nt,
            double delay, double foc_delay, double steer_delay,
            double v, double cp, double weight,
            double *h, int err_level)
{
  double t;
  double xsmin, ysmin, xsmax, ysmax, ds, r;
  dream_idx_type it;
  double x, y;
  int err = NONE;

  ds = dx * dy;

  xsmin = xs - a/2;
  xsmax = xs + a/2;

  ysmin = ys - b/2;
  ysmax = ys + b/2;

  y = ysmin + dy/2.0;
  while (y <= ysmax) {

    x = xsmin + dx/2.0;
    while (x <= xsmax) {

      distance(xo, yo, zo, x, y, zs, &r);
      t = r * 1.0e3/cp; // Propagation delay in micro seconds.
      it = (dream_idx_type) rint((t - delay + foc_delay + steer_delay)/dt);

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

      x += dx;
    }
    y += dy;
  }

  return err;
}


int rect_ab(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
            double xo, double yo, double zo,
            double xs, double ys, double zs,
            double a, double b,
            double dx, double dy, double dt,
            dream_idx_type nt,
            double delay, double foc_delay, double steer_delay,
            double v, double cp, double weight,
            double *h, int err_level)
{
  double t;
  double xsmin, ysmin, xsmax, ysmax, ds, r;
  dream_idx_type it;
  double x, y;
  int err = NONE;

  ds = dx * dy;

  xsmin = xs - a/2;
  xsmax = xs + a/2;

  ysmin = ys - b/2;
  ysmax = ys + b/2;

  y = ysmin + dy/2.0;
  while (y <= ysmax) {

    x = xsmin + dx/2.0;
    while (x <= xsmax) {

      distance(xo, yo, zo, x, y, zs, &r);
      t = r * 1.0e3/cp; // Propagation delay in micro seconds.
      it = (dream_idx_type) rint((t - delay + foc_delay + steer_delay)/dt);

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

        if ( (err_level == PARALLEL_STOP) || (err_level == STOP) )
          return err; // Bail out.
      }

      x += dx;
    }
    y += dy;
  }

  return err;
}
