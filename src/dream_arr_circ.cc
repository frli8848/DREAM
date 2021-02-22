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

#include "dream_arr_circ.h"
#include "arr_functions.h"
#include "dream_error.h"

//
// Function prototypes.
//

int circ_arr(double xo, double yo, double zo,
             double x, double y,
             double R,
             double dx, double dy, double dt, dream_idx_type nt,
             double delay, double foc_delay, double steer_delay,
             double v, double cp,
             double weight,
             double *h, int err_level);

int circ_arr(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
             double xo, double yo, double zo,
             double x, double y,
             double R,
             double dx, double dy, double dt, dream_idx_type nt,
             double delay, double foc_delay, double steer_delay,
             double v, double cp,
             double weight,
             double *h, int err_level);


/***
 *
 * dream_arr_cir - 2D array with circular elements.
 *
 ***/

int dream_arr_circ(double xo, double yo, double zo,
                   double R,
                   double dx, double dy, double dt, dream_idx_type nt,
                   double delay,
                   double v, double cp,
                   int  num_elements, double *gx, double *gy, double *gz,
                   int foc_type, double *focal,
                   int steer_type, double theta, double phi, double *apod, bool do_apod,
                   int apod_type, double param, double *h, int err_level)
{
  dream_idx_type i;
  double ramax, xamax, yamax;
  double foc_delay, steer_delay, weight;
  int err = NONE, out_err = NONE;

  for (i=0; i<nt; i++) {
    h[i] = 0.0;
  }

  foc_delay   = 0.0;
  steer_delay = 0.0;
  weight      = 1.0;

  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, num_elements);

  for (i=0; i<num_elements; i++) {

    if (foc_type != FOCUS_UD) {
      focusing(foc_type, focal[0], gx[i], gy[i], xamax, yamax, ramax, cp, &foc_delay);
    } else {
      focusing(foc_type, focal[i], gx[i], gy[i], xamax, yamax, ramax, cp, &foc_delay);
    }

    beamsteering(steer_type, theta, phi, gx[i], gy[i], xamax, yamax, ramax, cp, &steer_delay);

    if (do_apod) {
      apodization(apod_type, i, apod, &weight, gx[i], gy[i], ramax, param);
    }

    // Compute the response for the i:th element and add it to the impulse response vector h.
    err = circ_arr(xo, yo, zo,
                   gx[i], gy[i],
                   R,
                   dx, dy, dt, nt,
                   delay, foc_delay, steer_delay,
                   v, cp, weight,
                   h, err_level);

    if (err != NONE)
      out_err = err;
  }

  return out_err;
}

int dream_arr_circ(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                   double xo, double yo, double zo,
                   double R,
                   double dx, double dy, double dt, dream_idx_type nt,
                   double delay,
                   double v, double cp,
                   int  num_elements, double *gx, double *gy, double *gz,
                   int foc_type, double *focal,
                   int steer_type, double theta, double phi, double *apod, bool do_apod,
                   int apod_type, double param, double *h, int err_level)
{
  dream_idx_type i;
  double ramax, xamax, yamax;
  double foc_delay, steer_delay, weight;
  int err = NONE, out_err = NONE;

  for (i=0; i<nt; i++) {
    h[i] = 0.0;
  }

  foc_delay   = 0.0;
  steer_delay = 0.0;
  weight      = 1.0;

  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, num_elements);

  for (i=0; i<num_elements; i++) {

    if (foc_type != FOCUS_UD) {
      focusing(foc_type, focal[0], gx[i], gy[i], xamax, yamax, ramax, cp, &foc_delay);
    } else {
      focusing(foc_type, focal[i], gx[i], gy[i], xamax, yamax, ramax, cp, &foc_delay);
    }

    beamsteering(steer_type, theta, phi, gx[i], gy[i], xamax, yamax, ramax, cp, &steer_delay);

    if (do_apod) {
      apodization(apod_type, i, apod, &weight, gx[i], gy[i], ramax, param);
    }

    // Compute the response for the i:th element and add it to the impulse response vector h.
    err = circ_arr(att, xc_vec, x_vec,
                   xo, yo, zo,
                   gx[i], gy[i],
                   R,
                   dx, dy, dt, nt,
                   delay, foc_delay, steer_delay,
                   v, cp, weight,
                   h, err_level);

    if (err != NONE)
      out_err = err;
  }

  return out_err;
}



/***
 *
 * circ_arr
 *
 * NB. We add (super impose) the response to impulse response vector h!
 *
 ***/

int circ_arr(double xo, double yo, double zo,
             double x, double y, // Element center position
             double R,
             double dx, double dy, double dt, dream_idx_type nt,
             double delay, double foc_delay, double steer_delay,
             double v, double cp,
             double weight,
             double *h, int err_level)
{
  double t;
  double x_min, ysmin, x_max, ysmax, ai, ds;
  dream_idx_type  it;
  double xs, ys, zs;
  int err = NONE;

  ds = dx * dy;
  zs = (double) 0.0;
  ysmin = -R + ys;
  ysmax =  R + ys;

  // Loop over all surface elements (xs, ys)

  ys = ysmin + dy/2.0;
  while (ys <= ysmax) {

    // Compute the x-axis integration limits.
    double rs = sqrt(R*R - (y-ys)*(y-ys));
    x_min = -rs + x;
    x_max =  rs + x;

    xs = x_min + dx/2.0;
    while (xs <= x_max) {

      double r;
      distance(xo, yo, zo, xs, ys, zs, &r);
      ai = weight * v * ds / (2*M_PI * r);
      ai /= dt;
      ai *= 1.0e3; // Convert to SI units.

      // Propagation delay in micro seconds.
      t = r * 1.0e3/cp;
      it = (dream_idx_type) rint((t - delay + foc_delay + steer_delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {

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

int circ_arr(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
             double xo, double yo, double zo,
             double x, double y, // Element center position
             double R,
             double dx, double dy, double dt, dream_idx_type nt,
             double delay, double foc_delay, double steer_delay,
             double v, double cp,
             double weight,
             double *h, int err_level)
{
  double t;
  double x_min, ysmin, x_max, ysmax, ai;
  dream_idx_type  it;
  double xs, ys, zs;
  int err = NONE;

  double ds = dx * dy;
  zs = (double) 0.0;
  ysmin = -R + ys;
  ysmax =  R + ys;

  // Loop over all surface elements (xs, ys)

  ys = ysmin + dy/2.0;
  while (ys <= ysmax) {

    // Compute the x-axis integration limits.
    double rs = sqrt(R*R - (y-ys)*(y-ys));
    x_min = -rs + x;
    x_max =  rs + x;

    xs = x_min + dx/2.0;
    while (xs <= x_max) {

      double r;
      distance(xo, yo, zo, xs, ys, zs, &r);
      ai = weight * v * ds / (2*M_PI * r);
      ai /= dt;
      ai *= 1.0e3; // Convert to SI units.

      // Propagation delay in micro seconds.
      t = r * 1.0e3/cp;
      it = (dream_idx_type) rint((t - delay + foc_delay + steer_delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {

        att.att(xc_vec, x_vec, r, it, h, ai);

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
