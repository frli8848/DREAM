/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2014,2019 Fredrik Lingvall
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

#include "dream_error.h"
#include "dream_arr_rect.h"
#include "att.h"
#include "arr_functions.h"

#if defined(_MSC_VER) || defined(__LCC__)
#include "msvc_rint.h"
#endif

//
// Function prototypes
//


int rect_ab(double xo, double yo, double zo, double xs, double ys, double zs, double a, double b,
            double dx, double dy, double dt, dream_idx_type nt, double delay, double retfoc,
            double retsteer, double v, double cp, double alpha, double weight, double *RESTRICT h, int err_level);


/***
 *
 * Subroutine dream_arr_rect
 *
 ***/

int dream_arr_rect(double xo, double yo, double zo, double a, double b, double dx, double dy, double dt,
                  dream_idx_type nt, double delay, double v, double cp, double alpha, int num_elements,
                  double *RESTRICT gx, double *RESTRICT gy, double *RESTRICT gz, int foc_type, double focal, int ister,
                  double theta, double phi, double *RESTRICT apod, bool do_apod, int apod_type, double param, double *RESTRICT ha,
                  int err_level)
{
  double retsteer;
  dream_idx_type i;
  double ramax, xamax, yamax;
  double xs, ys, zs, retfoc, weight;
  int err = NONE, out_err = NONE;

  for (i=0; i<nt; i++) {
    ha[i] = 0.0;
  }

  retfoc   = 0.0;
  retsteer = 0.0;
  weight   = 1.0;

  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, num_elements);

  for (i=0; i<num_elements; i++) {
    center_pos(&xs, &ys, &zs, i, gx, gy, gz);
    focusing(foc_type, focal, xs, ys, xamax, yamax, ramax, cp, &retfoc);
    beamsteering(ister, theta, phi, xs, ys, xamax, yamax, ramax, cp, &retsteer);
    if (do_apod) {
      apodization(apod_type, i, apod, &weight, xs, ys, ramax, param);
    }
    // Compute the response for the i:th elemen and add it to the impulse response vector ha.
    err = rect_ab(xo, yo, zo, xs, ys, zs, a, b, dx, dy, dt, nt,
                  delay, retfoc, retsteer, v, cp, alpha, weight, ha, err_level);
    if (err != NONE)
      out_err = err;
  }

  return out_err;
}

/***
 *
 * Subroutine dream_arr_rect_ud - user defined focusing.
 *
 ***/

int dream_arr_rect_ud(double xo, double yo, double zo, double a, double b, double dx, double dy, double dt,
                  dream_idx_type nt, double delay, double v, double cp, double alpha, int num_elements,
                  double *RESTRICT gx, double *RESTRICT gy, double *RESTRICT gz, int foc_type, double *RESTRICT focal, int ister,
                  double theta, double phi, double *RESTRICT apod, bool do_apod, int apod_type, double param, double *RESTRICT ha,
                  int err_level)
{
  double retsteer;
  double *RESTRICT h;
  dream_idx_type i;
  double ramax, xamax, yamax;
  double xs, ys, zs, retfoc, weight;
  int err = NONE, out_err = NONE;

  for (i=0; i<nt; i++) {
    ha[i] = (double) 0.0;
  }

  retfoc   = (double) 0.0;
  retsteer = (double) 0.0;
  weight   = (double) 1.0;

  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, num_elements);

  for (i=0; i<num_elements; i++) {
    center_pos(&xs, &ys, &zs, i, gx, gy, gz);
    focusing(foc_type, focal[i], xs, ys, xamax, yamax, ramax, cp, &retfoc);  // Note foc_type must be 6 here!
    beamsteering(ister, theta, phi, xs, ys, xamax, yamax, ramax, cp, &retsteer);
    if (do_apod) {
      apodization(apod_type, i, apod, &weight, xs, ys, ramax, param);
    }
    // Compute the response for the i:th element and add it to the impulse response vector ha.
    err = rect_ab(xo, yo, zo, xs, ys, zs, a, b, dx, dy, dt, nt,
                  delay, retfoc, retsteer, v, cp, alpha, weight, ha, err_level);

    if (err != NONE)
      out_err = err;
  }

  return out_err;
}


/********************************************************************/

/***
 *
 * rect_ab
 *
 * Computes the impulse respone of one rectangular element for use within an array.
 *
 * NB. We add (super impose) the response to impulse response vector h!
 *
 ***/

int rect_ab(double xo, double yo, double zo, double xs, double ys, double zs, double a, double b,
            double dx, double dy, double dt, dream_idx_type nt, double delay, double retfoc,
            double retsteer, double v, double cp, double alpha, double weight, double *RESTRICT h, int err_level)
{
  dream_idx_type i;
  double t;
  double xsmin, ysmin, xsmax, ysmax, ds, pi, r;
  dream_idx_type it;
  double x, y;
  int err = NONE;

  double decal = retfoc + retsteer;
  pi = 4.0 * atan(1.0);
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
      t = r * 1000.0/cp; // Propagation delay in micro seconds.
      it = (dream_idx_type) rint((t - delay + decal)/dt);
      if ((it < nt) && (it >= 0)) {

        double ai = weight * v * ds / (2*pi*r);
        ai /= dt;
        ai *= 1000.0;      // Convert to SI units.

        if (alpha == 0.0 ) {
          h[it] += ai; // No attenuation.
        } else {
          att(alpha, r, it, dt, cp, h, nt, ai); // Attenuation compensation.
        }
      } else {

        if  (it >= 0)
          err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
        else
          err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

        if ( (err_level == PARALLEL_STOP) || (err_level == STOP) )
          return err; // Bail out.
      }

      x += dx;
    } // while
    y += dy;
  } // while

  return err;
} // rect_ab
