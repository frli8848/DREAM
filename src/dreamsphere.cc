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

#include "dreamsphere.h"
#include "dream_error.h"

//
// Function prototypes.
//

double sphere_f(double xs, double ys,
                double Rcurv,
                double xo, double yo, double zo);

double sphere_d(double xs, double ys,
                double Rcurv,
                double xo, double yo, double zo);

/***
 *
 *  dreamsphere
 *
 *  Computes the (spatial) impulse response for a focused (concave) spherical aperture.
 *
 ***/

int dreamsphere(double xo, double yo, double zo,
                double R, double Rcurv,
                double dx, double dy, double dt,
                dream_idx_type nt, double delay, double v, double cp,
                double *h, int err_level)
{
  dream_idx_type i, it;
  double t, ai, ds, r;
  double xsmin, ysmin, xsmax, ysmax;
  double x, y;
  int err = NONE;

  ds = dx * dy;

  ysmin = -R;
  ysmax =  R;

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  y = ysmin + dy / 2.0;
  while (y <= ysmax) {

    // Compute the x-axis integration limits.
    double rs = sqrt(R*R - y*y);
    xsmin = -rs;
    xsmax =  rs;

    x = xsmin + dx / 2.0;
    while (x <= xsmax) {

      if (Rcurv >= 0.0) {        // Concave/focused
        r = sphere_f(x, y,
                     Rcurv,
                     xo, yo, zo);
      } else {                  // Convex/de-focused
        r = sphere_d(x, y,
                     -Rcurv,
                     xo, yo, zo);
      }

      ai = v * ds / (2*M_PI * r);
      ai /= dt;
      ai *= 1.0e3;              // Convert to SI units.

      // Propagation delay in micro seconds.
      t = r * 1.0e3/cp;
      it = (dream_idx_type) rint((t - delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {
        h[it] += ai;
      } else {

        if  (it >= 0)
          err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
        else
          err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

        if ( (err_level == PARALLEL_STOP) || (err_level == STOP) ) {
          return err; // Bail out.
        }
      }
      x += dx;
    }
    y += dy;
  }

  return err;
}

int dreamsphere(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                double xo, double yo, double zo,
                double R, double Rcurv,
                double dx, double dy, double dt,
                dream_idx_type nt, double delay, double v, double cp,
                double *h, int err_level)
{
  dream_idx_type i, it;
  double t, ai, ds, r;
  double xsmin, ysmin, xsmax, ysmax;
  double x, y;
  int err = NONE;

  ds = dx * dy;

  ysmin = -R;
  ysmax =  R;

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  y = ysmin + dy / 2.0;
  while (y <= ysmax) {

    double rs = sqrt(R*R - y*y);
    xsmin = -rs;
    xsmax =  rs;

    x = xsmin + dx / 2.0;
    while (x <= xsmax) {

      if (Rcurv >= 0.0) {        // Concave/focused
        r = sphere_f(x, y,
                     Rcurv,
                     xo, yo, zo);
      } else {                  // Convex/de-focused
        r = sphere_d(x, y,
                     -Rcurv,
                     xo, yo, zo);
      }

      ai = v * ds / (2*M_PI * r);
      ai /= dt;
      ai *= 1.0e3;              // Convert to SI units.

      // Propagation delay in micro seconds.
      t = r * 1.0e3/cp;
      it = (dream_idx_type) rint((t - delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {
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
      x += dx;
    }
    y += dy;
  }

  return err;
}


/***
 *
 * sphere_f - focused (concave) sphere
 *
 ***/

double sphere_f(double xs, double ys,
                double Rcurv,
                double xo, double yo, double zo)
{
  double rx, ry, rz;

  double zs = Rcurv - sqrt(Rcurv*Rcurv - xs*xs - ys*ys);

  rx = xo - xs;
  ry = yo - ys;
  rz = zo - zs;
  double r = sqrt(rx*rx + ry*ry + rz*rz);

  return r;
}

/***
 *
 * sphere_d - de-focused (convex) sphere
 *
 ***/

double sphere_d(double xs, double ys,
                double Rcurv,
                double xo, double yo, double zo)
{
  double rx, ry, rz;

  double zs = Rcurv - sqrt(Rcurv*Rcurv - xs*xs - ys*ys);

  rx = xo - xs;
  ry = yo - ys;
  rz = zo + zs;                 // Change the sign of zs for defocused.
  double r = sqrt(rx*rx + ry*ry + rz*rz);

  return r;
}
