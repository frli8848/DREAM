/***
*
* Copyright (C) 2003,2006,2007,2008,2009,2014,2019,2021 Fredrik Lingvall
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

#include <string.h>
#include <math.h>

#include "dreamcirc_f.h"
#include "dream_error.h"

//
//  Function prototypes.
//

double focusing_circ_f(int foc_type, double focal,
                     double xs, double ys,
                     double ramax,
                     double cp);
/***
 *
 *  dreamcirc_f - Focused circular transducer.
 *
 ***/

int dreamcirc_f(double xo, double yo, double zo,
                double R, int foc_type, double focal,
                double dx, double dy, double dt,
                dream_idx_type nt, double delay, double v, double cp,
                double *h,
                int err_level)
{
  dream_idx_type i, it;
  double t, ai, ds, r;
  double xsmin, ysmin, xsmax, ysmax;
  double x, y, foc_delay;
  int err = NONE;

  ds = dx * dy;

  ysmin = -R;
  ysmax =  R;

  for (i = 0; i < nt; i++) {
    h[i] = 0.0 ;
  }

  y = ysmin + dy/2;
  while (y <= ysmax) {

    //xlimit_circ_f(y, R, xs, ys, &xsmin, &xsmax);
    double rs = sqrt(R*R - y*y);
    xsmin = -rs;
    xsmax = rs;

    double ry = yo - y;

    x = xsmin + dx/2;
    while (x <= xsmax) {

      // Compute the distance (length) from an observation point (xo,yo,zo)
      // to a point (x,y) on the transducer surface.
      //distance_circ_f(xo, yo, zo, x, y, &r);
      double rx = xo - x;
      r = sqrt(rx*rx + ry*ry + zo*zo);

      foc_delay = focusing_circ_f(foc_type, focal,
                                  x, y,
                                  r,
                                  cp);

      ai = v * ds / (2*M_PI * r);
      ai /= dt;
      ai *= 1.0e3; // Convert to SI units.

      // Propagation delay in micro seconds.
      t = r * 1.0e3/cp;
      it = (dream_idx_type) rint((t - delay + foc_delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {
        h[it] += ai;
      } else  {

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

int dreamcirc_f(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                double xo, double yo, double zo,
                double R, int foc_type, double focal,
                double dx, double dy, double dt,
                dream_idx_type nt, double delay, double v, double cp,
                double *h,
                int err_level)
{
  dream_idx_type i, it;
  double t, ai;
  double xsmin, ysmin, xsmax, ysmax, ds, r;
  double x, y, foc_delay;
  int err = NONE;

  ds = dx * dy;
  ysmin = -R;
  ysmax =  R;

  for (i = 0; i < nt; i++) {
    h[i] = 0.0 ;
  }

  y = ysmin + dy/2;
  while (y <= ysmax) {

    // Compute the x-axis integration limits.
    //xlimit_circ_f(y, R, xs, ys, &xsmin, &xsmax);
    double rs = sqrt(R*R - y*y);
    xsmin = -rs;
    xsmax = rs;

    double ry = yo - y;

    x = xsmin + dx/2;
    while (x <= xsmax) {

      // Compute the distance (length) from an observation point (xo,yo,zo)
      // to a point (x,y) on the transducer surface.
      double rx = xo - x;
      r = sqrt(rx*rx + ry*ry + zo*zo);

      foc_delay = focusing_circ_f(foc_type, focal,
                                  x, y,
                                  r,
                                  cp);

      ai = v * ds / (2*M_PI * r);
      ai /= dt;
      ai *= 1.0e3; // Convert to SI units.

      // Propagation delay in micro seconds.
      t = r * 1.0e3/cp;
      it = (dream_idx_type) rint((t - delay + foc_delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {
        att.att(xc_vec, x_vec, r, it, h, ai);
      } else  {
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
 *  Computes the focus delay
 *
 ***/

double focusing_circ_f(int foc_type, double focal,
             double xs, double ys,
             double ramax,
             double cp)
{
  double diff, rmax, x_delay, y_delay;
  double foc_delay = 0.0;

  // foc_type: 1 no foc, 2 foc x, 3 foc y, 4 foc xy, 5 foc x+y
  switch (foc_type) {

  case 1:
    break;

  case 2:
    rmax = sqrt(ramax*ramax + focal*focal);
    diff = rmax - sqrt(xs*xs + focal*focal);
    foc_delay = diff * 1.0e3 / cp;
    break;

  case 3:
    rmax = sqrt(ramax*ramax + focal*focal);
    diff = rmax - sqrt(ys*ys + focal*focal);
    foc_delay = diff * 1.0e3 / cp;
    break;

  case 4:
    rmax = sqrt(ramax*ramax + focal*focal);
    diff = rmax - sqrt(xs*xs + ys*ys + focal*focal);
    foc_delay = diff * 1.0e3 / cp;
    break;

  case 5:
    rmax = sqrt(ramax*ramax + focal*focal);
    x_delay = sqrt(xs*xs + focal*focal);
    y_delay = sqrt(ys*ys + focal*focal);
    diff = rmax - (x_delay + y_delay);
    foc_delay = diff * 1.0e3 / cp;

  default:
    break;

  }

  return foc_delay;
}
