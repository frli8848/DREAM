/***
*
* Copyright (C) 2002,2003,2005,2006,2007,2008,2009,2014,2019,2021 Fredrik Lingvall
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
#include "dreamcirc.h"
#include "att.h"
#include "dream_error.h"

/***
 *
 *  dreamcirc
 *
 * Computes the spatial impulse response of a circular aperture.
 *
 ***/

int dreamcirc(double xo, double yo, double zo,
               double R,
               double dx, double dy, double dt,
               dream_idx_type nt,
               double delay,
               double v, double cp,
               double *h,
               int err_level)
{
  dream_idx_type i, it;
  double t, ai;
  double xsmin, ysmin, xsmax, ysmax, r;
  double rx, ry;
  int    err = NONE;

  double ds = dx * dy;
  ysmin = -R;
  ysmax =  R;

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  double ys = ysmin + dy/2;
  while (ys <= ysmax) {

    double rs = sqrt(R*R - ys*ys);
    xsmin = -rs;
    xsmax = rs;

    ry = yo - ys;

    double xs = xsmin + dx / 2.0;
    while (xs <= xsmax) {

      rx = xo - xs;
      r = sqrt(rx*rx + ry*ry + zo*zo);

      ai = v * ds / (2*M_PI * r);
      ai /= dt;
      ai *= 1.0e3; // Convert to SI units.

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

        if ( (err_level == PARALLEL_STOP) || (err_level == STOP) )
          return err; // Bail out.
      }

      xs += dx;
    }

    ys += dy;
  }

  return err;
}


int dreamcirc(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
              double xo, double yo, double zo,
               double R,
               double dx, double dy, double dt,
               dream_idx_type nt,
               double delay,
               double v, double cp,
               double *h,
               int err_level)
{
  dream_idx_type i, it;
  double t, ai;
  double xsmin, ysmin, xsmax, ysmax, r;
  double rx, ry;
  int    err = NONE;

  double ds = dx * dy;
  ysmin = -R;
  ysmax =  R;

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  double ys = ysmin + dy/2;
  while (ys <= ysmax) {

    double rs = sqrt(R*R - ys*ys);
    xsmin = -rs;
    xsmax = rs;

    ry = yo - ys;

    double xs = xsmin + dx / 2.0;
    while (xs <= xsmax) {

      rx = xo - xs;
      r = sqrt(rx*rx + ry*ry + zo*zo);

      ai = v * ds / (2*M_PI * r);
      ai /= dt;
      ai *= 1000; // Convert to SI units.

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

        if ( (err_level == PARALLEL_STOP) || (err_level == STOP) )
          return err; // Bail out.
      }

      xs += dx;
    }

    ys += dy;
  }

  return err;
}
