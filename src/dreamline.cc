/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2012,2014,2021 Fredrik Lingvall
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

#include "dreamline.h"

/***
 *
 *  dreamline - Spatial impulse respone of line (slit)
 *
 ***/

ErrorLevel dreamline(double xo, double yo, double zo, double a,
                     double dx, double dy, double dt, dream_idx_type nt, double delay, double v,
                     double cp, double *h, ErrorLevel err_level)
{
  dream_idx_type i, it;
  double t;
  double ai;
  double r;
  ErrorLevel err=ErrorLevel::none;
  double xsmax = a/2.0;
  double xsmin = -a/2.0;

  // dy = width;
  double ds = dx * dy;

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  double xs = xsmin + dx / 2.0;

  while (xs <= xsmax) {

    r = sqrt((xo-xs)*(xo-xs) + yo*yo + zo*zo);

    ai = v * ds / (2*M_PI * r);
    ai /= dt;
    ai *= 1000.0;     // Convert to SI units.

    // Propagation delay in micro seconds.
    t = r * 1000.0/cp;
    it = (dream_idx_type) rint((t - delay)/dt);

    // Check if index is out of bounds.
    if ((it < nt) && (it >= 0)) {
      h[it] += ai;
    } else  {
      if  (it >= 0)
        err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
      else
        err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

      if ( (err_level == ErrorLevel::parallel_stop) || (err_level == ErrorLevel::stop) )
        return err; // Bail out.
    }
    xs += dx;
  }

  return err;
}

ErrorLevel dreamline(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                     double xo, double yo, double zo, double a,
                     double dx, double dy, double dt, dream_idx_type nt, double delay, double v,
                     double cp, double *h, ErrorLevel err_level)
{
  dream_idx_type i, it;
  double t;
  double ai;
  double r;
  ErrorLevel err=ErrorLevel::none;
  double xsmax = a/2.0;
  double xsmin = -a/2.0;

  // dy = width;
  double ds = dx * dy;

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  double xs = xsmin + dx / 2.0;

  while (xs <= xsmax) {

    r = sqrt((xo-xs)*(xo-xs) + yo*yo + zo*zo);

    ai = v * ds / (2*M_PI * r);
    ai /= dt;
    ai *= 1000.0;     // Convert to SI units.

    // Propagation delay in micro seconds.
    t = r * 1000.0/cp;
    it = (dream_idx_type) rint((t - delay)/dt);

    // Check if index is out of bounds.
    if ((it < nt) && (it >= 0)) {

      att.att(xc_vec, x_vec, r, it, h, ai);

    } else  {
      if  (it >= 0)
        err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
      else
        err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

      if ( (err_level == ErrorLevel::parallel_stop) || (err_level == ErrorLevel::stop) )
        return err; // Bail out.
    }
    xs += dx;
  }

  return err;
}
