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

#include <math.h>

#include "dreamrect_f.h"

//
//  Function prototypes.
//

double focus_delay_rect(FocusMet foc_met, double focal,
                double xs, double ys,
                double xamax, double yamax, double ramax,
                double cp);

/***
 *
 * dreamrect_f - Focused rectangular transducer.
 *
 ***/

ErrorLevel dreamrect_f(double xo, double yo, double zo,
                       double a, double b,
                       FocusMet foc_met, double focal,
                       double dx, double dy, double dt, dream_idx_type nt,
                       double delay,
                       double v, double cp,
                       double *h,
                       ErrorLevel err_level)
{
  ErrorLevel err = ErrorLevel::none;

  // x-dim and y-dim integration bounds
  double xsmin = -a/2;
  double xsmax = a/2;
  double ysmin = -b/2;
  double ysmax = b/2;

  double ds = dx * dy;
  double c = sqrt(a*a + b*b);

  for (dream_idx_type i = 0; i < nt; i++) {
    h[i] = 0.0;
  }

  double ys = ysmin + dy/2.0;
  while (ys <= ysmax) {

    double ry = yo - ys;

    double xs = xsmin + dx/2.0;
    while (xs <= xsmax) {

      // Compute the distance (length) from an observation point (xo,yo,zo)
      // to a point (xs,ys) on the transducer surface.
      double rx = xo - xs;
      double r = sqrt(rx*rx + ry*ry + zo*zo);

      double foc_delay = focus_delay_rect(foc_met, focal, xs, ys, a, b, c, cp);

      double ai = v * ds / (2*M_PI * r);
      ai /= dt;
      ai *= 1.0e3; // Convert to SI units.

      // Propagation delay in micro seconds.
      double t = r * 1.0e3/cp;
      dream_idx_type it = (dream_idx_type) rint((t - delay + foc_delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {
        h[it] += ai;
      } else {

        if (it >= 0) {
          err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
        } else {
          err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);
        }

        if ( (err_level == ErrorLevel::parallel_stop) || (err_level == ErrorLevel::stop) ) {
          return err; // Bail out.
        }
      }

      xs += dx;
    }
    ys += dy;
  }

  return err;
}

ErrorLevel dreamrect_f(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                       double xo, double yo, double zo,
                       double a, double b, FocusMet foc_met, double focal,
                       double dx, double dy, double dt,
                       dream_idx_type nt,
                       double delay,
                       double v, double cp,
                       double *h,
                       ErrorLevel err_level)
{
  // x-dim and y-dim integration bounds
  double xsmin = -a/2;
  double xsmax = a/2;
  double ysmin = -b/2;
  double ysmax = b/2;
  ErrorLevel err = ErrorLevel::none;

  double ds = dx * dy;
  double c = sqrt(a*a + b*b);

  for (dream_idx_type i = 0; i < nt; i++) {
    h[i] = 0.0;
  }

  double ys = ysmin + dy/2.0;

  while (ys <= ysmax) {

    double ry = yo - ys;

    double xs = xsmin + dx/2.0;

    while (xs <= xsmax) {

      // Compute the distance (length) from an observation point (xo,yo,zo)
      // to a point (xs,ys) on the transducer surface.
      double rx = xo - xs;
      double r = sqrt(rx*rx + ry*ry + zo*zo);

      double foc_delay = focus_delay_rect(foc_met, focal, xs, ys, a, b, c, cp);

      double ai = v * ds / (2*M_PI * r);
      ai /= dt;
      ai *= 1.0e3; // Convert to SI units.

      // Propagation delay in micro seconds.
      double t = r * 1.0e3/cp;
      dream_idx_type it = (dream_idx_type) rint((t - delay + foc_delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {
        att.att(xc_vec, x_vec, r, it, h, ai);
      } else {

        if (it >= 0) {
          err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
        } else {
          err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);
        }

        if ( (err_level == ErrorLevel::parallel_stop) || (err_level == ErrorLevel::stop) ) {
          return err; // Bail out.
        }
      }

      xs += dx;
    }
    ys += dy;
  }

  return err;
}

/***
 *
 *  Computes the focus delay
 *
 ***/

double focus_delay_rect(FocusMet foc_met, double focal,
                        double xs, double ys,
                        double xamax, double yamax, double ramax,
                        double cp)
{
  double foc_delay = 0.0;

  switch (foc_met) {

  case FocusMet::none:
    break;

  case FocusMet::x:
    {
      double rmax = sqrt(xamax*xamax + focal*focal);
      double diff = rmax - sqrt(xs*xs + focal*focal);
      foc_delay = diff * 1.0e3 / cp;
    }
    break;

  case FocusMet::y:
    {
      double rmax = sqrt(yamax*yamax + focal*focal);
      double diff = rmax - sqrt(ys*ys + focal*focal);
      foc_delay = diff * 1.0e3 / cp;
    }
    break;

  case FocusMet::xy:
    {
      double rmax = sqrt(ramax*ramax + focal*focal);
      double diff = rmax - sqrt(xs*xs + ys*ys + focal*focal);
      foc_delay = diff * 1.0e3 / cp;
    }
    break;

  case FocusMet::x_y:
    {
      double rmax = sqrt(ramax*ramax + focal*focal);
      double delay_x = sqrt(xs*xs + focal*focal);
      double delay_y = sqrt(ys*ys + focal*focal);
      double diff = rmax - (delay_x + delay_y);
      foc_delay = diff * 1.0e3 / cp;
    }
    break;

  default:
    break;

  }

  return foc_delay;             // [us]
}
