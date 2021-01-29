/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2012,2014 Fredrik Lingvall
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
#include "att.h"
#include "dream_error.h"

#if defined(_MSC_VER) || defined(__LCC__)
#include "msvc_rint.h"
#endif

/***
 *
 *  dreamline - Spatial impulse respone of line (slit)
 *
 ***/

int dreamline(double xo, double yo, double zo, double a,
              double dx, double dy, double dt, dream_idx_type nt, double delay, double v,
              double cp, double alpha, double *RESTRICT h, int err_level)
{
  dream_idx_type i, it;
  double t;
  double ai;
  double ds, pi, ri;
  int    err=NONE;
  double x;
  double xsmax = a/2.0;
  double xsmin = -a/2.0;

  pi = 4.0 * atan(1.0);
  // dy = width;
  ds = dx * dy;

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  x = xsmin + dx / 2.0;

  while (x <= xsmax) {

    //distance(xo, yo, zo, x, ys, zs, &ri);
    ri = sqrt((xo-x)*(xo-x) + yo*yo + zo*zo);

    ai = v * ds / (2*pi * ri);
    ai /= dt;
    ai *= 1000.0;     // Convert to SI units.

    // Propagation delay in micro seconds.
    t = ri * 1000.0/cp;
    it = (dream_idx_type) rint((t - delay)/dt);

    // Check if index is out of bounds.
    if ((it < nt) && (it >= 0)) {
      // Check if absorbtion is present.
      if (alpha == (double) 0.0) {
        h[it] += ai;
      }
      else {
        att(alpha,ri,it,dt,cp,h,nt,ai);
      }
    }
    else  {
      if  (it >= 0)
        err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
      else
        err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

      if ( (err_level == PARALLEL_STOP) || (err_level == STOP) )
        return err; // Bail out.
    }
    x += dx;
  }

  return err;
}
