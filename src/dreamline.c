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
 *      Spatial impulse respone of line (slit)
 *
 ***/

/***
 *
 *  dreamline
 *
 ***/

int dreamline(double xo, double yo, double zo, double a,
              double dx, double dy, double dt, dream_idx_type nt, double delay, double v,
              double cp, double alfa, double *RESTRICT h, int err_level)
{
  dream_idx_type i, it;
  double t;
  double ai;
  double ds, pi, ri;
  int    err=NONE;
  double xsi;
  double xsmax = a/2.0;
  double xsmin = -a/2.0;

  pi = atan((double) 1.0) * 4.0;
  // dy = width;
  ds = dx * dy;

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  //i = 0;
  //i++;
  //xsi = xsmin + (i-1) * dx + dx / 2;
  xsi = xsmin + dx / 2.0;
  while (xsi <= xsmax) {

    //modri(xo, yo, zo, xsi, ys, zs, &ri);
    //rx = xo - xsi;
    //ry = yo - ys;
    //rz = zo - zs;
    //ri = sqrt(rx*rx + yo*yo + zo*zo);
    ri = sqrt((xo-xsi)*(xo-xsi) + yo*yo + zo*zo);

    ai = v * ds / (2*pi * ri);
    ai /= dt;
    // Convert to SI units.
    ai *= 1000.0;
    // Propagation delay in micro seconds.
    t = ri * 1000.0/cp;
    it = (dream_idx_type) rint((t - delay)/dt);

    // Check if index is out of bounds.
    if ((it < nt) && (it >= 0)) {
      // Check if absorbtion is present.
      if (alfa == (double) 0.0) {
        h[it] += ai;
      }
      else {
        att(alfa,ri,it,dt,cp,h,nt,ai);
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

    //i++;
    //xsi = xsmin + (i-1) * dx + dx/2;
    xsi += dx;
  }

  return err;
} /* dreamline */
