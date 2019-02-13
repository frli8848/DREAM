/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2014 Fredrik Lingvall
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

// $Revision: 815 $ $Date: 2015-04-28 16:59:02 +0200 (Tue, 28 Apr 2015) $ $LastChangedBy: frli8848 $

#include <math.h>
#include <stdio.h>
#include "dreamrect.h"
#include "att.h"
#include "dream_error.h"

#if defined(_MSC_VER) || defined(__LCC__)
#include "msvc_rint.h"
#endif

/***
 *
 * subroutine dreamrect
 *
 ***/

int dreamrect(double xo,
               double yo,
               double zo,
               double a,
               double b,
               double dx,
               double dy,
               double dt,
               dream_idx_type nt,
               double delay,
               double v,
               double cp,
               double alfa,
               double *RESTRICT h,
               int err_level)
{
  dream_idx_type i, it, m, n, M, N;
  double t;
  double ai, ds, pi;
  double ri, xsi, ysj;
  double xsmin, xsmax, ysmin, ysmax;
  int err = NONE;
  double rx,ry,rz;

  xsmin = -a/2.0;
  xsmax =  a/2.0;
  ysmin = -b/2.0;
  ysmax =  b/2.0;

  pi = atan( (double) 1.0) * 4.0;
  ds = dx * dy;

  for (i = 0; i < nt; i++)
    h[i] = (double) 0.0;

  M = (dream_idx_type) (b/dy);
  N = (dream_idx_type) (a/dx);

  // Check if absorbtion is present.
  if (alfa == (double) 0.0) {

    rz = zo;
    ysj = ysmin + dy / 2.0;

    //for(m=0; m<M; m++) {
      while (ysj <= ysmax) {
      ry = yo - ysj;
      xsi = xsmin + dx / 2.0;

      //for(n=0; n<N; n++) {
        while (xsi <= xsmax) {

        rx = xo - xsi;
        ri = sqrt(rx*rx + ry*ry + rz*rz);

        ai = v * ds / (2*pi * ri);
        ai /= dt;
        ai *= 1000.0;		// Convert to SI units.

        t = ri * 1000.0/cp;	// Propagation delay in micro seconds.
        it = (dream_idx_type) rint((t - delay)/dt); // Sample index.

        // Check if index is out of bounds.
        if ( (it < nt) && (it >= 0) ) {
          h[it] += ai;
        }
        else {
          if  (it >= 0)
            err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
          else
            err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

          if ( (err_level == PARALLEL_STOP) || (err_level == STOP) )
            return err; // Bail out.
        }
        xsi += dx;
      }
      ysj += dy;
    }

  } else { // Absorbtion.

    rz = zo;
    ysj = ysmin + dy / 2.0;
    while (ysj <= ysmax) {
      ry = yo - ysj;
      xsi = xsmin + dx / 2.0;
      while (xsi <= xsmax) {

        rx = xo - xsi;
        ri = sqrt(rx*rx + ry*ry + rz*rz);

        ai = v * ds / (2*pi * ri);
        ai /= dt;
        ai *= 1000.0;		// Convert to SI units.

        t = ri * 1000.0/cp;	// Propagation delay in micro seconds.
        it = (dream_idx_type) rint((t - delay)/dt); // Sample index.

        // Check if index is out of bounds.
        if ( (it < nt) && (it >= 0) ) {
          att(alfa,ri,it,dt,cp,h,nt,ai);
        }
        else {
          if  (it >= 0)
            err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
          else
            err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

          if ( (err_level == PARALLEL_STOP) || (err_level == STOP) )
            return err; // Bail out.
        }
        xsi += dx;
      }
      ysj += dy;
    }

  } // if (alfa == ...)

  return err;
} /* dreamrect */
