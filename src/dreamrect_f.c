/***
*
* Copyright (C) 2003,2006,2007,2008,2009,2014 Fredrik Lingvall
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
#include "dreamrect.h"
#include "att.h"
#include "dream_error.h"

#if defined(_MSC_VER) || defined(__LCC__)
#include "msvc_rint.h"
#endif

//
//  Function prototypes.
//

void distance(double xo,
           double yo,
           double zo,
           double x,
           double y,
           double *RESTRICT ri);

void focusing(int ifoc, double focal,
              double xs, double ys,
              double xamax, double yamax, double ramax,
              double cp, double *RESTRICT retfoc);

/***
 *
 * subroutine dreamrect_f - Focused rectangular transducer.
 *
 ***/

int dreamrect_f(double xo,
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
                 double alpha,
                 int ifoc, double focal,
                 double *RESTRICT h,
                 int err_level)
{
  dream_idx_type i, j, it;
  double t;
  double ai, ds, pi, c, retfoc;
  double ri, xsi, ysj;
  double xsmin = -a/2;
  double xsmax = a/2;
  double ysmin = -b/2;
  double ysmax = b/2;
  int err = NONE;

  pi = atan( (double) 1.0) * 4.0;
  ds = dx * dy;

  //b = (-ysmin + ysmax) / 2;
  //a = (-xsmin + xsmax) / 2;
  c = sqrt(a*a + b*b);

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  j = 0;
  j++;
  ysj = ysmin + (j-1) * dy + dy / 2;
  while (ysj <= ysmax) {

    i = 0;
    i++;
    xsi = xsmin + (i-1) * dx + dx / 2;
    while (xsi <= xsmax) {

      distance(xo, yo, zo, xsi, ysj,  &ri);
      focusing(ifoc, focal, xsi, ysj, a, b, c, cp, &retfoc);

      ai = v * ds / (2*pi * ri);
      ai /= dt;
      // Convert to SI units.
      ai *= 1000;
      // Propagation delay in micro seconds.
      t = ri * 1000/cp;
      it = (dream_idx_type) rint((t - delay + retfoc)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {

        // Check if absorbtion is present.
        if (alpha == (double) 0.0) {
          h[it] += ai;
        } else {
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

      i++;
      xsi = xsmin + (i-1)*dx + dx/2;
    }

    j++;
    ysj = ysmin + (j-1)*dy + dy / 2;
  }

  return err;
} /* dreamrect */


/* ****************************************************** */

/***
 *
 * subrutine distance(xi,xs,hs,ri,rx,rz) pour trouver le longeur du vecteur
 *
 ***/

void distance(double xo,
           double yo,
           double zo,
           double x,
           double y,
           double *RESTRICT ri)
{
  double rx, ry, rz;

  ry = yo - y;
  rx = xo - x;
  rz = zo;
  *ri = sqrt(rx*rx + rz*rz + ry*ry);

  return;
} /* distance */


/***
 *
 * subroutine focussing gives le retard retfoc du au focussing.
 *
 ***/

void focusing(int ifoc, double focal,
             double xs, double ys,
             double xamax, double yamax, double ramax,
             double cp, double *RESTRICT retfoc)
{
  double diff, rmax, retx, rety;

  /* ifoc =1 - no foc, 2 foc x ,3 foc y, 4 foc xy 5 foc x+y */
  switch (ifoc) {

  case 1:
    return;

  case 2:
    rmax = sqrt(xamax*xamax + focal*focal);
    diff = rmax - sqrt(xs*xs + focal*focal);
    *retfoc = diff * 1000 / cp;
    break;

  case 3:
    rmax = sqrt(yamax*yamax + focal*focal);
    diff = rmax - sqrt(ys*ys + focal*focal);
    *retfoc = diff * 1000 / cp;
    break;

  case 4:
    rmax = sqrt(ramax*ramax + focal*focal);
    diff = rmax - sqrt(xs*xs + ys*ys + focal*focal);
    *retfoc = diff * 1000 / cp;
    break;

  case 5:
    rmax = sqrt(ramax*ramax + focal*focal);
    retx = sqrt(xs*xs + focal*focal);
    rety = sqrt(ys*ys + focal*focal);
    diff = rmax - (retx + rety);
    *retfoc = diff * 1000 / cp;

  default:
    break;

  } // switch

  return;
} /* focusing */