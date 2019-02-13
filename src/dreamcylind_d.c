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


#include <math.h>
#include <stdio.h>
#include "dreamcylind_d.h"
#include "att.h"
#include "dream_error.h"

#if defined(_MSC_VER) || defined(__LCC__)
#include "msvc_rint.h"
#endif

//
//  Function prototypes.
//

void cyl_cart(double x,double y,double R,double haut,double xo,double yo,double zo,double *RESTRICT rj, double *RESTRICT du);

/***
 *
 * subroutine dreamcylind_d
 *
 ***/

int dreamcylind_d(double xo,
                  double yo,
                  double zo,
                  double a,
                  double b,
                  double R,
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
  double phisj, haut;
  dream_idx_type i, it;
  double t, xsmin, xsmax, ai, phi, ds, pi, du, ri;
  double phismin, phismax, dphi, xsi, ysj;
  int err = NONE;

  if (b > 2*R) {
    dream_err_msg("Error in dreamcylind_d: the y-size, b, must be less than the diameter 2R!\n");
  }

  b /=2.0;
  pi = 4.0 * atan( (double) 1.0);
  haut = R - sqrt(R*R - b*b);

  xsmin = -a/2.0;
  xsmax = a/2.0;

  phi = 2.0 * atan(b/(R-haut));
  phismin = -phi/2.0;
  phismax = phi/2.0;

  // Disc on y in rad.
  dphi = asin(dy/R);
  //dphi = dy/r;
  ds = R * dx * dphi;
  //ds = dx * dy; // Approx the same as  ds = R * dx * dphi.

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  //j = 0;
  //j++;
  //phisj = phismin + (j-1)*dphi + dphi/2;
  phisj = phismin + dphi/2;
  ysj = R * sin(phisj);
  while (phisj <= phismax) {

    //i = 0;
    //i++;
    //xsi = xsmin + (i-1)*dx + dx/2;
    xsi = xsmin + dx/2;
    while (xsi <= xsmax) {

      // Compute ri and ds.
      cyl_cart(xsi, ysj, R, haut, xo, yo, zo, &ri, &du);
      ai = v * ds * du/(2*pi * ri);
      ai /= dt;
      // Convert to SI units.
      ai *= 1000;
      // Propagation delay in micro seconds.
      t = ri * 1000/cp;
      it = (dream_idx_type) rint((t - delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {

        // Check if absorbtion is present.
        if (alfa == (double) 0.0) {
          h[it] += ai;
        } else {
          att(alfa,ri,it,dt,cp,h,nt,ai);
        }
      }
      else   {
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

    //j++;
    //phisj = phismin + (j-1) * dphi + dphi/2;
    phisj += dphi;
    ysj = R * sin(phisj);
  }

  return err;
}  /* dreamcylind_d */


/***
 *
 * cyl_cart.
 *
 ***/

void cyl_cart(double x,double y,double R,double haut,double xo,double yo,double zo,double *RESTRICT rj, double *RESTRICT du)
{
  double  z, rx, ry, rz;
  double cotetj;

  *du = (double) 1.0;

  z = R - sqrt(R*R - y*y);
  z = -z;
  rx = xo - x;
  ry = yo - y;
  rz = zo - z;
  *rj = sqrt(rx*rx + ry*ry + rz*rz);

  // Cos theta = -[rx ry rz]*[x y (z-r)]' / (R * *rj) (Scalar prod.). (focused)
  // Cos theta = [rx ry rz]*[x y (z+r)]' / (R * *rj) (Scalar prod.). (defocused)
  cotetj = (rx*x + ry*y + rz*(z+R)) / (R * *rj);

  if (cotetj < (double) 0.0) {
    *du = (double) 0.0;
  }

  return;
} /* cyl_cart */
