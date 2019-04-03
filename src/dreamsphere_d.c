/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2014,2019 Fredrik Lingvall
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
#include "dreamsphere_d.h"
#include "att.h"
#include "dream_error.h"

#if defined(_MSC_VER) || defined(__LCC__)
#include "msvc_rint.h"
#endif

//
// Function prototypes.
//

void sph_cart_d(double xi, double yi, double dx, double dy, double R, double haut,
              double xoi, double yoi, double zoi, double *RESTRICT rj, double *RESTRICT cotetj, double *RESTRICT du);
void xlimit_sphere_d(double yi, double r, double xs, double ys, double *RESTRICT xsmin, double *RESTRICT xsmax);

/***
 *
 *  dreamsphere_d
 *
 *  Computes the (spatial) impulse response for a de-focused (convex) spherical aperture.
 *
 ***/

int dreamsphere_d(double xo, double yo, double zo, double r, double R, double
                 dx, double dy, double dt, dream_idx_type nt,
                 double  delay, double v, double cp, double alpha, double  *h, int err_level)
{
  double haut;
  dream_idx_type i;
  double t, cotet, xsmin, ysmin, xsmax, ysmax, ai, ds, pi, ri;
  dream_idx_type it;
  double xs, ys;
  double x, y;
  int err = NONE;

  pi = atan( (double) 1.0) * 4.0;
  haut = R - sqrt(R*R - r*r);
  xs = (double) 0.0;
  ys = (double) 0.0;
  ysmin = -r + ys;
  ysmax =  r + ys;

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  y = ysmin + dy / 2.0;
  while (y <= ysmax) {

    xlimit_sphere_d(y, r, xs, ys, &xsmin, &xsmax);

    x = xsmin + dx / 2.0;
    while (x <= xsmax) {

      sph_cart_d(x,y,dx,dy,R,haut,xo,yo,zo,&ri,&cotet,&ds);

      ai = v * ds / (2*pi * ri);
      ai /= dt;
      ai *= 1000; // Convert to SI units.

      // Propagation delay in micro seconds.
      t = ri * 1000/cp;
      it = (dream_idx_type) rint((t - delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {

        // Check if absorbtion is present.
        if (alpha == (double) 0.0) {
          h[it] += ai;
        } else {
          att(alpha,ri,it,dt,cp,h,nt,ai);
        }
      } else  {

        if  (it >= 0)
          err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
        else
          err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

        if ( (err_level == PARALLEL_STOP) || (err_level == STOP) )
          return err; // Bail out.
      }

      x += dx;
    }

    y += dy;
  }

  return err;
}

/***
 *
 * xlimit_sphere_d
 *
 * Computes the x-axis intregration limits.
 *
 ***/

void  xlimit_sphere_d(double yi, double r, double xs, double ys, double *RESTRICT xsmin, double *RESTRICT xsmax)
{
  double rs;

  rs = r*r - (ys-yi) * (ys-yi);
  rs = sqrt(rs);
  *xsmin = -rs + xs;
  *xsmax =  rs + xs;

  return;
}

/***
 *
 * sph_cart_d.
 *
 ***/

void sph_cart_d(double xi, double yi, double dx, double dy, double R, double haut,
              double xoi, double yoi, double zoi, double *RESTRICT rj, double *RESTRICT cotetj, double *RESTRICT du)
{
  double a, b, d, e, x, y, z, z1, z2, z3, z4, rx, ry, rz, az1, az2;
  double az3, az4;

  *cotetj = (double) 0.0;
  *rj = (double) 0.0;
  *du = (double) 0.0;

  x = xi;
  y = yi;
  z = sqrt(R * R - x*x - y*y);
  z = R - z;

  z = -z;

  rx = xoi - x;
  ry = yoi - y;
  rz = zoi - z;
  *rj = sqrt(rx * rx + ry * ry + rz * rz);

  z1 = R - sqrt(R*R - xi*xi - yi*yi);
  z2 = R - sqrt(R*R - (xi + dx) * (xi + dx) - yi*yi);
  z3 = R - sqrt(R*R - xi*xi - (yi + dy) * (yi + dy));
  z4 = R - sqrt(R*R - (xi + dx) * (xi + dx) - (yi + dy) * (yi + dy));

  z1 = -z1;
  z2 = -z2;
  z3 = -z3;
  z4 = -z4;

  az1 = fabs(z2 - z1);
  az2 = fabs(z3 - z1);
  az3 = fabs(z4 - z3);
  az4 = fabs(z4 - z2);
  a = sqrt(dx*dx + az1*az1);
  b = sqrt(dy*dy + az2*az2);
  e = sqrt(dx*dx + az3*az3);
  d = sqrt(dy*dy + az4*az4);
  *du = a * b / (double) 2.0 + e * d / (double) 2.0;

  *cotetj = (rx*x + ry*y + rz * (z + R)) / (R * *rj);

  if (*cotetj < (double)0.) {
    *du = (double) 0.0;
  }

  return;
}
