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
#include "dreamsphere_f.h"
#include "att.h"
#include "dream_error.h"

#if defined(_MSC_VER) || defined(__LCC__)
#include "msvc_rint.h"
#endif

//
// Function prototypes.
//

void sph_cart(double xi, double yi, double dx, double dy, double r, double haut,
              double xoi, double yoi, double zoi, double *RESTRICT rj, double *RESTRICT cotetj, double *RESTRICT du);
void xlimit(double yi, double r, double xs, double ys, double *RESTRICT xsmin, double *RESTRICT xsmax);

/***
 *
 *  spher_f - pour calculer pulse response of a sperical focussing aperture.
 *
 ***/

int dreamsphere_f(double xo, double yo, double zo, double r, double R, double dx,
                 double dy, double dt, dream_idx_type nt,
                 double  delay, double v, double cp, double alpha ,double  *h, int err_level)
{
  double haut;
  dream_idx_type i; // j;
  double t, cotet, xsmin, ysmin, xsmax, ysmax, ai, ds, pi, ri;
  dream_idx_type it;
  double xs, ys;
  double xsi, ysj;
  int  err = NONE;

  pi = atan( (double) 1.0) * 4.0;
  haut = R - sqrt(R*R - r*r);

  // Position transducer at (0,0).
  xs = (double) 0.0;
  ys = (double) 0.0;

  ysmin = -r;
  ysmax =  r;

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  //j=0;
  //j++;
  //ysj = ysmin + (j-1) * dy + dy / 2;
  ysj = ysmin + dy / 2.0;
  while (ysj <= ysmax) {

    xlimit(ysj,r,xs,ys, &xsmin, &xsmax);

    //i=0;
    //i++;
    //xsi = xsmin + (i-1) * dx + dx / 2;
    xsi = xsmin + dx / 2.0;
    while (xsi <= xsmax) {

      sph_cart(xsi,ysj,dx,dy,R,haut,xo,yo,zo,&ri,&cotet,&ds);

      ai = v * ds / (2*pi * ri);
      ai /= dt;
      // Convert to SI units.
      ai *= 1000;
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
      //xsi = xsmin + (i-1)*dx + dx/2;
      xsi += dx;
    }

    //j++;
    //ysj = ysmin + (j-1)*dy + dy / 2;
    ysj += dy;
  }

  return err;
} /* dreamspher_f */


/***
 *
 * call xlimit(ysj,a,xs,ys,xsmin,xsmax)
 * subroutine xlimit - pour definir les limits d integration en x
 *
 ***/

void  xlimit(double yi, double r, double xs, double ys, double *RESTRICT xsmin, double *RESTRICT xsmax)
{
  double rs;

  rs = sqrt(r*r - (ys-yi) * (ys-yi));
  *xsmin = -rs + xs;
  *xsmax =  rs + xs;

  return;
} /* xlimit */


/***
 *
 * sph_cart.
 *
 ***/

void sph_cart(double xi, double yi, double dx, double dy, double R, double haut,
              double xoi, double yoi, double zoi, double *RESTRICT rj, double *RESTRICT cotetj, double *RESTRICT du)
{
  double a, b, d, e, x, y, z, z1, z2, z3, z4, rx, ry, rz, az1, az2;
  double az3, az4, dis1=0, dis2=0;

  *cotetj = (double) 0.0;
  *rj = (double) 0.0;
  *du = (double) 0.0;

  if (zoi <= haut) {
    dis1 = sqrt(xoi*xoi + yoi*yoi);
    dis2 = sqrt(zoi * (R * (double) 2.0 - zoi));
  }

  x = xi;
  y = yi;
  z = sqrt(R*R - x*x - y*y);
  z = R - z;

  rx = xoi - x;
  ry = yoi - y;
  rz = zoi - z;
  *rj = sqrt(rx*rx + ry*ry + rz*rz);

  z1 = R - sqrt(R*R - xi*xi - yi*yi);
  z2 = R - sqrt(R*R - (xi+dx)*(xi+dx) - yi*yi);
  z3 = R - sqrt(R*R - xi*xi - (yi+dy)*(yi+dy));
  z4 = R - sqrt(R*R - (xi+dx)*(xi+dx) - (yi+dy)*(yi+dy));

  az1 = fabs(z2 - z1);
  az2 = fabs(z3 - z1);
  az3 = fabs(z4 - z3);
  az4 = fabs(z4 - z2);
  a = sqrt(dx*dx + az1*az1);
  b = sqrt(dy*dy + az2*az2);
  e = sqrt(dx*dx + az3*az3);
  d = sqrt(dy*dy + az4*az4);
  *du = a * b / (double) 2.0 + e * d / (double) 2.0;

  *cotetj = -(rx*x + ry*y + rz*(z-R)) / (R * *rj);

  if (zoi <= haut) {
    if (*cotetj < (double) 0.0 || dis1 > dis2) {
      *du = (double) 0.0;
    }
  } else if (zoi > haut) {
    if (*cotetj < (double) 0.0) {
      *du = (double) 0.0;
    }
  }

  return;
} /* sph_cart */