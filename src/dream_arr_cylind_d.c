/***
*
* Copyright (C) 2002,2006,2007,2008,2009,2014,2019 Fredrik Lingvall
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
#include <stdlib.h>

#include "dream_arr_cylind_d.h"
#include "att.h"
#include "arr_functions.h"
#include "dream_error.h"

#if defined(_MSC_VER) || defined(__LCC__)
#include "msvc_rint.h"
#endif

/***
 *
 * Array with cylindrical convex elements.
 *
 ***/

//
// Function prototypes
//

int cylind_d(double xo, double yo, double zo, double xs, double ys, double zs, double R, double a, double b,
             double dx, double dy, double dt, dream_idx_type nt, double delay, double retfoc, double retsteer,
             double v, double cp, double alpha,  double weight, double *RESTRICT h, int err_level);

void cyl_arr_d(double xs, double ys, double zs, double R, double haut,
               double xo, double yo, double zo,
               double *RESTRICT r, double *RESTRICT du);

/***
 *
 * Subroutine dream_arr_cylind_d
 *
 ***/

int dream_arr_cylind_d(double xo, double yo, double zo, double a, double b, double R,
                       double dx, double dy, double dt,
                       dream_idx_type nt, double delay, double v, double cp, double alpha, int num_elements,
                       double *RESTRICT gx, double *RESTRICT gy, double *RESTRICT gz, int foc_type, double focal, int ister,
                       double theta, double phi, double *RESTRICT apod, bool do_apod, int apod_type, double param,
                       double *RESTRICT ha, int err_level)
{
  double retsteer;
  double *RESTRICT h;
  dream_idx_type i;
  double ramax, xamax, yamax;
  double xs, ys, zs, retfoc, weight;
  int err = NONE, out_err = NONE;

  h = (double*) malloc(nt*sizeof(double));

  for (i=0; i<nt; i++) {
    ha[i] = (double) 0.0;
  }

  retfoc   = (double) 0.0;
  retsteer = (double) 0.0;
  weight   = (double) 1.0;

  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, num_elements);

  for (i=0; i<num_elements; i++) {
    center_pos(&xs, &ys, &zs, i, gx, gy, gz);
    focusing(foc_type, focal, xs, ys, xamax, yamax, ramax, cp, &retfoc);
    beamsteering(ister, theta, phi, xs, ys, xamax, yamax, ramax, cp, &retsteer);
    if (do_apod) {
      apodization(apod_type, i, apod, &weight, xs, ys, ramax, param);
    }
    err = cylind_d(xo,yo,zo,xs,ys,zs,R,a,b,dx,dy,dt,nt,delay,retfoc,retsteer,v,cp,alpha,weight,h,err_level);
    if (err != NONE)
      out_err = err;

    superpos(h, ha, nt);
  }

  free(h);

  return out_err;
} /* dream_arr_cylind_f */

/***
 *
 * Subroutine dream_arr_cylind_udd - user defined focusing.
 *
 ***/

int dream_arr_cylind_udd(double xo, double yo, double zo, double a, double b, double R,
                         double dx, double dy, double dt,
                         dream_idx_type nt, double delay, double v, double cp, double alpha, int num_elements,
                         double *RESTRICT gx, double *RESTRICT gy, double *RESTRICT gz, int foc_type, double *RESTRICT focal,
                         int ister, double theta, double phi,
                         double *RESTRICT apod, bool do_apod, int apod_type, double param, double *RESTRICT ha, int err_level)
{
  double retsteer;
  double *RESTRICT h;
  dream_idx_type i;
  double ramax, xamax, yamax;
  double xs, ys, zs, retfoc, weight;
  int err = NONE, out_err = NONE;

  h = (double*) malloc(nt*sizeof(double));

  for (i=0; i<nt; i++) {
    ha[i] = (double) 0.0;
  }

  retfoc   = (double) 0.0;
  retsteer = (double) 0.0;
  weight   = (double) 1.0;

  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, num_elements);

  for (i=0; i<num_elements; i++) {
    center_pos(&xs, &ys, &zs, i, gx, gy, gz);
    focusing(foc_type, focal[i], xs, ys, xamax, yamax, ramax, cp, &retfoc); // Note foc_type must be 6 here!
    beamsteering(ister, theta, phi, xs, ys, xamax, yamax, ramax, cp, &retsteer);
    if (do_apod) {
      apodization(apod_type, i, apod, &weight, xs, ys, ramax, param);
    }
    err = cylind_d(xo,yo,zo,xs,ys,zs,R,a,b,dx,dy,dt,nt,delay,retfoc,retsteer,v,cp,alpha,weight,h,err_level);
    if (err != NONE)
      out_err = err;

    superpos(h, ha, nt);
  }

  free(h);

  return out_err;
} /* dream_arr_cylind_udd */

/********************************************************************/


/***
 *
 * cylind_d
 *
 ***/

int cylind_d(double xo, double yo, double zo, double xs, double ys, double zs, double R, double a, double b,
             double dx, double dy, double dt, dream_idx_type nt, double delay, double retfoc, double retsteer,
             double v, double cp, double alpha,  double weight, double *RESTRICT h, int err_level)
{
  double haut;
  dream_idx_type i, it;
  double t, xsmin, xsmax, ai, phi, ds, pi, du, r;
  double phismin, phismax, dphi, x, y;
  double decal;
  int err = NONE;

  if (b > 2*R)
    dream_err_msg("Error in dream_arr_cylind_d: the y-width, b, must be less than the diameter 2R!");

  b /=2.0;
  decal = retfoc + retsteer;
  pi = 4.0 * atan( (double) 1.0);
  haut = R - sqrt(R*R - b*b);

  xsmin = xs - a/2.0;
  xsmax = xs + a/2.0;

  phi = 2.0 * atan(b/(R-haut));
  phismin = -phi/2.0;
  phismax = phi/2.0;

  // Disc on y in rad.
  dphi = asin(dy/R);
  //dphi = dy/R;
  ds = R * dx * dphi;
  //ds = dx * dy; // Approx the same as  ds = R * dx * dphi.

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  phi = phismin + dphi/2.0;
  y = R * sin(phi);
  while (phi <= phismax) {

    x = xsmin + dx/2.0;
    while (x <= xsmax) {

      // Compute r and ds.
      cyl_arr_d(x, y, zs, R, haut, xo, yo, zo, &r, &du);
      ai = v * ds * du/(2*pi * r);
      ai /= dt;
      ai *= 1000; // Convert to SI units.

      // Propagation delay in micro seconds.
      t = r * 1000/cp;
      it = (dream_idx_type) rint((t - delay + decal)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {

        // Check if absorbtion is present.
        if (alpha == 0.0) {
          h[it] += ai;
        } else {
          att(alpha,r,it,dt,cp,h,nt,ai);
        }

      } else {

        if  (it >= 0)
          err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
        else
          err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

        if ( (err_level == PARALLEL_STOP) || (err_level == STOP) )
          return err; // Bail out.
      }

      x += dx;
    }

    phi += dphi;
    y = R * sin(phi);
  }

  return err;
} /* cylind_d */


/***
 *
 * cyl_arr_d.
 *
 ***/

void cyl_arr_d(double xs, double ys, double zs, double R, double haut,double xo,
              double yo,double zo,double *RESTRICT r, double *RESTRICT du)
{
  double z, rx, ry, rz;
  double cotetj;

  *du = 1.0;

  z = zs + R - sqrt(R*R - ys*ys);
  z = -z;
  rx = xo - xs;
  ry = yo - ys;
  rz = zo - z;
  *r = sqrt(rx*rx + ry*ry + rz*rz);

  // Cos theta = -[rx ry rz]*[x y (z-r)]' / (r * *r) (Scalar prod.). (focused)
  // Cos theta = [rx ry rz]*[x y (z+r)]' / (r * *r) (Scalar prod.). (defocused)
  cotetj = (rx*xs + ry*ys + rz*(z+R)) / (R * *r);

  if (zo <= haut) {

    double d1 = sqrt(xo*xo + yo*yo);
    double d2 = sqrt(zo * (R * 2.0 - zo));

    if ( (cotetj < 0.0) || (d1 > d2)) {
      *du = 0.0;
    }
  } else {
    if (cotetj < 0.0) {
      *du = 0.0;
    }
  }
  return;
}
