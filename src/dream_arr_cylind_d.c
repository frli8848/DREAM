/***
*
* Copyright (C) 2002,2006,2007,2008,2009,2014 Fredrik Lingvall
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

// $Revision: 767 $ $Date: 2014-05-15 20:44:36 +0200 (Thu, 15 May 2014) $ $LastChangedBy: frli8848 $

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
 *  2D Array with cylindrical concave elements.
 *
 ***/

//
// Function prototypes
//

int cylind_d(double xo, double yo, double zo, double xs, double ys, double zs, double R, double a, double b,
             double dx, double dy, double dt, dream_idx_type nt, double delay, double retfoc, double retsteer,
             double v, double cp, double alfa,  double weight, double *RESTRICT h, int err_level);

void cyl_cart(double xs, double ys, double zs, double r,double haut,double
              xo,double yo,double zo,double *RESTRICT rj, double *RESTRICT du);

/***
 *
 * Subroutine dream_arr_cylind_d
 *
 ***/

int dream_arr_cylind_d(double xo, double yo, double zo, double a, double b, double R,
                       double dx, double dy, double dt,
                       dream_idx_type nt, double delay, double v, double cp, double alfa, int isize,
                       double *RESTRICT gx, double *RESTRICT gy, double *RESTRICT gz, int ifoc, double focal, int ister,
                       double theta, double phi, double *RESTRICT apod, int iweight, int iapo, double param,
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

  maxdimarr(&xamax, &yamax, &ramax, gx, gy, gz, isize);

  for (i=0; i<isize; i++) {
    center_pos(&xs, &ys, &zs, i, gx, gy, gz);
    focusing(ifoc, focal, xs, ys, xamax, yamax, ramax, cp, &retfoc);
    beamsteering(ister, theta, phi, xs, ys, xamax, yamax, ramax, cp, &retsteer);
    weighting(iweight, iapo, i, apod, &weight, xs, ys, ramax, param, isize);

    err = cylind_d(xo,yo,zo,xs,ys,zs,R,a,b,dx,dy,dt,nt,delay,retfoc,retsteer,v,cp,alfa,weight,h,err_level);
    if (err != NONE)
      out_err = err;

    superpoz(h, ha, nt);
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
                         dream_idx_type nt, double delay, double v, double cp, double alfa, int isize,
                         double *RESTRICT gx, double *RESTRICT gy, double *RESTRICT gz, int ifoc, double *RESTRICT focal,
                         int ister, double theta, double phi,
                         double *RESTRICT apod, int iweight, int iapo, double param, double *RESTRICT ha, int err_level)
{
  double retsteer;
  double *RESTRICT h;
  dream_idx_type i;
  double ramax, xamax, yamax;
  double xs, ys, zs, retfoc, weight;
  int err = NONE;

  h = (double*) malloc(nt*sizeof(double));

  for (i=0; i<nt; i++) {
    ha[i] = (double) 0.0;
  }

  retfoc   = (double) 0.0;
  retsteer = (double) 0.0;
  weight   = (double) 1.0;

  maxdimarr(&xamax, &yamax, &ramax, gx, gy, gz, isize);

  for (i=0; i<isize; i++) {
    center_pos(&xs, &ys, &zs, i, gx, gy, gz);
    focusing(ifoc, focal[i], xs, ys, xamax, yamax, ramax, cp, &retfoc); // Note ifoc must be 6 here!
    beamsteering(ister, theta, phi, xs, ys, xamax, yamax, ramax, cp, &retsteer);
    weighting(iweight, iapo, i, apod, &weight, xs, ys, ramax, param, isize);

    err = cylind_d(xo,yo,zo,xs,ys,zs,R,a,b,dx,dy,dt,nt,delay,retfoc,retsteer,v,cp,alfa,weight,h,err_level);

    superpoz(h, ha, nt);
  }

  free(h);

  return err;
} /* dream_arr_cylind_udd */

/********************************************************************/


/***
 *
 * cylind_d
 *
 ***/

int cylind_d(double xo, double yo, double zo, double xs, double ys, double zs, double R, double a, double b,
             double dx, double dy, double dt, dream_idx_type nt, double delay, double retfoc, double retsteer,
             double v, double cp, double alfa,  double weight, double *RESTRICT h, int err_level)
{
  double phisj, haut;
  dream_idx_type i, it;
  double t, xsmin, xsmax, ai, phi, ds, pi, du, ri;
  double phismin, phismax, dphi, xsi, ysj;
  double decal;
  int err = NONE;

  if (b > 2*R)
    dream_err_msg("Error in dream_arr_cylind_d: the y-width, b, must be less than the diameter 2R!");

  b /=2;
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

  for (i = 0; i < nt; i++)
    h[i] = (double) 0.0 ;

  //j = 0;
  //j++;
  //phisj = phismin + (j-1)*dphi + dphi/2;
  phisj = phismin + dphi/2.0;
  ysj = R * sin(phisj);
  while (phisj <= phismax) {

    //i = 0;
    //i++;
    //xsi = xsmin + (i-1)*dx + dx/2;
    xsi = xsmin + dx/2.0;
    while (xsi <= xsmax) {

      // Compute ri and ds.
      cyl_cart(xsi, ysj, zs, R, haut, xo, yo, zo, &ri, &du);
      ai = v * ds * du/(2*pi * ri);
      ai /= dt;
      // Convert to SI units.
      ai *= 1000;
      // Propagation delay in micro seconds.
      t = ri * 1000/cp;
      it = (dream_idx_type) rint((t - delay + decal)/dt);

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
} /* cylind_d */


/***
 *
 * cyl_cart.
 *
 ***/

void cyl_cart(double xs, double ys, double zs, double R,double haut,double xo,
              double yo,double zo,double *RESTRICT rj, double *RESTRICT du)
{
  double z, rx, ry, rz;
  double cotetj;

  *du = (double) 1.0;

  z = zs + R - sqrt(R*R - ys*ys);
  z = -z;
  rx = xo - xs;
  ry = yo - ys;
  rz = zo - z;
  *rj = sqrt(rx*rx + ry*ry + rz*rz);

  // Cos theta = -[rx ry rz]*[x y (z-r)]' / (r * *rj) (Scalar prod.). (focused)
  // Cos theta = [rx ry rz]*[x y (z+r)]' / (r * *rj) (Scalar prod.). (defocused)
  cotetj = (rx*xs + ry*ys + rz*(z+R)) / (R * *rj);

  if (cotetj < (double) 0.0) {
    *du = (double) 0.0;
  }

  return;
} /* cyl_cart */
