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
#include <stdlib.h>

#include "dream_error.h"
#include "dream_arr_rect.h"
#include "att.h"
#include "arr_functions.h"

#if defined(_MSC_VER) || defined(__LCC__)
#include "msvc_rint.h"
#endif

//
// Function prototypes
//


int rect_ab(double xo, double yo, double zo, double xs, double ys, double zs, double a, double b,
            double dx, double dy, double dt, dream_idx_type nt, int icheck, double delay, double retfoc,
            double retsteer, double v, double cp, double alfa,  double weight, double *RESTRICT h, int err_level);


/***
 *
 * Subroutine dream_arr_rect
 *
 ***/

int dream_arr_rect(double xo, double yo, double zo, double a, double b, double dx, double dy, double dt,
                  dream_idx_type nt, double delay, double v, double cp, double alfa, int isize,
                  double *RESTRICT gx, double *RESTRICT gy, double *RESTRICT gz, int ifoc, double focal, int ister,
                  double theta, double phi, double *RESTRICT apod, int iweight, int iapo, double param, double *RESTRICT ha,
                  int err_level)
{
  double retsteer;
  double *RESTRICT h;
  dream_idx_type i;
  double ramax, xamax, yamax;
  int icheck;
  double xs, ys, zs, retfoc, weight;
  int err = NONE, out_err = NONE;

  h = (double*) malloc(nt*sizeof(double));

  for (i=0; i<nt; i++) {
    ha[i] = (double) 0.0;
  }

  retfoc   = (double) 0.0;
  retsteer = (double) 0.0;
  weight   = (double) 1.0;
  icheck   = 1;

  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, isize);

  for (i=0; i<isize; i++) {
    center_pos(&xs, &ys, &zs, i, gx, gy, gz);
    focusing(ifoc, focal, xs, ys, xamax, yamax, ramax, cp, &retfoc);
    beamsteering(ister, theta, phi, xs, ys, xamax, yamax, ramax, cp, &retsteer);
    weighting(iweight, iapo, i, apod, &weight, xs, ys, ramax, param, isize);

    err = rect_ab(xo,yo,zo,xs,ys,zs,a,b,dx,dy,dt,nt,icheck,
                 delay,retfoc,retsteer,v,cp,alfa,weight,h,err_level);
    if (err != NONE)
      out_err = err;

    if (icheck == 2) {
      break;
    }
    superpos(h, ha, nt);
  }

  free(h);

  return out_err;
} /* dream_arr_rect */

/***
 *
 * Subroutine dream_arr_rect_ud - user defined focusing.
 *
 ***/

int dream_arr_rect_ud(double xo, double yo, double zo, double a, double b, double dx, double dy, double dt,
                  dream_idx_type nt, double delay, double v, double cp, double alfa, int isize,
                  double *RESTRICT gx, double *RESTRICT gy, double *RESTRICT gz, int ifoc, double *RESTRICT focal, int ister,
                  double theta, double phi, double *RESTRICT apod, int iweight, int iapo, double param, double *RESTRICT ha,
                  int err_level)
{
  double retsteer;
  double *RESTRICT h;
  dream_idx_type i;
  double ramax, xamax, yamax;
  int icheck;
  double xs, ys, zs, retfoc, weight;
  int err = NONE, out_err = NONE;

  h = (double*) malloc(nt*sizeof(double));

  for (i=0; i<nt; i++) {
    ha[i] = (double) 0.0;
  }

  retfoc   = (double) 0.0;
  retsteer = (double) 0.0;
  weight   = (double) 1.0;
  icheck   = 1;

  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, isize);

  for (i=0; i<isize; i++) {
    center_pos(&xs, &ys, &zs, i, gx, gy, gz);
    focusing(ifoc, focal[i], xs, ys, xamax, yamax, ramax, cp, &retfoc);  // Note ifoc must be 6 here!
    beamsteering(ister, theta, phi, xs, ys, xamax, yamax, ramax, cp, &retsteer);
    weighting(iweight, iapo, i, apod, &weight, xs, ys, ramax, param, isize);

    err = rect_ab(xo,yo,zo,xs,ys,zs,a,b,dx,dy,dt,nt,icheck,
                 delay,retfoc,retsteer,v,cp,alfa,weight,h,err_level);
    if (err != NONE)
      out_err = err;

    if (icheck == 2) {
      break;
    }
    superpos(h, ha, nt);
  }

  free(h);

  return out_err;
} /* dream_arr_rect_ud */


/********************************************************************/

/***
 *
 * rect_ab
 *
 * for calculatung the impulse respone of one rectangular element for use within an array.
 *
 ***/

int rect_ab(double xo, double yo, double zo, double xs, double ys, double zs, double a, double b,
            double dx, double dy, double dt, dream_idx_type nt, int icheck, double delay, double retfoc,
            double retsteer, double v, double cp, double alfa, double weight, double *RESTRICT h, int err_level)
{
  dream_idx_type i;
  double t, decal;
  double xsmin, ysmin, xsmax, ysmax, ai, ds, pi, ri;
  dream_idx_type it;
  double xsi, ysj;
  int err = NONE;

  decal = retfoc + retsteer;
  pi = atan( (double) 1.0) * 4.0;
  ds = dx * dy;

  xsmin = xs - a/2;
  xsmax = xs + a/2;

  ysmin = ys - b/2;
  ysmax = ys + b/2;

  for (i = 0; i < nt; i++)
    h[i] = (double) 0.0;

  ysj = ysmin + dy/2.0;
  while (ysj <= ysmax) {

    xsi = xsmin + dx/2.0;
    while (xsi <= xsmax) {

      distance(xo, yo, zo, xsi, ysj, zs, &ri);
      ai = weight * v * ds / (2*pi*ri);
      ai /= dt;
      ai *= 1000;      // Convert to SI units.

      t = ri * 1000/cp; // Propagation delay in micro seconds.
      it = (dream_idx_type) rint((t - delay + decal)/dt);

      if ((it < nt) && (it >= 0)) {

        if (alfa == (double) 0.0 )
          h[it] += ai; // No attenuation.
        else
          att(alfa, ri, it, dt, cp, h, nt, ai); // Attenuation compensation.
      }
      else   {
        if  (it >= 0)
          err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
        else
          err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

        if ( (err_level == PARALLEL_STOP) || (err_level == STOP) )
          return err; // Bail out.
      }

      xsi += dx;
    } // while
    ysj += dy;
  } // while

  return err;
} // rect_ab
