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

#include "dream_arr_circ.h"
#include "att.h"
#include "arr_functions.h"
#include "dream_error.h"

#if defined(_MSC_VER) || defined(__LCC__)
#include "msvc_rint.h"
#endif

//
// Function prototypes.
//

int circ_arr(double xo, double yo, double zo, double xs, double ys, double r, double dx, double dy, double dt,
              dream_idx_type nt, double delay, double retfoc, double retsteer, double v, double cp, double alpha,
              double weight, double *RESTRICT h, int err_level);
void xlimit(double yi, double r, double xs, double ys, double *RESTRICT xsmin, double *RESTRICT xsmax);

/***
 *
 * dream_arr_cir - 2D array with circular elements.
 *
 ***/

int dream_arr_circ(double xo, double yo, double zo, double r, double dx, double dy, double dt, dream_idx_type nt,
                    double delay, double v, double cp, double alpha, int  isize,
                    double *RESTRICT gx, double *RESTRICT gy, double *RESTRICT gz, int ifoc, double focal,
                    int ister, double theta, double phi, double *RESTRICT apod, int iweight,
                    int iapo, double param, double *RESTRICT ha, int err_level)
{
  double retsteer;
  double *RESTRICT h;
  dream_idx_type i;
  double ramax, xamax, yamax;
  double xs, ys, zs, retfoc, weight;
  int err = NONE, out_err = NONE;

  h = (double*) malloc(nt*sizeof(double));

  for (i=0; i<nt; i++)
    ha[i] = (double) 0.0;

  retfoc   = (double) 0.0;
  retsteer = (double) 0.0;
  weight   = (double) 1.0;

  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, isize);

  for (i=0; i<isize; i++) {
    center_pos(&xs, &ys, &zs, i, gx, gy, gz);
    focusing(ifoc, focal, xs, ys, xamax, yamax, ramax, cp, &retfoc);
    beamsteering(ister, theta, phi, xs, ys, xamax, yamax, ramax, cp, &retsteer);
    apodization(iweight, iapo, i, apod, &weight, xs, ys, ramax, param, isize);

    err = circ_arr(xo,yo,zo,xs,ys,r,dx,dy,dt,nt,delay,retfoc,retsteer,v,cp,alpha,weight,h,err_level);
    if (err != NONE)
      out_err = err;

    superpos(h, ha, nt);
  }

  free(h);

  return out_err;
} /* dream_arr_circ */

/***
 *
 * dream_arr_circ_ud - 2D array with circular elements - user defined focusing.
 *
 ***/

int dream_arr_circ_ud(double xo, double yo, double zo, double r, double dx, double dy, double dt, dream_idx_type nt,
                    double delay, double v, double cp, double alpha, int  isize,
                    double *RESTRICT gx, double *RESTRICT gy, double *RESTRICT gz, int ifoc, double *RESTRICT focal,
                    int ister, double theta, double phi, double *RESTRICT apod, int iweight,
                    int iapo, double param, double *RESTRICT ha, int err_level)
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

  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, isize);

  for (i=0; i<isize; i++) {
    center_pos(&xs, &ys, &zs, i, gx, gy, gz);
    focusing(ifoc, focal[i], xs, ys, xamax, yamax, ramax, cp, &retfoc);   // Note ifoc must be 6 here!
    beamsteering(ister, theta, phi, xs, ys, xamax, yamax, ramax, cp, &retsteer);
    apodization(iweight, iapo, i, apod, &weight, xs, ys, ramax, param, isize);

    err = circ_arr(xo,yo,zo,xs,ys,r,dx,dy,dt,nt,delay,retfoc,retsteer,v,cp,alpha,weight,h,err_level);
    if (err != NONE)
      out_err = err;

    superpos(h, ha, nt);
  }

  free(h);

  return out_err;
} /* dream_arr_circ_ud */

/********************************************************************/


/***
 *
 * circ_arr
 *
 ***/

int circ_arr(double xo, double yo, double zo, double xs, double ys, double r, double dx, double dy, double dt,
              dream_idx_type nt, double delay, double retfoc, double retsteer, double v, double cp, double alpha,
              double weight, double *RESTRICT h, int err_level)
{
  dream_idx_type    i;
  double t, decal;
  double xsmin, ysmin, xsmax, ysmax, ai, ds, pi, ri;
  dream_idx_type    it;
  double zs, xsi, ysj;
  int err = NONE;

  pi = atan((double) 1.0) * 4.0;
  decal = retfoc + retsteer;
  ds = dx * dy;
  zs = (double) 0.0;
  ysmin = -r + ys;
  ysmax =  r + ys;

  for (i = 0; i < nt; i++)
    h[i] = (double) 0.0 ;

  //j = 0;
  //j++;
  //ysj = ysmin + (j-1) * dy + dy/2;
  ysj = ysmin + dy/2.0;
  while (ysj <= ysmax) {

    xlimit(ysj, r, xs, ys, &xsmin, &xsmax);

    //i = 0;
    //i++;
    //xsi = xsmin + (i-1) * dx + dx/2;
    xsi = xsmin + dx/2.0;
    while (xsi <= xsmax) {

      //distance(xo, yo, zo, xsi, ysj, zs, &ri, &rx, &ry, &rz);
      distance(xo, yo, zo, xsi, ysj, zs, &ri);
      ai = weight * v * ds / (2*pi*ri);
      ai /= dt;
      // Convert to SI units.
      ai *= 1000;
      // Propagation delay in micro seconds.
      t = ri * 1000/cp;
      it = (dream_idx_type) rint((t - delay + decal)/dt);

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
    } // while
    //j++;
    //ysj = ysmin + (j-1)*dy + dy/2;
    ysj += dy;
  } // while

  return err;
} /* circ_arr */

/***
 *
 * subroutine xlimit - pour definir les limits d integration en x
 *
 ***/

void xlimit(double yi, double r, double xs, double ys, double *RESTRICT xsmin, double *RESTRICT xsmax)
{
  double rs;

  rs = r*r - (ys-yi)*(ys-yi);
  rs = sqrt(rs);
  *xsmin = -rs + xs;
  *xsmax =  rs + xs;

  return;
} /* xlimit */