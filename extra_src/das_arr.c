/***
*
* Copyright (C) 2003,2006,2006,2007,2008,2009,2014 Fredrik Lingvall
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

#include "das_arr.h"
#include "arr_functions.h"
#include "dream_error.h"

#if defined(_MSC_VER) || defined(__LCC__)
#include "msvc_rint.h"
#endif

//
// Function prototypes.
//

int delay_arr(double xo, double yo, double zo, double xs, double ys, double zs, double dt,
              dream_idx_type nt, double delay, double retfoc, double retsteer, double cp,
              double weight, double *RESTRICT h, int err_level);

int centroid(double *RESTRICT h,dream_idx_type nt);

/***
 *
 * das_arr - Delay-and-sum (DAS) array processing.
 *
 ***/

int das_arr(double xo, double yo, double zo, double dt, dream_idx_type nt,
            double delay, double cp, int  isize,
            double *RESTRICT gx, double *RESTRICT gy, double *RESTRICT gz, int ifoc, double focal,
            int ister, double theta, double phi, double *RESTRICT apod, int iweight,
            int iapo, double param, double *RESTRICT ha,int err_level)
{
  double retsteer;
  double *RESTRICT h;
  dream_idx_type i, i_c;
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
    focusing(ifoc, focal, xs, ys, xamax, yamax, ramax, cp, &retfoc);
    beamsteering(ister, theta, phi, xs, ys, xamax, yamax, ramax, cp, &retsteer);
    apodization(iweight, iapo, i, apod, &weight, xs, ys, ramax, param, isize);

    err = delay_arr(xo,yo,zo,xs,ys,zs,dt,nt,delay,retfoc,retsteer,cp,weight,h,err_level);
    if (err != NONE) {
      out_err = err;
      if ( (err == PARALLEL_STOP) || (err == STOP) ) {
        return err; // Bail out.
      }
    }

    superpos(h, ha, nt);
  }

  i_c = centroid(ha,nt);
  for (i = 0; i < nt; i++) {
    ha[i] = (double) 0.0;
  }

  if ((i_c < nt) && (i_c >= 0)) {
    ha[i_c] += 1.0;
  }
  else {
    if  (i_c >= 0)
      err = dream_out_of_bounds_err("Centroid out of bounds",i_c-nt+1,err_level);
    else
      err = dream_out_of_bounds_err("Centroid out of bounds",i_c,err_level);

    if ( (err == PARALLEL_STOP) || (err == STOP) ) {
      free(h);
      return err; // Bail out.
    }

  }

  free(h);

  return out_err;
} /* dream_arr_circ */

/***
 *
 * das_arr_ud - Delay-and-sum array processing - user defined focusing.
 *
 ***/

int das_arr_ud(double xo, double yo, double zo, double dt, dream_idx_type nt,
                    double delay, double cp, int  isize,
                      double *RESTRICT gx, double *RESTRICT gy, double *RESTRICT gz, int ifoc, double *RESTRICT focal,
                      int ister, double theta, double phi, double *RESTRICT apod, int iweight,
                      int iapo, double param, double *RESTRICT ha,int err_level)
{
  double retsteer;
  double *RESTRICT h;
  dream_idx_type i, i_c;
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

    err = delay_arr(xo,yo,zo,xs,ys,zs,dt,nt,delay,retfoc,retsteer,cp,weight,h,err_level);
    if (err != NONE) {
      out_err = err;
      if ( (err == PARALLEL_STOP) || (err == STOP) ) {
        return err; // Bail out.
      }
    }

    superpos(h, ha, nt);
  }


  i_c = centroid(ha,nt);

  for (i = 0; i < nt; i++) {
    ha[i] = (double) 0.0;
  }

  if ((i_c < nt) && (i_c >= 0)) {
    ha[i_c] += 1.0;
  }
  else {
    if  (i_c >= 0)
      err = dream_out_of_bounds_err("Centroid out of bounds",i_c-nt+1,err_level);
    else
      err = dream_out_of_bounds_err("Centroid out of bounds",i_c,err_level);

    if ( (err == PARALLEL_STOP) || (err == STOP) ) {
      free(h);
      return err; // Bail out.
    }

  }

  free(h);

  return out_err;
} /* das_arr_ud */

/********************************************************************/


/***
 *
 * Compute the time-of-flight (delay) from array element to observation point.
 *
 ***/

int delay_arr(double xo, double yo, double zo, double xs, double ys, double zs, double dt,
              dream_idx_type nt, double delay, double retfoc, double retsteer, double cp,
              double weight, double *RESTRICT h, int err_level)
{
  dream_idx_type    i;
  double t,decal,tt;
  double ri;
  dream_idx_type    it;
  double qan;
  int err = NONE;

  decal = retfoc + retsteer;

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0;
  }

  distance(xo, yo, zo, xs, ys, zs, &ri);
  t = ri * 1000/cp;
  tt = t - delay + decal;
  qan = tt / dt;
  it = (dream_idx_type) rint(qan);

  // Check if index is out of bounds.
  if ((it < nt) && (it >= 0)) {
    h[it] += 1.0;
  }
  else  {
    if  (it >= 0)
      err = dream_out_of_bounds_err("Delay out of bounds",it-nt+1,err_level);
    else
      err = dream_out_of_bounds_err("Delay out of bounds",it,err_level);

    if ( (err == PARALLEL_STOP) || (err == STOP) )
      return err; // Bail out.
  }

  return err;
} /* delay_arr */

/***
 *
 * Compute center-of-mass.
 *
 ***/

int centroid(double *RESTRICT h, dream_idx_type nt)
{
  dream_idx_type n;
  double sum=0.0, t_c = 0.0;

  for (n=0; n<nt; n++) {
    sum += h[n];
    t_c += ((double) n) * h[n];
  }

  t_c /= sum;

  return rint(t_c);
}
