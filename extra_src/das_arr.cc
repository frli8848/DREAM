/***
*
* Copyright (C) 2003,2006,2006,2007,2008,2009,2014,2021 Fredrik Lingvall
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

//
// Function prototypes.
//

int delay_arr(double xo, double yo, double zo, double xs, double ys, double zs, double dt,
              dream_idx_type nt, double delay, double foc_delay, double steer_delay, double cp,
              double weight, double *h, int err_level);

int centroid(double *h,dream_idx_type nt);

/***
 *
 * das_arr - Delay-and-sum (DAS) array processing.
 *
 ***/

int das_arr(double xo, double yo, double zo, double dt, dream_idx_type nt,
            double delay, double cp, int  num_elements,
            double *gx, double *gy, double *gz, int ifoc, double focal,
            int ister, double theta, double phi, double *apod, bool do_apod,
            int apod_type, double param, double *ha,int err_level)
{
  double steer_delay;
  double *h;
  dream_idx_type i, i_c;
  double ramax, xamax, yamax;
  double xs, ys, zs, foc_delay, weight;
  int err = NONE, out_err = NONE;

  h = (double*) malloc(nt*sizeof(double));

  for (i=0; i<nt; i++) {
    ha[i] = (double) 0.0;
  }

  foc_delay   = (double) 0.0;
  steer_delay = (double) 0.0;
  weight   = (double) 1.0;

  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, num_elements);

  for (i=0; i<num_elements; i++) {
    center_pos(&xs, &ys, &zs, i, gx, gy, gz);
    focusing(ifoc, focal, xs, ys, xamax, yamax, ramax, cp, &foc_delay);
    beamsteering(ister, theta, phi, xs, ys, xamax, yamax, ramax, cp, &steer_delay);
    if (do_apod) {
      apodization(apod_type, i, apod, &weight, xs, ys, ramax, param);
    }
    err = delay_arr(xo,yo,zo,xs,ys,zs,dt,nt,delay,foc_delay,steer_delay,cp,weight,h,err_level);
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
                    double delay, double cp, int  num_elements,
                      double *gx, double *gy, double *gz, int ifoc, double *focal,
                      int ister, double theta, double phi, double *apod, bool do_apod,
                      int apod_type, double param, double *ha,int err_level)
{
  double steer_delay;
  double *h;
  dream_idx_type i, i_c;
  double ramax, xamax, yamax;
  double xs, ys, zs, foc_delay, weight;
  int err = NONE, out_err = NONE;

  h = (double*) malloc(nt*sizeof(double));

  for (i=0; i<nt; i++) {
    ha[i] = (double) 0.0;
  }

  foc_delay   = (double) 0.0;
  steer_delay = (double) 0.0;
  weight   = (double) 1.0;

  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, num_elements);

  for (i=0; i<num_elements; i++) {
    center_pos(&xs, &ys, &zs, i, gx, gy, gz);
    focusing(ifoc, focal[i], xs, ys, xamax, yamax, ramax, cp, &foc_delay);   // Note ifoc must be 6 here!
    beamsteering(ister, theta, phi, xs, ys, xamax, yamax, ramax, cp, &steer_delay);
    if (do_apod) {
      apodization(apod_type, i, apod, &weight, xs, ys, ramax, param);
    }
    err = delay_arr(xo,yo,zo,xs,ys,zs,dt,nt,delay,foc_delay,steer_delay,cp,weight,h,err_level);
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
              dream_idx_type nt, double delay, double foc_delay, double steer_delay, double cp,
              double weight, double *h, int err_level)
{
  int err = NONE;

  for (dream_idx_type it = 0; it < nt; it++) {
    h[it] = 0.0;
  }

  double ri;
  distance(xo, yo, zo, xs, ys, zs, &ri);
  double t = ri * 1.0e3/cp;
  dream_idx_type it = (dream_idx_type) rint((t - delay + foc_delay + steer_delay)/dt);

  // Check if index is out of bounds.
  if ((it < nt) && (it >= 0)) {
    h[it] += 1.0;
  } else  {
    if  (it >= 0) {
      err = dream_out_of_bounds_err("Delay out of bounds",it-nt+1,err_level);
    } else {
      err = dream_out_of_bounds_err("Delay out of bounds",it,err_level);
    }

    if ( (err == PARALLEL_STOP) || (err == STOP) ) {
      return err; // Bail out.
    }
  }

  return err;
} /* delay_arr */

/***
 *
 * Compute center-of-mass.
 *
 ***/

int centroid(double *h, dream_idx_type nt)
{
  double sum=0.0, t_c = 0.0;

  for (dream_idx_type n=0; n<nt; n++) {
    sum += h[n];
    t_c += ((double) n) * h[n];
  }

  t_c /= sum;

  return rint(t_c);
}
