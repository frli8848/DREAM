/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2014,2019,2021 Fredrik Lingvall
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

//
// Function prototypes.
//

int circ_arr(double xo, double yo, double zo, double xs, double ys, double r, double dx, double dy, double dt,
              dream_idx_type nt, double delay, double retfoc, double retsteer, double v, double cp, double alpha,
              double weight, double *h, int err_level);
void xlimit_arr_circ(double y, double r, double xs, double ys, double *x_min, double *x_max);

/***
 *
 * dream_arr_cir - 2D array with circular elements.
 *
 ***/

int dream_arr_circ(double xo, double yo, double zo, double r, double dx, double dy, double dt, dream_idx_type nt,
                    double delay, double v, double cp, double alpha, int  num_elements,
                    double *gx, double *gy, double *gz, int foc_type, double focal,
                    int ister, double theta, double phi, double *apod, bool do_apod,
                    int apod_type, double param, double *ha, int err_level)
{
  double retsteer;
  dream_idx_type i;
  double ramax, xamax, yamax;
  double xs, ys, zs, retfoc, weight;
  int err = NONE, out_err = NONE;

  for (i=0; i<nt; i++) {
    ha[i] = 0.0;
  }

  retfoc   = 0.0;
  retsteer = 0.0;
  weight   = 1.0;

  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, num_elements);

  for (i=0; i<num_elements; i++) {
    center_pos(&xs, &ys, &zs, i, gx, gy, gz);
    focusing(foc_type, focal, xs, ys, xamax, yamax, ramax, cp, &retfoc);
    beamsteering(ister, theta, phi, xs, ys, xamax, yamax, ramax, cp, &retsteer);
    if (do_apod) {
      apodization(apod_type, i, apod, &weight, xs, ys, ramax, param);
    }
    // Compute the response for the i:th element and add it to the impulse response vector ha.
    err = circ_arr(xo,yo,zo,xs,ys,r,dx,dy,dt,nt,delay,retfoc,retsteer,v,cp,alpha,weight,ha,err_level);
    if (err != NONE)
      out_err = err;
  }

  return out_err;
}

/***
 *
 * dream_arr_circ_ud - 2D array with circular elements - user defined focusing.
 *
 ***/

int dream_arr_circ_ud(double xo, double yo, double zo, double r, double dx, double dy, double dt, dream_idx_type nt,
                    double delay, double v, double cp, double alpha, int  num_elements,
                    double *gx, double *gy, double *gz, int foc_type, double *focal,
                    int ister, double theta, double phi, double *apod, bool do_apod,
                    int apod_type, double param, double *ha, int err_level)
{
  double retsteer;
  dream_idx_type i;
  double ramax, xamax, yamax;
  double xs, ys, zs, retfoc, weight;
  int err = NONE, out_err = NONE;

  for (i=0; i<nt; i++) {
    ha[i] = (double) 0.0;
  }

  retfoc   = (double) 0.0;
  retsteer = (double) 0.0;
  weight   = (double) 1.0;

  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, num_elements);

  for (i=0; i<num_elements; i++) {
    center_pos(&xs, &ys, &zs, i, gx, gy, gz);
    focusing(foc_type, focal[i], xs, ys, xamax, yamax, ramax, cp, &retfoc);   // Note foc_type must be 6 here!
    beamsteering(ister, theta, phi, xs, ys, xamax, yamax, ramax, cp, &retsteer);
    if (do_apod) {
      apodization(apod_type, i, apod, &weight, xs, ys, ramax, param);
    }
    // Compute the response for the i:th element and add it to the impulse response vector ha.
    err = circ_arr(xo,yo,zo,xs,ys,r,dx,dy,dt,nt,delay,retfoc,retsteer,v,cp,alpha,weight,ha,err_level);
    if (err != NONE)
      out_err = err;
  }

  return out_err;
}


/***
 *
 * circ_arr
 *
 * NB. We add (super impose) the response to impulse response vector h!
 *
 ***/

int circ_arr(double xo, double yo, double zo, double xs, double ys, double r, double dx, double dy, double dt,
              dream_idx_type nt, double delay, double retfoc, double retsteer, double v, double cp, double alpha,
              double weight, double *h, int err_level)
{
  dream_idx_type i;
  double t, decal;
  double x_min, ysmin, x_max, ysmax, ai, ds, pi, ri;
  dream_idx_type    it;
  double zs, x, y;
  int err = NONE;

  pi = 4.0 * atan(1.0);
  decal = retfoc + retsteer;
  ds = dx * dy;
  zs = (double) 0.0;
  ysmin = -r + ys;
  ysmax =  r + ys;

  y = ysmin + dy/2.0;
  while (y <= ysmax) {

    xlimit_arr_circ(y, r, xs, ys, &x_min, &x_max);

    x = x_min + dx/2.0;
    while (x <= x_max) {

      //distance(xo, yo, zo, x, y, zs, &ri, &rx, &ry, &rz);
      distance(xo, yo, zo, x, y, zs, &ri);
      ai = weight * v * ds / (2*pi*ri);
      ai /= dt;
      ai *= 1000; // Convert to SI units.

      // Propagation delay in micro seconds.
      t = ri * 1000/cp;
      it = (dream_idx_type) rint((t - delay + decal)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {

        // Check if absorbtion is present.
        if (alpha == 0.0) {
          h[it] += ai;
        } else {
          att(alpha,ri,it,dt,cp,h,nt,ai);
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
    y += dy;
  }

  return err;
}

/***
 *
 * xlimit_arr_circ
 *
 * Computes the x-axis integration limits.
 *
 ***/

void xlimit_arr_circ(double y, double r, double xs, double ys, double *x_min, double *x_max)
{
  double rs;

  rs = r*r - (ys-y)*(ys-y);
  rs = sqrt(rs);
  *x_min = -rs + xs;
  *x_max =  rs + xs;

  return;
}
