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
#include <stdlib.h>

#include "arr_functions.h"

/***
 *
 * Subroutines for common for array transducers.
 *
 *
 *
 ***/

/***
*
*  center_pos
*
*  Return the center point of an array element.
*
***/

void center_pos(double *RESTRICT xs, double *RESTRICT ys, double *RESTRICT zs, int i,
                double *RESTRICT gx, double *RESTRICT gy, double *RESTRICT gz)
{
  *xs = gx[i];
  *ys = gy[i];
  *zs = gz[i];

  return;
}

/***
 *
 *  max_dim_arr
 *
 *  Computes the maximum aperture of the array.
 *
 ***/

void max_dim_arr(double *RESTRICT x_max, double *RESTRICT y_max, double *RESTRICT ramax,
                 double *RESTRICT gx, double *RESTRICT gy, double *RESTRICT gz, int num_elements)
{
  int i;
  double ret;

  *x_max = fabs(gx[0]);
  for (i=1; i<num_elements; i++) {
    ret =  fabs(gx[i]);
    if (ret > *x_max) {
      *x_max = ret;
    }
  }

  *y_max = fabs(gy[0]);
  for (i=1; i<num_elements; i++) {
    ret = fabs(gy[i]);
    if (ret > *y_max) {
      *y_max = ret;
    }
  }
  *ramax = sqrt(*x_max * *x_max + *y_max * *y_max);

  return;
}

/***
 *
 * focusing
 *
 * Computes the focusing delay, retfoc, for point (xs,ys) on the
 * transducer surface (?).
 *
 ***/

void focusing(int foc_type, double focal, double xs, double ys,
              double x_max, double y_max, double ramax, double cp, double *RESTRICT retfoc)
{
  double diff, rmax, retx, rety;

  //
  // foc_type = 1 No foc, 2 Foc x ,3 Foc y, 4 Foc xy 5 Foc x+y, 6 ud (user defined) */
  //

  switch(foc_type) {

  case NO_FOCUS:
    return;

  case FOCUS_X:
    rmax = sqrt(x_max*x_max + focal*focal);
    diff = rmax - sqrt(xs*xs + focal*focal);
    *retfoc = diff*1000/cp;
    break;

  case FOCUS_Y:
    rmax = sqrt(y_max*y_max + focal*focal);
    diff  = rmax - sqrt(ys*ys + focal*focal);
    *retfoc = diff*1000/cp;
    break;

  case FOCUS_XY:
    rmax = sqrt(ramax*ramax + focal*focal);
    diff = rmax - sqrt(xs*xs + ys*ys + focal*focal);
    *retfoc = diff*1000/cp;
    break;

  case FOCUS_X_Y:
    rmax = sqrt(ramax*ramax + focal*focal);
    retx = sqrt(xs*xs + focal*focal);
    rety = sqrt(ys*ys + focal*focal);
    diff = rmax - (retx + rety);
    *retfoc = diff*1000/cp;
    break;

  case FOCUS_UD:
    *retfoc = focal; // Here focal is the user defined time delay in [us] (not the focal depth).
    break;

  default:
    break;

  }

  return;
} /* focusing */



/***
 *
 * Subroutine for beam steering
 *
 ***/

void beamsteering(int ister, double theta, double phi, double xs, double ys,
                  double x_max, double y_max, double ramax, double cp, double *RESTRICT retsteer)
{
  double diff, rmax, sinx, siny, retsteerx, retsteery;

  double pi = 4.0 * atan(1.0);
  double pii = pi / (double) 180.0;

  //
  // ister = 1 No steering, 2 Steer x ,3 Steer y, 4 Steer xy.
  //

  switch(ister) {

  case 0:
  case NO_STEER:
    break;

  case STEER_X:
    sinx = sin(theta * pii);
    rmax = x_max * sinx;
    diff =  rmax + xs*sinx;
    *retsteer = diff*1000/cp;
    break;

  case STEER_Y:
    siny = sin(phi * pii);
    rmax = y_max * siny;
    diff = rmax + ys * siny;
    *retsteer = diff*1000/cp;
    break;

  case STEER_XY:
    sinx = sin(theta * pii);
    rmax = x_max * sinx;
    diff = rmax + xs*sinx;
    retsteerx = diff*1000/cp;
    siny = sin(phi * pii);
    rmax = y_max * siny;
    diff = rmax + ys * siny;
    retsteery = diff*1000/cp;
    *retsteer = retsteerx + retsteery;
    break;

  default:
    break;
  }

  return;
} /* beamsteering */

/***
 *
 *  Aperture weighting/Apodization
 *
 * param=input parameter
 *
 * apod_type = 0 User defined
 * apod_type = 1 Traingle
 * apod_type = 2 Gauss
 * apod_type = 3 Rised cosine
 * apod_type = 4 Simply supported
 * apod_type = 5 Clamped
 *
 ***/

void apodization(int apod_type, int i, double *RESTRICT apod_vec, double *RESTRICT weight,
                 double xs, double ys, double ramax, double param)
{
  double pi = 4.0 * atan(1.0);
  double r = sqrt(xs*xs + ys*ys);

  switch(apod_type) {

  case APOD_UD:
    *weight = apod_vec[i];
    break;

  case APOD_TRIANGLE:
    *weight = 1.0 - fabs(r) / ramax;
    break;

  case APOD_GAUSS:
    *weight = exp(-(param * r*r) / (ramax*ramax));
    break;

  case APOD_RISED_COSINE:
    *weight = param + cos(r*pi/ramax);
    break;

  case APOD_SIMPLY_SUPPORTED:
    *weight = 1.0 - r*r / (ramax*ramax);
    break;

  case APOD_CLAMPED:
    *weight = (1.0 - r*r / (ramax*ramax)) * (1.0  - r*r / (ramax*ramax));

  default:
    break;
  }

  return;
}

/***
 *
 * distance
 *
 * Compute the distance bentween two points (vector length/magnitude). Normally used here
 * to compute the distance (length) from an observation point (xo,yo,zo)
 * to a point (xs,ys,zs) on the transducer surface.
 *
 ***/

void distance(double xo, double yo, double zo,double xs,double ys, double zs, double *RESTRICT ri)
{
  double rx, ry, rz;

  rx = xo - xs;
  ry = yo - ys;
  rz = zo - zs;
  *ri = sqrt(rx*rx + rz*rz + ry*ry);

  return;
}

/***
 *
 * superpos
 *
 * Add contributions from an array element.
 *
 * ha = output response
 * h  = input responce of actual element
 *
 ***/

void superpos(double *RESTRICT h, double *RESTRICT ha, dream_idx_type nt)
{
  dream_idx_type i;

  for (i=0; i< nt; i++)
    ha[i] += h[i];

  return;
}

/***
 *
 * check_delay
 *
 * Check that the delay is within bounds of the impulse response vector.
 *
 ***/

bool check_delay(dream_idx_type it, double t, dream_idx_type nt)
{
  bool retval=true;

  if (t < 0.0) {
    retval = false;
  }

  if (it > nt) {
    retval = false;
  }

  return retval;
}
