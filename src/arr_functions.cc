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

#include <cmath>

#include "arr_functions.h"

/***
 *
 * Subroutines for common for array transducers.
 *
 ***/

/***
*
*  center_pos
*
*  Return the center point of an array element.
*
***/

void center_pos(double *x, double *y, double *z, dream_idx_type i,
                double *gx, double *gy, double *gz)
{
  *x = gx[i];
  *y = gy[i];
  *z = gz[i];

  return;
}

/***
 *
 *  max_dim_arr
 *
 *  Computes the maximum aperture of the array.
 *
 ***/

void max_dim_arr(double *x_max, double *y_max, double *r_max,
                 double *gx, double *gy, double *gz, dream_idx_type num_elements)
{
  dream_idx_type i;
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
  *r_max = std::sqrt(*x_max * *x_max + *y_max * *y_max);

  return;
}

/***
 *
 * focusing
 *
 * Computes the focusing delay, foc_delay, for array element at (gx,gy).
 *
 ***/

void focusing(FocusMet foc_met, double focal, double gx, double gy,
              double x_max, double y_max, double r_max, double cp, double *foc_delay)
{
  switch(foc_met) {

  case FocusMet::none:
    *foc_delay = 0.0;
    return;

  case FocusMet::x:
    {
      double rmax = std::sqrt(x_max*x_max + focal*focal);
      double diff = rmax - std::sqrt(gx*gx + focal*focal);
      *foc_delay = diff*1.0e3/cp; // [us]
    }
    break;

  case FocusMet::y:
    {
      double rmax = std::sqrt(y_max*y_max + focal*focal);
      double diff = rmax - std::sqrt(gy*gy + focal*focal);
      *foc_delay = diff*1.0e3/cp; // [us]
    }
    break;

  case FocusMet::xy:
    {
      double rmax = std::sqrt(r_max*r_max + focal*focal);
      double diff = rmax - std::sqrt(gx*gx + gy*gy + focal*focal);
      *foc_delay = diff*1.0e3/cp; // [us]
    }
    break;

  case FocusMet::x_y:
    {
      double rmax = std::sqrt(r_max*r_max + focal*focal);
      double retx = std::sqrt(gx*gx + focal*focal);
      double rety = std::sqrt(gy*gy + focal*focal);
      double diff = rmax - (retx + rety);
      *foc_delay = diff*1.0e3/cp; // [us]
    }
    break;

  case FocusMet::ud:
    *foc_delay = focal; // Here focal is the user defined time delay in [us] (not the focal depth).
    break;

  default:
    *foc_delay = 0.0;
    break;

  }

  return;
}

/***
 *
 * Beam steering
 *
 ***/

void beamsteering(SteerMet steer_met, double theta, double phi, double gx, double gy,
                  double x_max, double y_max, double r_max, double cp, double *steer_delay)
{
  const double deg2rad = M_PI / (double) 180.0;

  switch(steer_met) {

  case SteerMet::none:
    *steer_delay = 0.0;
    break;

  case SteerMet::x:
    {
      double sinx = std::sin(theta * deg2rad);
      double rmax = x_max * sinx;
      double diff =  rmax + gx*sinx;
      *steer_delay = diff*1.0e3/cp;
    }
    break;

  case SteerMet::y:
    {
      double siny = std::sin(phi * deg2rad);
      double rmax = y_max * siny;
      double diff = rmax + gy * siny;
      *steer_delay = diff*1.0e3/cp;
    }
    break;

  case SteerMet::xy:
    {
      double sinx = std::sin(theta * deg2rad);
      double rmax = x_max * sinx;
      double diff = rmax + gx*sinx;
      double steer_x_delay = diff*1.0e3/cp;
      double siny = std::sin(phi * deg2rad);
      rmax = y_max * siny;
      diff = rmax + gy * siny;
      double steer_y_delay = diff*1.0e3/cp;
      *steer_delay = steer_x_delay + steer_y_delay;
    }
    break;

  default:
    *steer_delay = 0.0;
    break;
  }

  return;
}

/***
 *
 * Aperture weighting/Apodization
 *
 * apod_par: apodization shape parameter.
 *
 ***/

void apodization(ApodMet apod_met, dream_idx_type n, double *apod_vec, double *weight,
                 double gx, double gy, double r_max, double apod_par)
{
  double r = std::sqrt(gx*gx + gy*gy);

  switch(apod_met) {

  case ApodMet::ud:
    *weight = apod_vec[n];
    break;

  case ApodMet::triangle:
    *weight = 1.0 - fabs(r) / r_max;
    break;

  case ApodMet::gauss:
    *weight = exp(-(apod_par * r*r) / (r_max*r_max));
    break;

  case ApodMet::raised_cosine:
    *weight = apod_par + std::cos(r*M_PI/r_max);
    break;

  case ApodMet::hann:
    *weight = 0.5 * (1.0 - std::cos((r+r_max)*M_PI/r_max));
    break;

  case ApodMet::hamming:
    {
      double a0 = 25.0/46.0;
      *weight = a0 - (1.0 - a0) * std::cos((r+r_max)*M_PI/r_max);
    }
    break;

  case ApodMet::simply_supported:
    *weight = 1.0 - r*r / (r_max*r_max);
    break;

  case ApodMet::clamped:
    *weight = (1.0 - r*r / (r_max*r_max)) * (1.0  - r*r / (r_max*r_max));
    break;

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

void distance(double xo, double yo, double zo,
              double xs, double ys, double zs,
              double *ri)
{
  double rx, ry, rz;

  rx = xo - xs;
  ry = yo - ys;
  rz = zo - zs;
  *ri = std::sqrt(rx*rx + rz*rz + ry*ry);

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

void superpos(double *h, double *ha, dream_idx_type nt)
{
  dream_idx_type i;

  for (i=0; i< nt; i++) {
    ha[i] += h[i];
  }

  return;
}
