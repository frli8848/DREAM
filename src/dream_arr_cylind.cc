/***
*
* Copyright (C) 2002,2003,2005,2006,2007,2008,2009,2014,2019,2021 Fredrik Lingvall
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

#include "dream_arr_cylind.h"
#include "dreamcylind.h"
#include "dream_error.h"
#include "attenuation.h"
#include "arr_functions.h"

/***
 *
 *  Array with cylindrical concave or convex elements.
 *
 ***/

int dream_arr_cylind(double xo, double yo, double zo,
                     double a, double b, double Rcurv,
                     double dx, double dy, double dt, dream_idx_type nt,
                     double delay, double v, double cp,
                     dream_idx_type num_elements, double *gx, double *gy, double *gz,
                     FocusMet foc_met, double *focal,
                     SteerMet steer_met, double theta, double phi,
                     double *apod, bool do_apod, ApodMet apod_met, double param,
                     double *ha, int err_level)
{
  dream_idx_type i;
  double ramax, xamax, yamax;
  //double xs, ys, zs;
  int err = NONE, out_err = NONE;

  for (i=0; i<nt; i++) {
    ha[i] = 0.0;
  }

  double foc_delay   = 0.0;
  double steer_delay = 0.0;
  double weight = 1.0;

  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, num_elements);

  for (i=0; i<num_elements; i++) {

    //center_pos(&xs, &ys, &zs, i, gx, gy, gz);
    //focusing(foc_met, focal, xs, ys, xamax, yamax, ramax, cp, &retfoc);

    if (foc_met != FocusMet::ud) {
      focusing(foc_met, focal[0], gx[i], gy[i], xamax, yamax, ramax, cp, &foc_delay);
    } else {
      focusing(foc_met, focal[i], gx[i], gy[i], xamax, yamax, ramax, cp, &foc_delay);
    }

    //beamsteering(steer_met, theta, phi, xs, ys, xamax, yamax, ramax, cp, &retsteer);
    beamsteering(steer_met, theta, phi, gx[i], gy[i], xamax, yamax, ramax, cp, &steer_delay);

    if (do_apod) {
      //apodization(apod_met, i, apod, &weight, xs, ys, ramax, param);
      apodization(apod_met, i, apod, &weight, gx[i], gy[i], ramax, param);
    }
    // Compute the response for the i:th element and add it to the impulse response vector ha.
    //err = cylind_f(xo,yo,zo,xs,ys,zs,Rcurv,a,b,dx,dy,dt,nt,delay,retfoc,retsteer,v,cp,alpha,weight,ha,err_level);
    err = dreamcylind(gx[i], gy[i], 0.0,
                      a, b, Rcurv,
                      dx, dy, dt, nt, delay - foc_delay - steer_delay,
                      v, cp,
                      ha, err_level,
                      weight);

    if (err != NONE) {
      out_err = err;
    }
  }

  return out_err;
}

int dream_arr_cylind(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                     double xo, double yo, double zo,
                     double a, double b, double Rcurv,
                     double dx, double dy, double dt, dream_idx_type nt,
                     double delay, double v, double cp,
                     dream_idx_type num_elements, double *gx, double *gy, double *gz,
                     FocusMet foc_met, double *focal,
                     SteerMet steer_met, double theta, double phi,
                     double *apod, bool do_apod, ApodMet apod_met, double param,
                     double *ha, int err_level)
{
  dream_idx_type i;
  double ramax, xamax, yamax;
  double xs, ys, zs;
  int err = NONE, out_err = NONE;

  for (i=0; i<nt; i++) {
    ha[i] = 0.0;
  }

  double foc_delay   = 0.0;
  double steer_delay = 0.0;
  double weight = 1.0;

  max_dim_arr(&xamax, &yamax, &ramax, gx, gy, gz, num_elements);

  for (i=0; i<num_elements; i++) {

    //center_pos(&xs, &ys, &zs, i, gx, gy, gz);
    //focusing(foc_met, focal, xs, ys, xamax, yamax, ramax, cp, &retfoc);

    if (foc_met != FocusMet::ud) {
      focusing(foc_met, focal[0], gx[i], gy[i], xamax, yamax, ramax, cp, &foc_delay);
    } else {
      focusing(foc_met, focal[i], gx[i], gy[i], xamax, yamax, ramax, cp, &foc_delay);
    }

    //beamsteering(steer_met, theta, phi, xs, ys, xamax, yamax, ramax, cp, &retsteer);
    beamsteering(steer_met, theta, phi, gx[i], gy[i], xamax, yamax, ramax, cp, &steer_delay);

    if (do_apod) {
      //apodization(apod_met, i, apod, &weight, xs, ys, ramax, param);
      apodization(apod_met, i, apod, &weight, gx[i], gy[i], ramax, param);
    }
    // Compute the response for the i:th element and add it to the impulse response vector ha.
    //err = cylind_f(xo,yo,zo,xs,ys,zs,Rcurv,a,b,dx,dy,dt,nt,delay,retfoc,retsteer,v,cp,alpha,weight,ha,err_level);
    err = dreamcylind(att, xc_vec, x_vec,
                      gx[i], gy[i], 0.0,
                      a, b, Rcurv,
                      dx, dy, dt, nt, delay - foc_delay - steer_delay,
                      v, cp,
                      ha, err_level,
                      weight);

    if (err != NONE) {
      out_err = err;
    }
  }

  return out_err;
}
