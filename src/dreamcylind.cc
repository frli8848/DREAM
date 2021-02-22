/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2014,2021 Fredrik Lingvall
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

#include "dreamcylind.h"
#include "dream_error.h"

//
//  Function prototypes.
//

void cylind_f(double x, double y,
              double Rcurv,
              double z_Rcurv,
              double xo, double yo, double zo,
              double &r, double &du);

void cylind_d(double x, double y,
              double Rcurv, double z_Rcurv,
              double xo, double yo, double zo,
              double &r, double &du);

/***
 *
 * dreamcylind : cylindric concave or convex transducer
 *
 ***/

int dreamcylind(double xo, double yo, double zo,
                double a, double b, double Rcurv,
                double dx, double dy, double dt,
                dream_idx_type nt, double delay, double v, double cp,
                double *h,
                int err_level)
{
  dream_idx_type i, it;
  double t, xsmin, xsmax, ai, du, r;
  double phi, phi_min, phi_max;
  int err = NONE;

  if (b > 2*fabs(Rcurv)) {
    dream_err_msg("Error in dreamcylind: the y-size, b, must be less than the curvature diameter 2*Rcurv!\n");
  }

  b /= 2.0;
  double z_Rcurv = fabs(Rcurv) - sqrt(Rcurv*Rcurv - b*b);

  xsmin = -a/2.0;
  xsmax = a/2.0;

  phi = 2.0 * atan(b/(fabs(Rcurv)-z_Rcurv));
  phi_min = -phi/2;
  phi_max = phi/2;

  // dphi in y-dim [rad].
  double dphi = asin(dy/fabs(Rcurv));
  //dphi = dy/r;
  double ds = fabs(Rcurv) * dx * dphi;
  //ds = dx * dy; // Approx the same as  ds = Rcurv * dx * dphi.

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  phi = phi_min + dphi/2.0;
  double ys = fabs(Rcurv) * sin(phi);
  while (phi <= phi_max) {

    double xs = xsmin + dx/2.0;
    while (xs <= xsmax) {

      // Compute r and ds.
      if (Rcurv >= 0.0) {       // Focused
        cylind_f(xs, ys, Rcurv, z_Rcurv, xo, yo, zo, r, du);
      } else {                  // Defocused
        cylind_d(xs, ys, fabs(Rcurv), z_Rcurv, xo, yo, zo, r, du);
      }

      ai = v * ds * du/(2.0*M_PI * r);
      ai /= dt;
      // Convert to SI units.
      ai *= 1.0e3;
      // Propagation delay in micro seconds.
      t = r * 1.0e3/cp;
      it = (dream_idx_type) rint((t - delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {
        h[it] += ai;
      }
      else  {
        if  (it >= 0)
          err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
        else
          err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

        if ( (err_level == PARALLEL_STOP) || (err_level == STOP) ) {
          return err; // Bail out.
        }
      }

      xs += dx;
    }

    phi += dphi;
    ys = fabs(Rcurv) * sin(phi);
  }

  return err;
}

int dreamcylind(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                double xo, double yo, double zo,
                double a, double b, double Rcurv,
                double dx, double dy, double dt,
                dream_idx_type nt, double delay, double v, double cp,
                double *h,
                int err_level)
{
  dream_idx_type i, it;
  double t, xsmin, xsmax, ai, du, r;
  double phi, phi_min, phi_max;
  int err = NONE;

  if (b > 2*fabs(Rcurv)) {
    dream_err_msg("Error in dreamcylind: the y-size, b, must be less than the curvature diameter 2*Rcurv!\n");
  }

  b /= 2.0;
  double z_Rcurv = fabs(Rcurv) - sqrt(Rcurv*Rcurv - b*b);

  xsmin = -a/2.0;
  xsmax = a/2.0;

  phi = 2.0 * atan(b/(fabs(Rcurv)-z_Rcurv));
  phi_min = -phi/2;
  phi_max = phi/2;

  // dphi in y-dim [rad].
  double dphi = asin(dy/fabs(Rcurv));
  //dphi = dy/r;
  double ds = fabs(Rcurv) * dx * dphi;
  //ds = dx * dy; // Approx the same as  ds = Rcurv * dx * dphi.

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  phi = phi_min + dphi/2.0;
  double ys = fabs(Rcurv) * sin(phi);
  while (phi <= phi_max) {

    double xs = xsmin + dx/2.0;
    while (xs <= xsmax) {

      // Compute r and ds.
      if (Rcurv >= 0.0) {           // Focused
        cylind_f(xs, ys, Rcurv, z_Rcurv, xo, yo, zo, r, du);
      } else {                  // Defocused
        cylind_d(xs, ys, fabs(Rcurv), z_Rcurv, xo, yo, zo, r, du);
      }

      ai = v * ds * du/(2.0*M_PI * r);
      ai /= dt;
      // Convert to SI units.
      ai *= 1.0e3;
      // Propagation delay in micro seconds.
      t = r * 1.0e3/cp;
      it = (dream_idx_type) rint((t - delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {
        att.att(xc_vec, x_vec, r, it, h, ai);
      } else  {
        if  (it >= 0)
          err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
        else
          err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

        if ( (err_level == PARALLEL_STOP) || (err_level == STOP) ) {
          return err; // Bail out.
        }
      }

      xs += dx;
    }

    phi += dphi;
    ys = fabs(Rcurv) * sin(phi);
  }

  return err;
}


/***
 *
 * Focused
 *
 ***/

void cylind_f(double xs, double ys,
              double Rcurv, double z_Rcurv,
              double xo, double yo, double zo,
              double &r, double &du)
{
  du = 1.0;

  double zs = Rcurv - sqrt(Rcurv*Rcurv - ys*ys);

  double rx = xo - xs;
  double ry = yo - ys;
  double rz = zo - zs;
  r = sqrt(rx*rx + ry*ry + rz*rz);

  /* (**) FIXME: The code below protects about something : if we are inside the
     transducer? or perhaps if the angle is too large so we cannot reach
     the observation point? We disable these checks until we have figured out what
     they actually do!

  // Cos theta = -[rx ry rz]*[x y (z-r)]' / (R*r) (Scalar prd.). (focused)
  // Cos theta = [rx ry rz]*[x y (z+r)]' / (R*r) (Scalar prod.). (defocused)
  double cos_theta = -(rx*xs + ry*ys + rz*(zs-Rcurv)) / (Rcurv*r);

  if (zo <= z_Rcurv) {

    double dis1 = sqrt(xo*xo + yo*yo); // Horizontal distance from origo.
    double dis2 = sqrt(zo * (2.0*Rcurv - zo));

    if ( (cos_theta <  0.0) || (dis1 > dis2)) {
      du = (double) 0.0;
    }

  } else {

    if (cos_theta < (double) 0.0) {
      du = 0.0;
    }

  }
  */

  return;
}

/***
 *
 * Defocused
 *
 ***/


void cylind_d(double xs, double ys,
              double Rcurv, double z_Rcurv,
              double xo, double yo, double zo,
              double &r, double &du)
{
  du = 1.0;

  double zs = Rcurv - sqrt(Rcurv*Rcurv - ys*ys);

  double rx = xo - xs;
  double ry = yo - ys;
  double rz = zo + zs;          // zs is negative here.
  r = sqrt(rx*rx + ry*ry + rz*rz);

  /* FIXME see (**) above!
  // Cos theta = -[rx ry rz]*[x y (z-r)]' / (R*r) (Scalar prod.). (focused)
  // Cos theta = [rx ry rz]*[x y (z+r)]' / (R*r) (Scalar prod.). (defocused)
  double cos_theta = (rx*xs + ry*ys + rz*(z+Rcurv)) / (Rcurv*r);

  if (cos_theta < (double) 0.0) {
    du = 0.0;
  }
  */

  return;
}
