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

#include <cmath>

#include "dreamcylind.h"

//
//  Function prototypes.
//

double cylind_f(double xs, double ys,
                double Rcurv,
                double z_Rcurv,
                double xo, double yo, double zo,
                double &du);

double cylind_d(double xs, double ys,
                double Rcurv, double z_Rcurv,
                double xo, double yo, double zo,
                double &du);

/***
 *
 * dreamcylind : cylindric concave or convex transducer
 *
 ***/

ErrorLevel dreamcylind(double xo, double yo, double zo,
                       double a, double b, double Rcurv,
                       double dx, double dy, double dt,
                       dream_idx_type nt, double delay, double v, double cp,
                       double *h,
                       ErrorLevel err_level,
                       double weight)
{
  double t, xsmin, xsmax, ai, du, r;
  double phi, phi_min, phi_max;
  ErrorLevel err = ErrorLevel::none;

  if (b > 2.0*std::fabs(Rcurv)) {
    dream_err_msg("Error in dreamcylind: the y-size, b, must be less than the curvature diameter 2*Rcurv!\n");
  }

  b /= 2.0;
  double z_Rcurv = std::fabs(Rcurv) - std::sqrt(Rcurv*Rcurv - b*b);

  xsmin = -a/2.0;
  xsmax = a/2.0;

  phi = 2.0 * std::atan(b/(std::fabs(Rcurv)-z_Rcurv));
  phi_min = -phi/2.0;
  phi_max = phi/2.0;

  // dphi in y-dim [rad].
  double dphi = std::asin(dy/std::fabs(Rcurv));
  //dphi = dy/r;
  double ds = std::fabs(Rcurv) * dx * dphi;
  //ds = dx * dy; // Approx the same as  ds = Rcurv * dx * dphi.

  phi = phi_min + dphi/2.0;
  double ys = std::fabs(Rcurv) * std::sin(phi);
  while (phi <= phi_max) {

    double xs = xsmin + dx/2.0;
    while (xs <= xsmax) {

      // Compute r and ds.
      if (Rcurv >= 0.0) {       // Focused
        r = cylind_f(xs, ys, Rcurv, z_Rcurv, xo, yo, zo, du);
      } else {                  // Defocused
        r = cylind_d(xs, ys, std::fabs(Rcurv), z_Rcurv, xo, yo, zo, du);
      }

      ai = weight * v * ds * du/(2.0*M_PI * r);
      ai /= dt;
      // Convert to SI units.
      ai *= 1.0e3;
      // Propagation delay in micro seconds.
      t = r * 1.0e3/cp;
      dream_idx_type it = (dream_idx_type) std::rint((t - delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {
        h[it] += ai;
      } else  {
        if  (it >= 0) {
          err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
        } else {
          err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);
        }

        if ( (err_level == ErrorLevel::parallel_stop) || (err_level == ErrorLevel::stop) ) {
          return err; // Bail out.
        }
      }

      xs += dx;
    }

    phi += dphi;
    ys = std::fabs(Rcurv) * std::sin(phi);
  }

  return err;
}

ErrorLevel dreamcylind(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                       double xo, double yo, double zo,
                       double a, double b, double Rcurv,
                       double dx, double dy, double dt,
                       dream_idx_type nt, double delay, double v, double cp,
                       double *h,
                       ErrorLevel err_level,
                       double weight)
{
  dream_idx_type it;
  double t, xsmin, xsmax, ai, du, r;
  double phi, phi_min, phi_max;
  ErrorLevel err = ErrorLevel::none;

  if (b > 2.0*std::fabs(Rcurv)) {
    dream_err_msg("Error in dreamcylind: the y-size, b, must be less than the curvature diameter 2*Rcurv!\n");
  }

  b /= 2.0;
  double z_Rcurv = std::fabs(Rcurv) - std::sqrt(Rcurv*Rcurv - b*b);

  xsmin = -a/2.0;
  xsmax = a/2.0;

  phi = 2.0 * std::atan(b/(std::fabs(Rcurv)-z_Rcurv));
  phi_min = -phi/2.0;
  phi_max = phi/2.0;

  // dphi in y-dim [rad].
  double dphi = std::asin(dy/std::fabs(Rcurv));
  //dphi = dy/r;
  double ds = std::fabs(Rcurv) * dx * dphi;
  //ds = dx * dy; // Approx the same as  ds = Rcurv * dx * dphi.

  phi = phi_min + dphi/2.0;
  double ys = std::fabs(Rcurv) * std::sin(phi);
  while (phi <= phi_max) {

    double xs = xsmin + dx/2.0;
    while (xs <= xsmax) {

      // Compute r and ds.
      if (Rcurv >= 0.0) {           // Focused
        r = cylind_f(xs, ys, Rcurv, z_Rcurv, xo, yo, zo, du);
      } else {                  // Defocused
        r = cylind_d(xs, ys, std::fabs(Rcurv), z_Rcurv, xo, yo, zo, du);
      }

      ai = weight * v * ds * du/(2.0*M_PI * r);
      ai /= dt;
      // Convert to SI units.
      ai *= 1.0e3;
      // Propagation delay in micro seconds.
      t = r * 1.0e3/cp;
      it = (dream_idx_type) std::rint((t - delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {
        att.att(xc_vec, x_vec, r, it, h, ai);
      } else  {
        if  (it >= 0)
          err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
        else
          err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

        if ( (err_level == ErrorLevel::parallel_stop) || (err_level == ErrorLevel::stop) ) {
          return err; // Bail out.
        }
      }

      xs += dx;
    }

    phi += dphi;
    ys = std::fabs(Rcurv) * std::sin(phi);
  }

  return err;
}

/***
 *
 * Focused
 *
 ***/

double cylind_f(double xs, double ys,
                double Rcurv, double z_Rcurv,
                double xo, double yo, double zo,
                double &du)
{
  du = 1.0;

  double zs = Rcurv - std::sqrt(Rcurv*Rcurv - ys*ys);

  double rx = xo - xs;
  double ry = yo - ys;
  double rz = zo - zs;
  double r = std::sqrt(rx*rx + ry*ry + rz*rz);

  /* (**) FIXME: The code below protects about something : if we are inside the
     transducer? or perhaps if the angle is too large so we cannot reach
     the observation point? We disable these checks until we have figured out what
     they actually do!

  // Cos theta = -[rx ry rz]*[x y (z-r)]' / (R*r) (Scalar prd.). (focused)
  // Cos theta = [rx ry rz]*[x y (z+r)]' / (R*r) (Scalar prod.). (defocused)
  double cos_theta = -(rx*xs + ry*ys + rz*(zs-Rcurv)) / (Rcurv*r);

  if (zo <= z_Rcurv) {

    double d1 = std::sqrt(xo*xo + yo*yo); // Horizontal distance from origo.
    double d2 = std::sqrt(zo * (2.0*Rcurv - zo));

    if ( (cos_theta <  0.0) || (d1 > d2)) {
      du = (double) 0.0;
    }

  } else {

    if (cos_theta < 0.0) {
      du = 0.0;
    }

  }
  */

  return r;
}

/***
 *
 * Defocused
 *
 ***/

double cylind_d(double xs, double ys,
                double Rcurv, double z_Rcurv,
                double xo, double yo, double zo,
                double &du)
{
  du = 1.0;

  double zs = Rcurv - std::sqrt(Rcurv*Rcurv - ys*ys);

  double rx = xo - xs;
  double ry = yo - ys;
  double rz = zo + zs;          // zs is negative here.
  double r = std::sqrt(rx*rx + ry*ry + rz*rz);

  /* FIXME see (**) above!
  // Cos theta = -[rx ry rz]*[x y (z-r)]' / (R*r) (Scalar prod.). (focused)
  // Cos theta = [rx ry rz]*[x y (z+r)]' / (R*r) (Scalar prod.). (defocused)
  double cos_theta = (rx*xs + ry*ys + rz*(z+Rcurv)) / (Rcurv*r);

  if (cos_theta < (double) 0.0) {
    du = 0.0;
  }
  */

  return r;
}
