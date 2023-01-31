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
#include <atomic>

#include "dreamcylind.h"
#include "attenuation.h"
#include "affinity.h"

std::mutex err_mutex_cylind;
std::atomic<bool> running_cylind;

void Cylind::abort(int signum)
{
  running_cylind = false;
}

bool Cylind::is_running()
{
  return running_cylind;
}

// Thread data.
typedef struct
{
  dream_idx_type start;
  dream_idx_type stop;
  double *Ro;
  dream_idx_type No;
  double a;
  double b;
  double Rcurv;
  double dx;
  double dy;
  double dt;
  dream_idx_type nt;
  DelayType delay_type;
  double *delay;
  double v;
  double cp;
  Attenuation *att;
  double *h;
  ErrorLevel err_level;
} DATA_CYLIND;

/***
 *
 * Thread function.
 *
 ***/

void* Cylind::smp_dream_cylind(void *arg)
{
  ErrorLevel tmp_err=ErrorLevel::none, err=ErrorLevel::none;
  DATA_CYLIND D = *(DATA_CYLIND *)arg;
  double xo, yo, zo;
  double *h = D.h;
  double a=D.a, b=D.b, Rcurv=D.Rcurv, dx=D.dx, dy=D.dy, dt=D.dt;
  dream_idx_type n, No=D.No, nt=D.nt;
  ErrorLevel tmp_lev=ErrorLevel::none, err_level=D.err_level;
  double *delay=D.delay, *Ro=D.Ro, v=D.v, cp=D.cp;
  Attenuation *att = D.att;
  dream_idx_type start=D.start, stop=D.stop;

  // Buffers for the FFTs in the Attenuation
  std::unique_ptr<FFTCVec> xc_vec;
  std::unique_ptr<FFTVec> x_vec;
  if (att) {
    xc_vec = std::make_unique<FFTCVec>(nt);
    x_vec = std::make_unique<FFTVec>(nt);
  }

  // Let the thread finish and then catch the error.
  if (err_level == ErrorLevel::stop)
    tmp_lev = ErrorLevel::parallel_stop;
  else
    tmp_lev = err_level;

  for (n=start; n<stop; n++) {
    xo = Ro[n];
    yo = Ro[n+1*No];
    zo = Ro[n+2*No];

    double dlay = 0.0;
    if (D.delay_type == DelayType::single) {
      dlay = delay[0];
    } else { // DelayType::multiple.
      dlay = delay[n];
    }

    if (att == nullptr) {
      err = dreamcylind_serial(xo, yo, zo,
                               a, b, Rcurv,
                               dx, dy, dt,
                               nt, dlay, v, cp,
                               &h[n*nt],tmp_lev);
    } else {
      err = dreamcylind_serial(*att, *xc_vec.get(),*x_vec.get(),
                               xo, yo, zo,
                               a, b, Rcurv,
                               dx, dy, dt,
                               nt, dlay, v, cp,
                               &h[n*nt],tmp_lev);
    }

    if (err != ErrorLevel::none || m_out_err ==  ErrorLevel::parallel_stop) {
        tmp_err = err;
        if (err == ErrorLevel::parallel_stop || m_out_err ==  ErrorLevel::parallel_stop)
          break; // Jump out when a ErrorLevel::stop error occurs.
    }

    if (!running_cylind) {
      std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!\n";
      return(NULL);
      }

  }

  if ((tmp_err != ErrorLevel::none) && (m_out_err == ErrorLevel::none)) {
    std::lock_guard<std::mutex> lk(err_mutex_cylind);

    m_out_err = tmp_err;
  }

  return(NULL);
}

/***
 *
 * dreamcylind : Threaded cylindric concave or convex transducer function
 *
 ***/

ErrorLevel Cylind::dreamcylind(double alpha,
                               double *Ro, dream_idx_type No,
                               double a, double b, double Rcurv,
                               double dx, double dy, double dt, dream_idx_type nt,
                               DelayType delay_type, double *delay,
                               double v, double cp,
                               double *h, ErrorLevel err_level)
{
  std::thread *threads;
  dream_idx_type thread_n, nthreads;
  dream_idx_type start, stop;
  DATA_CYLIND *D;

  running_cylind = true;

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11).
  nthreads = std::thread::hardware_concurrency();

  // Read DREAM_NUM_THREADS env var
  if(const char* env_p = std::getenv("DREAM_NUM_THREADS")) {
    dream_idx_type dream_threads = (dream_idx_type) std::stoul(env_p);
    if (dream_threads < nthreads) {
      nthreads = dream_threads;
    }
  }

  // nthreads can't be larger then the number of observation points.
  if (nthreads > No) {
    nthreads = No;
  }

  // Check if we have attenuation
  Attenuation att(nt, dt, alpha);
  Attenuation *att_ptr = nullptr;
  if (alpha > std::numeric_limits<double>::epsilon() ) {
    att_ptr = &att;
  }
  // Allocate local data.
  D = (DATA_CYLIND*) malloc(nthreads*sizeof(DATA_CYLIND));

  // Allocate mem for the threads.
  threads = new std::thread[nthreads]; // Init thread data.

  for (thread_n = 0; thread_n < nthreads; thread_n++) {

    start = thread_n * No/nthreads;
    stop =  (thread_n+1) * No/nthreads;

    // Init local data.
    D[thread_n].start = start; // Local start index;
    D[thread_n].stop = stop; // Local stop index;
    D[thread_n].No = No;
    D[thread_n].Ro = Ro;
    D[thread_n].a = a;
    D[thread_n].b = b;
    D[thread_n].Rcurv = Rcurv;
    D[thread_n].dx = dx;
    D[thread_n].dy = dy;
    D[thread_n].dt = dt;
    D[thread_n].nt = nt;
    D[thread_n].delay_type = delay_type;
    D[thread_n].delay = delay;
    D[thread_n].v = v;
    D[thread_n].cp = cp;
    D[thread_n].att = att_ptr;
    D[thread_n].h = h;
    D[thread_n].err_level = err_level;

    if (nthreads>1) {
      // Start the threads.
      threads[thread_n] = cylind_thread(&D[thread_n]);
      set_dream_thread_affinity(thread_n, nthreads, threads);
    } else {
      smp_dream_cylind(&D[0]);
    }
  }

  // Wait for all threads to finish.
  if (nthreads>1) {
    for (thread_n = 0; thread_n < nthreads; thread_n++) {
      threads[thread_n].join();
    }
  }

  // Free memory.
  if (D) {
    free((void*) D);
  }

  return m_out_err;
}

ErrorLevel Cylind::dreamcylind_serial(double xo, double yo, double zo,
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

ErrorLevel Cylind::dreamcylind_serial(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
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

/***
 *
 * Focused
 *
 ***/

double Cylind::cylind_f(double xs, double ys,
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

double Cylind::cylind_d(double xs, double ys,
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
