/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2014,2019,2021,2023,2024 Fredrik Lingvall
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

#include "dreamsphere.h"
#include "attenuation.h"
#include "affinity.h"

std::mutex err_mutex;
std::atomic<bool> running;

void Sphere::abort(int signum)
{
  running = false;
}

bool Sphere::is_running()
{
  return running;
}

// Thread data.
typedef struct
{
  dream_idx_type start;
  dream_idx_type stop;
  double *Ro;
  dream_idx_type No;
  double R;
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
  SIRError err; // Output error.
} DATA;

/***
 *
 * Thread function.
 *
 ***/

void* Sphere::smp_dream_sphere(void *arg)
{
  DATA *D = (DATA *)arg;
  double xo, yo, zo;
  double *h = D->h;
  double R=D->R, Rcurv=D->Rcurv, dx=D->dx, dy=D->dy, dt=D->dt;
  dream_idx_type n, No=D->No, nt=D->nt;
  double *delay=D->delay, *Ro=D->Ro, v=D->v, cp=D->cp;
  Attenuation *att = D->att;
  dream_idx_type start=D->start, stop=D->stop;
  ErrorLevel err_level = D->err_level;
  SIRError err = SIRError::none;

  // Buffers for the FFTs in the Attenuation
  std::unique_ptr<FFTCVec> xc_vec;
  std::unique_ptr<FFTVec> x_vec;
  if (att) {
    xc_vec = std::make_unique<FFTCVec>(nt);
    x_vec = std::make_unique<FFTVec>(nt);
  }

  D->err = SIRError::none;  // Default to no output error.

  for (n=start; n<stop; n++) {
    xo = Ro[n];
    yo = Ro[n+1*No];
    zo = Ro[n+2*No];

    double dlay = 0.0;
    if (D->delay_type == DelayType::single) {
      dlay = delay[0];
    } else { // DelayType::multiple.
      dlay = delay[n];
    }

    if (att == nullptr) {
      err = dreamsphere_serial(xo, yo, zo,
                               R, Rcurv,
                               dx, dy, dt,
                               nt, dlay, v, cp,
                               &h[n*nt], err_level);
    } else {
      err = dreamsphere_serial(*att, *xc_vec.get(),*x_vec.get(),
                               xo, yo, zo,
                               R, Rcurv,
                               dx, dy, dt,
                               nt, dlay, v, cp,
                               &h[n*nt], err_level);
    }

    if (err == SIRError::out_of_bounds) {
      D->err = err; // Return the out-of-bounds error for this thread.
      running = false;          // Tell all threads to exit.
    }

    if (!running) {
      std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!\n";
      return(NULL);
    }

  }

  return(NULL);
}

/***
 *
 *  dreamsphere - Threaded  focused (concave) / defucused (convex) spherical transducer function
 *
 ***/

SIRError Sphere::dreamsphere(double alpha,
                             double *Ro, dream_idx_type No,
                             double R, double Rcurv,
                             double dx, double dy, double dt, dream_idx_type nt,
                             DelayType delay_type, double *delay,
                             double v, double cp,
                             double *h, ErrorLevel err_level)
{
  dream_idx_type thread_n;
  dream_idx_type start, stop;

  SIRError err = SIRError::none;

  running = true;

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11)
  dream_idx_type nthreads = std::thread::hardware_concurrency();

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
  DATA *D = (DATA*) malloc(nthreads*sizeof(DATA));

  // Allocate mem for the threads.
  std::thread *threads = new std::thread[nthreads]; // Init thread data.

  for (dream_idx_type thread_n = 0; thread_n < nthreads; thread_n++) {

    dream_idx_type start = thread_n * No/nthreads;
    dream_idx_type stop =  (thread_n+1) * No/nthreads;

    // Init local data.
    D[thread_n].start = start; // Local start index;
    D[thread_n].stop = stop; // Local stop index;
    D[thread_n].No = No;
    D[thread_n].Ro = Ro;
    D[thread_n].R = R;
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
    D[thread_n].err = SIRError::none;

    if (nthreads>1) {
      // Start the threads.
      threads[thread_n] = sphere_thread(&D[thread_n]);
      set_dream_thread_affinity(thread_n, nthreads, threads);
    } else {
      smp_dream_sphere(&D[0]);
    }
  }

  // Wait for all threads to finish.
  if (nthreads>1) {
    for (dream_idx_type thread_n = 0; thread_n < nthreads; thread_n++) {
      threads[thread_n].join();

      // Check if the current thread or a previous had an out-of-bounds error.
      if ( (err == SIRError::out_of_bounds) || (D[thread_n].err == SIRError::out_of_bounds) ) {
        err = SIRError::out_of_bounds;
      }

    }
  }

  // Free memory.
  if (D) {
    free((void*) D);
  }

  return err;
}

SIRError Sphere::dreamsphere_serial(double xo, double yo, double zo,
                                    double R, double Rcurv,
                                    double dx, double dy, double dt,
                                    dream_idx_type nt, double delay, double v, double cp,
                                    double *h, ErrorLevel err_level)
{
  SIRError err = SIRError::none;

  dream_idx_type i, it;
  double t, ai, r;
  double xsmin, ysmin, xsmax, ysmax;

  double ds = dx * dy;

  ysmin = -R;
  ysmax =  R;

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  double ys = ysmin + dy / 2.0;
  while (ys <= ysmax) {

    // Compute the x-axis integration limits.
    double rs = std::sqrt(R*R - ys*ys);
    xsmin = -rs;
    xsmax =  rs;

    double xs = xsmin + dx / 2.0;
    while (xs <= xsmax) {

      if (Rcurv >= 0.0) {        // Concave/focused
        r = sphere_f(xs, ys,
                     Rcurv,
                     xo, yo, zo);
      } else {                  // Convex/de-focused
        r = sphere_d(xs, ys,
                     -Rcurv,
                     xo, yo, zo);
      }

      ai = v * ds / (2*M_PI * r);
      ai /= dt;
      ai *= 1000.0;              // Convert to SI units.

      // Propagation delay in micro seconds.
      t = r * 1000.0/cp;
      it = (dream_idx_type) std::rint((t - delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {
        h[it] += ai;
      } else {

        if  (it >= 0) {
          err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
        } else {
          err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);
        }

        if (err == SIRError::out_of_bounds) {
          return err; // Bail out.
        }

      }
      xs += dx;
    }
    ys += dy;
  }

  return err;
}

SIRError Sphere::dreamsphere_serial(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                                    double xo, double yo, double zo,
                                    double R, double Rcurv,
                                    double dx, double dy, double dt,
                                    dream_idx_type nt, double delay, double v, double cp,
                                    double *h, ErrorLevel err_level)
{
  SIRError err = SIRError::none;

  dream_idx_type i, it;
  double t, ai, r;
  double xsmin, ysmin, xsmax, ysmax;

  double ds = dx * dy;

  ysmin = -R;
  ysmax =  R;

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  double ys = ysmin + dy / 2.0;
  while (ys <= ysmax) {

    double rs = std::sqrt(R*R - ys*ys);
    xsmin = -rs;
    xsmax =  rs;

    double xs = xsmin + dx / 2.0;
    while (xs <= xsmax) {

      if (Rcurv >= 0.0) {        // Concave/focused
        r = sphere_f(xs, ys,
                     Rcurv,
                     xo, yo, zo);
      } else {                  // Convex/de-focused
        r = sphere_d(xs, ys,
                     -Rcurv,
                     xo, yo, zo);
      }

      ai = v * ds / (2*M_PI * r);
      ai /= dt;
      ai *= 1000.0;              // Convert to SI units.

      // Propagation delay in micro seconds.
      t = r * 1000.0/cp;
      it = (dream_idx_type) std::rint((t - delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {
        att.att(xc_vec, x_vec, r, it, h, ai);
      } else {

        if (it >= 0) {
          err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
        } else {
          err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);
        }

        if (err == SIRError::out_of_bounds) {
          return err; // Bail out.
        }

      }
      xs += dx;
    }
    ys += dy;
  }

  return err;
}

/***
 *
 * sphere_f - focused (concave) sphere
 *
 ***/

double Sphere::sphere_f(double xs, double ys,
                        double Rcurv,
                        double xo, double yo, double zo)
{
  double rx, ry, rz;

  double zs = Rcurv - std::sqrt(Rcurv*Rcurv - xs*xs - ys*ys);

  rx = xo - xs;
  ry = yo - ys;
  rz = zo - zs;
  double r = std::sqrt(rx*rx + ry*ry + rz*rz);

  return r;
}

/***
 *
 * sphere_d - de-focused (convex) sphere
 *
 ***/

double Sphere::sphere_d(double xs, double ys,
                        double Rcurv,
                        double xo, double yo, double zo)
{
  double rx, ry, rz;

  double zs = Rcurv - std::sqrt(Rcurv*Rcurv - xs*xs - ys*ys);

  rx = xo - xs;
  ry = yo - ys;
  rz = zo + zs;                // Change the sign of zs for defocused.
  double r = std::sqrt(rx*rx + ry*ry + rz*rz);

  return r;
}
