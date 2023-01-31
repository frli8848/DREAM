/***
 *
 * Copyright (C) 2002,2003,2006,2007,2008,2009,2014,2019,2021,2023 Fredrik Lingvall
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

#include "dreamrect.h"
#include "attenuation.h"
#include "affinity.h"

std::mutex err_mutex_rect;
std::atomic<bool> running_rect;

void Rect::abort(int signum)
{
  running_rect = false;
}

bool Rect::is_running()
{
  return running_rect;
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
} DATA_RECT;

/***
 *
 * Thread function.
 *
 ***/

void* Rect::smp_dream_rect(void *arg)
{
  ErrorLevel tmp_err=ErrorLevel::none, err=ErrorLevel::none;
  DATA_RECT D = *(DATA_RECT *)arg;
  double xo, yo, zo;
  double *h = D.h;
  double a=D.a, b=D.b, dx=D.dx, dy=D.dy, dt=D.dt;
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
  if (err_level == ErrorLevel::stop) {
    tmp_lev = ErrorLevel::parallel_stop;
  } else {
    tmp_lev = err_level;
  }

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
      err = dreamrect_serial(xo, yo, zo,
                             a, b,
                             dx, dy, dt,
                             nt, dlay, v, cp,
                             &h[n*nt], tmp_lev);
    } else {
      err = dreamrect_serial(*att, *xc_vec.get(),*x_vec.get(),
                             xo, yo, zo,
                             a, b,
                             dx, dy, dt,
                             nt, dlay, v, cp,
                             &h[n*nt], tmp_lev);
    }

    if (err != ErrorLevel::none || m_out_err ==  ErrorLevel::parallel_stop) {
      tmp_err = err;
      if (err == ErrorLevel::parallel_stop || m_out_err ==  ErrorLevel::parallel_stop) {
        break; // Jump out when a ErrorLevel::stop error occurs.
      }
    }

    if (!running_rect) {
      std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!\n";
      return(NULL);
    }
  }

  if ((tmp_err != ErrorLevel::none) && (m_out_err == ErrorLevel::none)) {
    std::lock_guard<std::mutex> lk(err_mutex_rect);

    m_out_err = tmp_err;
  }

  return(NULL);
}

/***
 *
 * dreamrect - Threaded rectangular transducer function.
 *
 ***/

ErrorLevel Rect::dreamrect(double alpha,
                           double *Ro, dream_idx_type No,
                           double a, double b,
                           double dx, double dy, double dt, dream_idx_type nt,
                           DelayType delay_type, double *delay,
                           double v, double cp,
                           double *h, ErrorLevel err_level)
{
  std::thread *threads;
  dream_idx_type thread_n, nthreads;
  dream_idx_type start, stop;
  DATA_RECT *D;

  running_rect = true;

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
  if (nthreads >  No) {
    nthreads = No;
  }

  // Check if we have attenuation
  Attenuation att(nt, dt, alpha);
  Attenuation *att_ptr = nullptr;
  if (alpha > std::numeric_limits<double>::epsilon() ) {
    att_ptr = &att;
  }

  // Allocate local data.
  D = (DATA_RECT*) malloc(nthreads*sizeof(DATA_RECT));

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
      threads[thread_n] = rect_thread(&D[thread_n]);
      set_dream_thread_affinity(thread_n, nthreads, threads);
    } else {
      smp_dream_rect(&D[0]);
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

ErrorLevel dreamrect_serial(double xo, double yo, double zo,
                            double a, double b,
                            double dx, double dy, double dt,
                            dream_idx_type nt,
                            double delay,
                            double v, double cp,
                            double *h,
                            ErrorLevel err_level,
                            double weight)
{
  ErrorLevel err = ErrorLevel::none;

  double xsmin = -a/2.0;
  double xsmax =  a/2.0;
  double ysmin = -b/2.0;
  double ysmax =  b/2.0;

  double ds = dx * dy;
  double rz = zo;
  double ys = ysmin + dy / 2.0;

  while (ys <= ysmax) {

    double ry = yo - ys;
    double xs = xsmin + dx / 2.0;

    while (xs <= xsmax) {

      double rx = xo - xs;
      double r = std::sqrt(rx*rx + ry*ry + rz*rz);

      // FIXME: avoid division here.
      double ai = weight * v * ds / (2*M_PI * r);
      ai /= dt;
      ai *= 1000.0;		// Convert to SI units.

      double t = r * 1000.0/cp;	// Propagation delay in micro seconds.
      dream_idx_type it = (dream_idx_type) std::rint((t - delay)/dt); // Sample index.

      // Check if index is out of bounds.
      if ( (it < nt) && (it >= 0) ) {
        h[it] += ai;
      } else {

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
    ys += dy;
  }

  return err;
}

ErrorLevel dreamrect_serial(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                            double xo, double yo, double zo,
                            double a, double b,
                            double dx, double dy, double dt,
                            dream_idx_type nt,
                            double delay,
                            double v, double cp,
                            double *h,
                            ErrorLevel err_level,
                            double weight)
{
  ErrorLevel err = ErrorLevel::none;

  double xsmin = -a/2.0;
  double xsmax =  a/2.0;
  double ysmin = -b/2.0;
  double ysmax =  b/2.0;

  double ds = dx * dy;

  double rz = zo;
  double ys = ysmin + dy / 2.0;

  while (ys <= ysmax) {

    double ry = yo - ys;

    double xs = xsmin + dx / 2.0;

    while (xs <= xsmax) {

      double rx = xo - xs;
      double r = std::sqrt(rx*rx + ry*ry + rz*rz);

      double ai = weight* v * ds / (2*M_PI * r);
      ai /= dt;
      ai *= 1000.0;		// Convert to SI units.

      double t = r * 1000.0/cp;	// Propagation delay in micro seconds.
      dream_idx_type it = (dream_idx_type) std::rint((t - delay)/dt); // Sample index.

      // Check if index is out of bounds.
      if ( (it < nt) && (it >= 0) ) {
        att.att(xc_vec, x_vec, r, it, h, ai);
      }  else {
        if  (it >= 0)
          err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
        else
          err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

        if ( (err_level == ErrorLevel::parallel_stop) || (err_level == ErrorLevel::stop) )
          return err; // Bail out.
      }
      xs += dx;
    }
    ys += dy;
  }

  return err;
}
