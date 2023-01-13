/***
 *
 * Copyright (C) 2002,2003,2005,2006,2007,2008,2009,2014,2019,2021,2023 Fredrik Lingvall
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

#include "dreamcirc.h"
#include "attenuation.h"
#include "affinity.h"

std::mutex err_mutex_circ;
std::atomic<bool> running_circ;

void Circ::abort(int signum)
{
  running_circ = false;
}

bool Circ::is_running()
{
  return running_circ;
}

// Thread data.
typedef struct
{
  dream_idx_type no;
  dream_idx_type start;
  dream_idx_type stop;
  double *ro;
  double R;
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
} DATA_CIRC;

/***
 *
 * Thread function.
 *
 ***/

void* Circ::smp_dream_circ(void *arg)
{
  ErrorLevel tmp_err=ErrorLevel::none, err=ErrorLevel::none;
  DATA_CIRC D = *(DATA_CIRC *)arg;
  double xo, yo, zo;
  double *h = D.h;
  double R=D.R, dx=D.dx, dy=D.dy, dt=D.dt;
  dream_idx_type n, no=D.no, nt=D.nt;
  ErrorLevel tmp_lev=ErrorLevel::none, err_level=D.err_level;
  double *delay=D.delay, *ro=D.ro, v=D.v, cp=D.cp;
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
    xo = ro[n];
    yo = ro[n+1*no];
    zo = ro[n+2*no];

    double dlay = 0.0;
    if (D.delay_type == DelayType::single) {
      dlay = delay[0];
    } else { // DelayType::multiple.
      dlay = delay[n];
    }

    if (att == nullptr) {
      err = dreamcirc_serial(xo, yo, zo,
                             R,
                             dx, dy, dt,
                             nt, dlay, v, cp,
                             &h[n*nt], tmp_lev);
    } else {
      err = dreamcirc_serial(*att, *xc_vec.get(),*x_vec.get(),
                             xo, yo, zo,
                             R,
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

    if (!running_circ) {
      std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!\n";
      return(NULL);
    }
  }

  if ((tmp_err != ErrorLevel::none) && (m_out_err == ErrorLevel::none)) {
    std::lock_guard<std::mutex> lk(err_mutex_circ);

    m_out_err = tmp_err;
  }

  return(NULL);
}

/***
 *
 * dreamcirc - Threaded circular transducer function.
 *
 ***/

ErrorLevel Circ::dreamcirc(double alpha,
                           double *ro, dream_idx_type no,
                           double R,
                           double dx, double dy, double dt, dream_idx_type nt,
                           DelayType delay_type, double *delay,
                           double v, double cp,
                           double *h, ErrorLevel err_level)
{
  std::thread *threads;
  unsigned int thread_n, nthreads;
  dream_idx_type start, stop;
  DATA_CIRC *D;

  running_circ = true;

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11).
  nthreads = std::thread::hardware_concurrency();

  // Read DREAM_NUM_THREADS env var
  if(const char* env_p = std::getenv("DREAM_NUM_THREADS")) {
    unsigned int dream_threads = std::stoul(env_p);
    if (dream_threads < nthreads) {
      nthreads = dream_threads;
    }
  }

  // nthreads can't be larger then the number of observation points.
  if (nthreads > (unsigned int) no) {
    nthreads = no;
  }

  // Check if we have attenuation
  Attenuation att(nt, dt, alpha);
  Attenuation *att_ptr = nullptr;
  if (alpha > std::numeric_limits<double>::epsilon() ) {
    att_ptr = &att;
  }

  // Allocate local data.
  D = (DATA_CIRC*) malloc(nthreads*sizeof(DATA_CIRC));

  // Allocate mem for the threads.
  threads = new std::thread[nthreads]; // Init thread data.

  for (thread_n = 0; thread_n < nthreads; thread_n++) {

    start = thread_n * no/nthreads;
    stop =  (thread_n+1) * no/nthreads;

    // Init local data.
    D[thread_n].start = start; // Local start index;
    D[thread_n].stop = stop; // Local stop index;
    D[thread_n].no = no;
    D[thread_n].ro = ro;
    D[thread_n].R = R;
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

    // Start the threads.
    threads[thread_n] = circ_thread(&D[thread_n]);
    set_dream_thread_affinity(thread_n, nthreads, threads);
  }

  // Wait for all threads to finish.
  for (thread_n = 0; thread_n < nthreads; thread_n++) {
    threads[thread_n].join();
  }

  // Free memory.
  if (D) {
    free((void*) D);
  }

  return m_out_err;
}


/***
 *
 *  dreamcirc
 *
 * Computes the spatial impulse response of a circular aperture.
 *
 ***/

ErrorLevel dreamcirc_serial(double xo, double yo, double zo,
                            double R,
                            double dx, double dy, double dt,
                            dream_idx_type nt,
                            double delay,
                            double v, double cp,
                            double *h,
                            ErrorLevel err_level,
                            double weight)
{
  ErrorLevel err = ErrorLevel::none;

  double ds = dx * dy;

  // y-dim integration bounds
  double ysmin = -R;
  double ysmax =  R;

  double ys = ysmin + dy/2.0;
  while (ys <= ysmax) {

    // x-dim integration bounds
    double rxs = std::sqrt(R*R - ys*ys);
    double xsmin = -rxs;
    double xsmax = rxs;

    double ry = yo - ys;

    double xs = xsmin + dx/2.0;
    while (xs <= xsmax) {

      // Compute the distance (length) from an observation point (xo,yo,zo)
      // to a point (xs,ys) on the transducer surface.
      double rx = xo - xs;
      double r = std::sqrt(rx*rx + ry*ry + zo*zo);

      double ai = weight * v * ds / (2*M_PI * r);
      ai /= dt;
      ai *= 1000.0; // Convert to SI units.

      // Propagation delay in micro seconds.
      double t = r * 1000.0/cp;
      dream_idx_type it = (dream_idx_type) std::rint((t - delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {
        h[it] += ai;
      } else {

        if (it >= 0) {
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

ErrorLevel dreamcirc_serial(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                            double xo, double yo, double zo,
                            double R,
                            double dx, double dy, double dt,
                            dream_idx_type nt,
                            double delay,
                            double v, double cp,
                            double *h,
                            ErrorLevel err_level,
                            double weight)
{
  ErrorLevel err = ErrorLevel::none;

  double ds = dx * dy;

  // y-dim integration bounds
  double ysmin = -R;
  double ysmax =  R;

  double ys = ysmin + dy/2.0;
  while (ys <= ysmax) {

    // x-dim integration bounds
    double rxs = std::sqrt(R*R - ys*ys);
    double xsmin = -rxs;
    double xsmax = rxs;

    double ry = yo - ys;

    double xs = xsmin + dx/2.0;
    while (xs <= xsmax) {

      // Compute the distance (length) from an observation point (xo,yo,zo)
      // to a point (xs,ys) on the transducer surface.
      double rx = xo - xs;
      double r = std::sqrt(rx*rx + ry*ry + zo*zo);

      double ai = weight * v * ds / (2*M_PI * r);
      ai /= dt;
      ai *= 1000.0; // Convert to SI units.

      // Propagation delay in micro seconds.
      double t = r * 1000.0/cp;
      dream_idx_type it = (dream_idx_type) std::rint((t - delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {
        att.att(xc_vec, x_vec, r, it, h, ai);
      } else {

        if (it >= 0) {
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
