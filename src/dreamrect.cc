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

#include "dreamrect.h"
#include "attenuation.h"
#include "affinity.h"

#ifdef OCTAVE
#include <octave/oct.h>
#endif

// NB. We link this one in dream_arr_rect so we need unique names.
std::mutex err_mutex_rect;
std::atomic<bool> running_rect;
std::atomic<bool> verbose_err_msg_rect;

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
  SIRError err; // Output error.
} DATA_RECT;

/***
 *
 * Thread function.
 *
 ***/

void* Rect::smp_dream_rect(void *arg)
{
  DATA_RECT *D = (DATA_RECT *) arg;
  double xo, yo, zo;
  double *h = D->h;
  double a=D->a, b=D->b, dx=D->dx, dy=D->dy, dt=D->dt;
  dream_idx_type n, No=D->No, nt=D->nt;
  double *delay=D->delay, *Ro=D->Ro, v=D->v, cp=D->cp;
  Attenuation *att = D->att;
  dream_idx_type start=D->start, stop=D->stop;
  ErrorLevel err_level = D->err_level;
  SIRError err = SIRError::none;

  D->err = SIRError::none;  // Default to no output error.

  // Buffers for the FFTs in the Attenuation
  std::unique_ptr<FFTCVec> xc_vec;
  std::unique_ptr<FFTVec> x_vec;
  if (att) {
    xc_vec = std::make_unique<FFTCVec>(nt);
    x_vec = std::make_unique<FFTVec>(nt);
  }

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

#ifdef OCTAVE
    // Octave throws an exception when pressing CTRL-C
    // so catch it here and set running to false.
    try {
      OCTAVE_QUIT;
    }

    catch (octave::interrupt_exception &e) {
      running_rect = false;
    }

    catch (int &signum) {;}
#endif

    if (att == nullptr) {
      err = dreamrect_serial(xo, yo, zo,
                             a, b,
                             dx, dy, dt,
                             nt, dlay, v, cp,
                             &h[n*nt], err_level);
    } else {
      err = dreamrect_serial(*att, *xc_vec.get(),*x_vec.get(),
                             xo, yo, zo,
                             a, b,
                             dx, dy, dt,
                             nt, dlay, v, cp,
                             &h[n*nt], err_level);
    }

    if (err != SIRError::none) {
      D->err = err;
    }

    if (err == SIRError::out_of_bounds) {
      running_rect = false; // Tell all threads to exit.
    }

    if (!running_rect) {
      if (verbose_err_msg_rect) {
        std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!\n";
      }
      return(NULL);
    }

  }

  return(NULL);
}

/***
 *
 * dreamrect - Threaded rectangular transducer function.
 *
 ***/

SIRError Rect::dreamrect(double alpha,
                         double *Ro, dream_idx_type No,
                         double a, double b,
                         double dx, double dy, double dt, dream_idx_type nt,
                         DelayType delay_type, double *delay,
                         double v, double cp,
                         double *h, ErrorLevel err_level)
{
  SIRError err = SIRError::none;
  verbose_err_msg_rect = false;

  running_rect = true;

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11).
  dream_idx_type nthreads = std::thread::hardware_concurrency();

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
  DATA_RECT *D = (DATA_RECT*) malloc(nthreads*sizeof(DATA_RECT));

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
    D[thread_n].err = SIRError::none;

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
    for (dream_idx_type thread_n = 0; thread_n < nthreads; thread_n++) {
      threads[thread_n].join();

      // Check if one of the threads had an out-of-bounds event.
      if (D[thread_n].err != SIRError::none) {
        err = D[thread_n].err;
      }

    }
  }

  // Free memory.
  if (D) {
    free((void*) D);
  }

  return err;
}

SIRError dreamrect_serial(double xo, double yo, double zo,
                          double a, double b,
                          double dx, double dy, double dt,
                          dream_idx_type nt,
                          double delay,
                          double v, double cp,
                          double *h,
                          ErrorLevel err_level,
                          double weight)
{
  SIRError err = SIRError::none;

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
          err = dream_out_of_bounds_err("SIR out of bounds", it-nt+1, err_level);
        } else {
          err = dream_out_of_bounds_err("SIR out of bounds", it, err_level);
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

SIRError dreamrect_serial(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
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
  SIRError err = SIRError::none;

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

        if  (it >= 0) {
          err = dream_out_of_bounds_err("SIR out of bounds", it-nt+1, err_level);
        } else {
          err = dream_out_of_bounds_err("SIR out of bounds", it, err_level);
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
