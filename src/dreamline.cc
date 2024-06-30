/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2012,2014,2021,2023,2024 Fredrik Lingvall
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

#include "dreamline.h"
#include "attenuation.h"
#include "affinity.h"

std::mutex err_mutex;
std::atomic<bool> running;
std::atomic<bool> verbose_err_messages;

void Line::abort(int signum)
{
  running = false;
}

bool Line::is_running()
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
  double a;
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

void* Line::smp_dream_line(void *arg)
{
  DATA *D = (DATA *) arg;
  double xo, yo, zo;
  double *h = D->h;
  double a=D->a, dx=D->dx, dy=D->dy, dt=D->dt;
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

    if (att == nullptr) {
      err = dreamline_serial(xo,yo,zo,a,dx,dy,dt,nt,dlay,v,cp, &h[n*nt], err_level);
    } else {
      err = dreamline_serial(*att, *xc_vec.get(),*x_vec.get(),
                             xo,yo,zo,a,dx,dy,dt,nt,dlay,v,cp, &h[n*nt], err_level);
    }

    if (err != SIRError::none) {
      D->err = err;
    }

    if (err == SIRError::out_of_bounds) {
      running = false; // Tell all threads to exit.
    }

    if (!running) {
      if (verbose_err_messages) {
        std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!\n";
      }
      return(NULL);
    }
  }

  return(NULL);
}

/***
 *
 *  dreamline - Threaded line (slit) trandsducer function
 *
 ***/

SIRError Line::dreamline(double alpha,
                         double *Ro, dream_idx_type No,
                         double a,
                         double dx, double dy, double dt, dream_idx_type nt,
                         DelayType delay_type, double *delay,
                         double v, double cp,
                         double *h, ErrorLevel err_level)
{
  SIRError err = SIRError::none;
  verbose_err_messages = false;

  running = true;

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
    D[thread_n].a = a;
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
    //D[thread_n].err = SIRError::none;

    if (nthreads>1) {
      // Start the threads.
      threads[thread_n] = line_thread(&D[thread_n]);
      set_dream_thread_affinity(thread_n, nthreads, threads);
    } else {
      smp_dream_line(&D[0]);
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

SIRError Line::dreamline_serial(double xo, double yo, double zo,
                                double a,
                                double dx, double dy, double dt, dream_idx_type nt,
                                double delay, double v,
                                double cp, double *h, ErrorLevel err_level)
{
  SIRError err = SIRError::none;

  dream_idx_type i, it;
  double t;
  double ai;
  double r;
  double xsmax = a/2.0;
  double xsmin = -a/2.0;

  // dy = width;
  double ds = dx * dy;

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  double xs = xsmin + dx / 2.0;

  while (xs <= xsmax) {

    r = std::sqrt((xo-xs)*(xo-xs) + yo*yo + zo*zo);

    ai = v * ds / (2*M_PI * r);
    ai /= dt;
    ai *= 1000.0;     // Convert to SI units.

    // Propagation delay in micro seconds.
    t = r * 1000.0/cp;
    it = (dream_idx_type) std::rint((t - delay)/dt);

    // Check if index is out of bounds.
    if ((it < nt) && (it >= 0)) {
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

  return err;
}

SIRError Line::dreamline_serial(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                                double xo, double yo, double zo,
                                double a,
                                double dx, double dy, double dt, dream_idx_type nt,
                                double delay, double v,
                                double cp, double *h, ErrorLevel err_level)
{
  SIRError err = SIRError::none;

  dream_idx_type i, it;
  double t;
  double ai;
  double r;
  double xsmax = a/2.0;
  double xsmin = -a/2.0;

  // dy = width;
  double ds = dx * dy;

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }

  double xs = xsmin + dx / 2.0;

  while (xs <= xsmax) {

    r = std::sqrt((xo-xs)*(xo-xs) + yo*yo + zo*zo);

    ai = v * ds / (2*M_PI * r);
    ai /= dt;
    ai *= 1000.0;     // Convert to SI units.

    // Propagation delay in micro seconds.
    t = r * 1000.0/cp;
    it = (dream_idx_type) std::rint((t - delay)/dt);

    // Check if index is out of bounds.
    if ((it < nt) && (it >= 0)) {

      att.att(xc_vec, x_vec, r, it, h, ai);

    } else  {

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

  return err;
}
