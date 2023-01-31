/***
*
* Copyright (C) 2004,2006,2007,2008,2009,2014,2021,2023 Fredrik Lingvall
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
#include <iostream>
#include <string>

#include "rect_sir.h"
#include "affinity.h"

//std::mutex err_mutex;
std::atomic<bool> running;

void RectSir::abort(int signum)
{
  running = false;
}

bool RectSir::is_running()
{
  return running;
}

// Thread data.
typedef struct
{
  dream_idx_type no;
  dream_idx_type start;
  dream_idx_type stop;
  double *ro;
  double a;
  double b;
  double dt;
  dream_idx_type nt;
  DelayType delay_type;
  double *delay;
  double v;
  double cp;
  //double alpha;
  double *h;
  ErrorLevel err_level;
} DATA;

/***
 *
 * Thread function.
 *
 ***/

void* RectSir::smp_rect_sir(void *arg)
{
  //ErrorLevel tmp_err=ErrorLevel::none, err=ErrorLevel::none;
  DATA D = *(DATA *)arg;
  double xo, yo, zo;
  double *h = D.h;
  double a=D.a, b=D.b, dt=D.dt;
  dream_idx_type n, no=D.no, nt=D.nt;
  //int    tmp_lev, err_level=D.err_level;
  double *delay=D.delay, *ro=D.ro, v=D.v, cp=D.cp; // alfa=D.alfa;
  dream_idx_type start=D.start, stop=D.stop;

  // Let the thread finish and then catch the error.

  /*
  if (err_level == ErrorLevel::stop) {
    tmp_lev = ErrorLevel::parallel_stop;
  } else {
    tmp_lev = err_level;
  }
  */

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

    // NB. rect_sir_serial do not return any error codes!
    rect_sir_serial(xo,yo,zo,a,b,dt,nt,dlay,v,cp,&h[n*nt]); // TODO: Add attenuation.

    /*
    if (err != ErrorLevel::none || m_out_err == ErrorLevel::parallel_stop) {
      tmp_err = err;
      if (err == ErrorLevel::parallel_stop || m_out_err ==  ErrorLevel::parallel_stop) {
        break; // Jump out when a ErrorLevel::stop error occurs.
      }
    }
    */
    if (!running) {
      std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!\n";
      return(NULL);
    }
  }

  /*
  if ((tmp_err != ErrorLevel::none) && (m_out_err == ErrorLevel::none)) {
    std::lock_guard<std::mutex> lk(err_mutex);

    m_out_err = tmp_err;
  }
  */

  return(NULL);
}

/***
 *
 * rect_sir - Threaded analytic rectangular transducer function.
 *
 ***/

ErrorLevel RectSir::rect_sir(double *ro, dream_idx_type no,
                             double a, double b,
                             double dt, dream_idx_type nt,
                             DelayType delay_type, double *delay,
                             double v, double cp,
                             double *h)
{
  std::thread *threads;
  unsigned int thread_n, nthreads;
  dream_idx_type start, stop;
  DATA *D;

  running = true;

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

  // Allocate local data.
  D = (DATA*) malloc(nthreads*sizeof(DATA));

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
    D[thread_n].a = a;
    D[thread_n].b = b;
    D[thread_n].dt = dt;
    D[thread_n].nt = nt;
    D[thread_n].delay_type = delay_type;
    D[thread_n].delay = delay;
    D[thread_n].v = v;
    D[thread_n].cp = cp;
    //D[thread_n].alpha = alpha;
    D[thread_n].h = h;
    //D[thread_n].err_level = err_level;

    if (nthreads>1) {
      // Start the threads.
      threads[thread_n] = rect_sir_thread(&D[thread_n]);
      set_dream_thread_affinity(thread_n, nthreads, threads);
    } else {
      smp_rect_sir(&D[0]);
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

/***
 *
 *  Subroutine rect_sir - Continous spatial impulse response for a rectangular aperture.
 *
 *  See "High-speed method for computing the exact solution for the pressure variations
 *  in the nearfield of a baffled piston", J.C. Lockwood and J.G. Willette, J. Acoust. Soc. Am.,
 *  year = {1973},  volume =  {53},  number = {3},  pages =  {735--741}, month =  {March}.
 *
 ***/

void RectSir::rect_sir_serial(double xo_i,
                              double yo_i,
                              double zo_i,
                              double a_i,
                              double b_i,
                              double dt,
                              dream_idx_type nt,
                              double delay,
                              double v,
                              double cp,
                              double *h)
{
  dream_idx_type it, k;
  double t, t_z, xo, yo, zo;
  double tau_1, tau_2, tau_3, tau_4, a, b;
  double a_k=0, g_k=0, s_k=0, l_k=0;

  // Convert to [m].
  xo = fabs(xo_i) / 1000.0;	// Can take abs due to symmetry.
  yo = fabs(yo_i) / 1000.0;	// Can take abs due to symmetry.
  zo = zo_i / 1000.0;
  a = a_i / 1000.0;
  b = b_i / 1000.0;

  for (it = 0; it < nt; it++) {
    h[it] = (double) 0.0;
  }

  for  (k=1; k<=4; k++) { // loop over all 4 sub-rectangles.

    if ( (xo <= a/2.0) && (yo <= b/2.0) ) {// Inside the aperture.

      g_k = 1; // Add all.
      switch (k) {

      case 1:
        s_k = fabs(a/2.0 - xo);
        l_k = fabs(b/2.0 - yo);
        break;

      case 2:
        s_k = fabs(xo + a/2.0);
        l_k = fabs(b/2.0 - yo);
        break;

      case 3:
        s_k = fabs(xo + a/2.0);
        l_k = fabs(yo + b/2.0);
        break;

      case 4:
        s_k = fabs(a/2.0 - xo);
        l_k = fabs(b/2.0 - yo);
        break;

      default:
        break;
      }
    }

    if ( (xo <= a/2.0) && (yo > b/2.0) ) {// Inside a/2 but outside b/2.

      switch (k) {

      case 1:
        s_k = fabs(a/2.0 - xo);
        l_k = fabs(yo - b/2.0);
        g_k = -1; // Outside => subtract.
        break;

      case 2:
        s_k = fabs(a/2.0 + xo);
        l_k = fabs(yo - b/2.0);
        g_k = -1; // Outside => subtract.
        break;

      case 3:
        s_k = fabs(a/2.0 + xo);
        l_k = fabs(b/2.0 + yo);
        g_k = 1; // Inside => add.
        break;

      case 4:
        s_k = fabs(a/2.0 - xo);
        l_k = fabs(b/2.0 + yo);
        g_k = 1; // Inside => add.
        break;

      default:
        break;
      }
    } // if inside a/2 but outside b/2.

    if ( (xo > a/2.0) && (yo <= b/2.0) ) { // Inside b/2 but outside a/2.

      switch (k) {

      case 1:
        s_k = fabs(xo - a/2.0);
        l_k = fabs(b/2.0 - yo);
        g_k = -1; // Outside => subtract.
        break;

      case 2:
        s_k = fabs(a/2.0 + xo );
        l_k = fabs(b/2.0 - yo);
        g_k = 1; // Inside => add.
        break;

      case 3:
        s_k = fabs(a/2.0 + xo);
        l_k = fabs(b/2.0 + yo);
        g_k = 1; // Inside => add.
        break;

      case 4:
        s_k = fabs(xo - a/2.0);
        l_k = fabs(b/2.0 + yo);
        g_k = -1; // Outside => subtract.
        break;

      default:
        break;
      }
    } // if inside b/2 but outside a/2.

    if ( (xo > a/2.0) && (yo > b/2.0) ) {// Outside both  a/2 and b/2.
      switch (k) {

      case 1:
        s_k = fabs(xo - a/2.0);
        l_k = fabs(yo - b/2.0);
        g_k = 1; // Really outside but need to add since 1 is included both in 2 and 4.
        break;

      case 2:
        s_k = fabs(a/2.0 + xo);
        l_k = fabs(yo - b/2.0);
        g_k = -1; // Outside => subtract.
        break;

      case 3:
        s_k = fabs(a/2.0 + xo);
        l_k = fabs(b/2.0 + yo);
        g_k = 1; // Inside => add.
        break;

      case 4:
        s_k = fabs(xo - a/2.0);
        l_k = fabs(yo + b/2.0);
        g_k = -1; // Outside => subtract.
        break;

      default:
        break;
      }
    } // if outside both  a/2 and b/2.

    t_z = zo / cp;

    tau_1 = t_z;
    tau_2 = std::sqrt(zo*zo + s_k*s_k) / cp;
    tau_3 = std::sqrt(zo*zo + l_k*l_k) / cp;
    tau_4 = std::sqrt(zo*zo + l_k*l_k + s_k*s_k) / cp;

    for (it=0; it<nt; it++) {

      t = (((double) it) * dt + delay)/1.0e6; // in [s].

      a_k = 0;
      if ( (t >= tau_1) && (t <= tau_4) )
        a_k += cp/4;

      if ( (t >= tau_2) && (t <= tau_4) )
        a_k -= cp/(2*M_PI) * std::acos( s_k / (cp * std::sqrt(t*t - t_z*t_z)));

      if ( (t >= tau_3) && (t <= tau_4) )
        a_k -= cp/(2*M_PI) * std::acos( l_k / (cp * std::sqrt(t*t - t_z*t_z)));

      h[it] += g_k * a_k;

    } // for it
  } // for k

  return;
}
