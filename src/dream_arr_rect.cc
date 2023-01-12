/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2014,2019,2021 Fredrik Lingvall
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

#include "dream_arr_rect.h"
#include "dreamrect.h"
#include "arr_functions.h"
#include "affinity.h"

std::mutex err_mutex;
std::atomic<bool> running;

void ArrRect::abort(int signum)
{
  running = false;
}

bool ArrRect::is_running()
{
  return running;
}

/***
 *
 * Thread function.
 *
 ***/

void* ArrRect::smp_dream_arr_rect(void *arg)
{
  ErrorLevel tmp_err=ErrorLevel::none, err=ErrorLevel::none;
  DATA D = *(DATA *)arg;
  double xo, yo, zo;
  double *h = D.h;
  double a=D.a, b=D.b, dx=D.dx, dy=D.dy, dt=D.dt;
  dream_idx_type n, no=D.no, nt=D.nt;
  ErrorLevel tmp_lev=ErrorLevel::none, err_level=D.err_level;
  double *delay=D.delay, *ro=D.ro, v=D.v, cp=D.cp;
  Attenuation *att = D.att;
  dream_idx_type start=D.start, stop=D.stop;
  FocusMet foc_met=D.foc_met;
  SteerMet steer_met=D.steer_met;
  int do_apod = D.do_apod;
  ApodMet apod_met = D.apod_met;
  double *focal=D.focal, *apod=D.apod, theta=D.theta,phi=D.phi,apod_par=D.apod_par;
  dream_idx_type num_elements = D.num_elements;

  double *gx = D.G;               // First column in the matrix.
  double *gy = gx + num_elements; // Second column in the matrix.
  double *gz = gy + num_elements; // Third column in the matrix.

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
      err = dream_arr_rect_serial(xo, yo, zo,
                                  a, b,
                                  dx, dy, dt, nt,
                                  dlay, v, cp,
                                  num_elements, gx, gy, gz,
                                  foc_met, focal,
                                  steer_met, theta, phi,
                                  apod, do_apod, apod_met, apod_par,
                                  &h[n*nt], tmp_lev);

    } else {
      err = dream_arr_rect_serial(*att, *xc_vec, *x_vec,
                                  xo, yo, zo,
                                  a, b,
                                  dx, dy, dt, nt,
                                  dlay, v, cp,
                                  num_elements, gx, gy, gz,
                                  foc_met, focal,
                                  steer_met, theta, phi,
                                  apod, do_apod, apod_met, apod_par,
                                  &h[n*nt], tmp_lev);
    }

    if (err != ErrorLevel::none || m_out_err ==  ErrorLevel::parallel_stop) {
      tmp_err = err;
      if (err == ErrorLevel::parallel_stop || m_out_err ==  ErrorLevel::parallel_stop) {
        break; // Jump out when a ErrorLevel::stop error occurs.
      }
    }

    if (!running) {
      std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!" << std::endl;
      return(NULL);
    }

  }

  if ((tmp_err != ErrorLevel::none) && (m_out_err == ErrorLevel::none)) {
    std::lock_guard<std::mutex> lk(err_mutex);

    m_out_err = tmp_err;
  }

  return(NULL);
}

ErrorLevel ArrRect::dream_arr_rect_serial(double xo, double yo, double zo,
                                 double a, double b,
                                 double dx, double dy, double dt,
                                 dream_idx_type nt,
                                 double delay, double v, double cp,
                                 dream_idx_type num_elements, double *gx, double *gy, double *gz,
                                 FocusMet foc_met, double *focal,
                                 SteerMet steer_met, double theta, double phi,
                                 double *apod, bool do_apod, ApodMet apod_met, double apod_par,
                                 double *h, ErrorLevel err_level)
{
  ErrorLevel err = ErrorLevel::none, out_err = ErrorLevel::none;

  double r_max, x_max, y_max;
  max_dim_arr(&x_max, &y_max, &r_max, gx, gy, gz, num_elements);

  for (dream_idx_type n=0; n<num_elements; n++) {

    double foc_delay = 0.0;
    if (foc_met != FocusMet::ud) {
      focusing(foc_met, focal[0], gx[n], gy[n], x_max, y_max, r_max, cp, &foc_delay);
    } else {
      focusing(foc_met, focal[n], gx[n], gy[n], x_max, y_max, r_max, cp, &foc_delay);
    }

    double steer_delay = 0.0;
    beamsteering(steer_met, theta, phi, gx[n], gy[n], x_max, y_max, r_max, cp, &steer_delay);

    double weight = 1.0;
    if (do_apod) {
      apodization(apod_met, n, apod, &weight, gx[n], gy[n], r_max, apod_par);
    }

    // Compute the response for the n:th elemen and add it to the impulse response vector h.
    err = dreamrect(xo - gx[n], yo - gy[n], zo - gz[n],
                    a, b, dx, dy, dt, nt,
                    delay - foc_delay - steer_delay,
                    v, cp,
                    h,
                    err_level,
                    weight);

    if (err != ErrorLevel::none) {
      out_err = err;
    }
  }

  return out_err;
}

ErrorLevel ArrRect::dream_arr_rect_serial(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                                 double xo, double yo, double zo,
                                 double a, double b,
                                 double dx, double dy, double dt,
                                 dream_idx_type nt, double delay, double v, double cp,
                                 dream_idx_type num_elements, double *gx, double *gy, double *gz,
                                 FocusMet foc_met, double *focal,
                                 SteerMet steer_met, double theta, double phi,
                                 double *apod, bool do_apod, ApodMet apod_met, double apod_par,
                                 double *h, ErrorLevel err_level)
{
  ErrorLevel err = ErrorLevel::none, out_err = ErrorLevel::none;

  for (dream_idx_type i=0; i<nt; i++) {
    h[i] = 0.0;
  }


  double r_max, x_max, y_max;
  max_dim_arr(&x_max, &y_max, &r_max, gx, gy, gz, num_elements);

  for (dream_idx_type n=0; n<num_elements; n++) {

    double foc_delay = 0.0;
    if (foc_met != FocusMet::ud) {
      focusing(foc_met, focal[0], gx[n], gy[n], x_max, y_max, r_max, cp, &foc_delay);
    } else {
      focusing(foc_met, focal[n], gx[n], gy[n], x_max, y_max, r_max, cp, &foc_delay);
    }

    double steer_delay = 0.0;
    beamsteering(steer_met, theta, phi, gx[n], gy[n], x_max, y_max, r_max, cp, &steer_delay);

    double weight = 1.0;
    if (do_apod) {
      apodization(apod_met, n, apod, &weight, gx[n], gy[n], r_max, apod_par);
    }

    // Compute the response for the n:th elemen and add it to the impulse response vector h.
    err = dreamrect(att, xc_vec, x_vec,
                    xo - gx[n], yo - gy[n], zo - gz[n],
                    a, b, dx, dy, dt, nt,
                    delay - foc_delay - steer_delay,
                    v, cp,
                    h,
                    err_level,
                    weight);
    if (err != ErrorLevel::none) {
      out_err = err;
    }
  }

  return out_err;
}

/***
 *
 * dream_arr_rect
 *
 ***/

ErrorLevel ArrRect::dream_arr_rect(double alpha,
                                   double *ro, dream_idx_type no,
                                   double a, double b,
                                   double dx, double dy, double dt,
                                   dream_idx_type nt,
                                   DelayType delay_type, double *delay,
                                   double v, double cp,
                                   dream_idx_type num_elements, double *G,
                                   FocusMet foc_met, double *focal,
                                   SteerMet steer_met, double theta, double phi,
                                   double *apod, bool do_apod, ApodMet apod_met, double apod_par,
                                   double *h, ErrorLevel err_level)
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

  Attenuation att(nt, dt, alpha);
  Attenuation *att_ptr = nullptr;
  if (alpha > std::numeric_limits<double>::epsilon() ) {
    att_ptr = &att;
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
    D[thread_n].dx = dx;
    D[thread_n].dy = dy;
    D[thread_n].dt = dt;
    D[thread_n].nt = nt;
    D[thread_n].delay_type = delay_type;
    D[thread_n].delay = delay;
    D[thread_n].v = v;
    D[thread_n].cp = cp;
    D[thread_n].att = att_ptr;
    D[thread_n].num_elements = num_elements;
    D[thread_n].G = G;
    D[thread_n].foc_met = foc_met;
    D[thread_n].steer_met = steer_met;
    D[thread_n].do_apod = do_apod;
    D[thread_n].apod_met = apod_met;
    D[thread_n].focal = focal;
    D[thread_n].apod = apod;
    D[thread_n].theta = theta;
    D[thread_n].phi = phi;
    D[thread_n].apod_par = apod_par;
    D[thread_n].h = h;
    D[thread_n].err_level = err_level;

    if (nthreads > 1) {
      // Start the threads.
      //threads[thread_n] = std::thread(smp_dream_arr_rect, &D[thread_n]); // Start the threads.
      threads[thread_n] = arr_rect_thread(&D[thread_n]); // Start the threads.
      set_dream_thread_affinity(thread_n, nthreads, threads);
    } else {
      smp_dream_arr_rect(&D[0]);
    }

  }

  // Wait for all threads to finish.
  if (nthreads > 1) {
    for (thread_n = 0; thread_n < nthreads; thread_n++) {
      threads[thread_n].join();
    }
  }

  // Free memory.
  free((void*) D);

  return  m_out_err;
}
