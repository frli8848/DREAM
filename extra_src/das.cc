/***
*
* Copyright (C) 2003,2006,2007,2008,2009,2014,2021 Fredrik Lingvall
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

#include "das.h"
#include "affinity.h"

std::mutex err_mutex;
std::atomic<bool> running;

void DAS::abort(int signum)
{
  running = false;
}

bool DAS::is_running()
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
  double *gt;                   // Transmit
  dream_idx_type num_t_elements;
  double *gr;                   // Recieve
  dream_idx_type num_r_elements;
  double dt;
  DelayType delay_type;
  double *delay;
  double cp;
  double *Y;
  dream_idx_type a_scan_len;
  double *Im;
  DASType das_type;
  ErrorLevel err_level;
} DATA;

/***
 *
 * Thread function.
 *
 ***/

void* DAS::smp_das(void *arg)
{
  DATA D = *(DATA *)arg;
  double *Y = D.Y;
  dream_idx_type a_scan_len=D.a_scan_len;
  double *Im = D.Im;
  double *gt = D.gt;
  dream_idx_type num_t_elements = D.num_t_elements;
  double *gr = D.gr;
  dream_idx_type num_r_elements  = D.num_r_elements;
  double dt=D.dt;
  dream_idx_type n, no=D.no;
  double *delay=D.delay, *ro=D.ro, cp=D.cp;
  DASType das_type=D.das_type;
  dream_idx_type start=D.start, stop=D.stop;
  ErrorLevel tmp_lev=ErrorLevel::none, err_level=D.err_level;
  ErrorLevel tmp_err=ErrorLevel::none, err=ErrorLevel::none;

  // Let the thread finish and then catch the error.
  if (err_level == ErrorLevel::stop) {
    tmp_lev = ErrorLevel::parallel_stop;
  } else {
    tmp_lev = err_level;
  }

  for (n=start; n<stop; n++) {
    double xo = ro[n];
    double yo = ro[n+1*no];
    double zo = ro[n+2*no];

    double dlay = 0.0;
    if (D.delay_type == DelayType::single) {
      dlay = delay[0];
    } else { // DelayType::multiple.
      dlay = delay[n];
    }

    if (das_type == DASType::saft) {
      err = das_saft_serial(Y, a_scan_len,
                            gt, num_t_elements, // Transmit = receive here.
                            xo, yo, zo,
                            dt, dlay,
                            cp, Im[n], tmp_lev);
    }

    if (das_type == DASType::tfm) {
      err = das_tfm_serial(Y, a_scan_len,
                           gt, num_t_elements,
                           gr, num_r_elements,
                           xo, yo, zo,
                           dt, dlay,
                           cp, Im[n], tmp_lev);
    }

    if (err != ErrorLevel::none || m_out_err ==  ErrorLevel::parallel_stop) {
      tmp_err = err;
      if (err == ErrorLevel::parallel_stop || m_out_err ==  ErrorLevel::parallel_stop) {
        break; // Jump out when an ErrorLevel::stop error occurs.
      }
    }

    if (!running) {
      std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!\n";
      return(NULL);
    }
  }

  if ((tmp_err != ErrorLevel::none) && (m_out_err == ErrorLevel::none)) {
    std::lock_guard<std::mutex> lk(err_mutex);

    m_out_err = tmp_err;
  }

  return(NULL);
}

/***
 *
 * Delay-and-sum (DAS) processing using SAFT or TFM methods.
 *
 ***/

ErrorLevel DAS::das(double *Y, dream_idx_type a_scan_len,
                    double *ro, dream_idx_type no,
                    double *gt, dream_idx_type num_t_elements,
                    double *gr, dream_idx_type num_r_elements, // SAFT if num_r_elements = 0;
                    double dt,
                    DelayType delay_type, double *delay,
                    double cp,
                    double *Im,
                    ErrorLevel err_level)
{
  std::thread *threads;
  unsigned int thread_n, nthreads;
  dream_idx_type start, stop;
  DATA *D;

  DASType das_type = DASType::saft;
  if (num_r_elements > 0) {
    das_type = DASType::tfm;
  }

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
    D[thread_n].Y = Y;
    D[thread_n].a_scan_len = a_scan_len;
    D[thread_n].Im = Im;
    D[thread_n].no = no;
    D[thread_n].ro = ro;
    D[thread_n].gt = gt;
    D[thread_n].num_t_elements = num_t_elements;
    D[thread_n].gr = gr;
    D[thread_n].num_r_elements = num_r_elements;
    D[thread_n].dt = dt;
    D[thread_n].delay_type = delay_type;
    D[thread_n].delay = delay;
    D[thread_n].cp = cp;
    D[thread_n].das_type = das_type;
    D[thread_n].err_level = err_level;

    if (nthreads>1) {
      // Start the threads.
      threads[thread_n] = das_thread(&D[thread_n]);
      set_dream_thread_affinity(thread_n, nthreads, threads);
    } else {
      smp_das(&D[0]);
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
 *  das-saft - Delay-and-sum synthetic focusing method (SAFT)
 *
 ***/

ErrorLevel DAS::das_saft_serial(double *Y, // Size: a_scan_len x num_elements
                                dream_idx_type a_scan_len,
                                double *g, dream_idx_type num_elements,
                                double xo, double yo, double zo,
                                double dt, double delay,
                                double cp, double &im, ErrorLevel err_level)
{
  ErrorLevel err = ErrorLevel::none;
  const double Fs_khz = (1.0/dt)*1000.0;
  double one_over_cp = 1.0/cp;

  //
  //  Transmit with one element - receive with one element.
  //

  im = 0.0;
  double *y_p = Y;
  for (size_t n_tr=0; n_tr<num_elements; n_tr++) {

    // Transmit/Receive
    double gx_tr = g[n_tr] - xo;
    double gy_tr = g[n_tr + 1*num_elements] - yo;
    double gz_tr = g[n_tr + 2*num_elements] - zo;
    double t_tr =  std::sqrt(gx_tr*gx_tr + gy_tr*gy_tr + gz_tr*gz_tr) * one_over_cp;

    double t_dp = 2.0*t_tr; // Double-path travel time.
    auto k = dream_idx_type(t_dp*Fs_khz - delay);

    if ((k < a_scan_len) && (k >= 0)) {
      im += y_p[k];
    } else {
      if (k >= 0) {
        err = dream_out_of_bounds_err("DAS out of bounds",k-a_scan_len+1,err_level);
      } else {
        err = dream_out_of_bounds_err("DAS out of bounds",k,err_level);
      }
    }

    y_p += a_scan_len; // Jump to the next A-scan
  }

  return err;
};

/***
 *
 *  das-tfm - Delay-and-sum total focusing method (TFM)
 *
 ***/

ErrorLevel DAS::das_tfm_serial(double *Y, // Size: a_scan_len x num_t_elements*num_r_elements (=FMC)
                               dream_idx_type a_scan_len,
                               double *gt, dream_idx_type num_t_elements,
                               double *gr, dream_idx_type num_r_elements,
                               double xo, double yo, double zo,
                               double dt, double delay,
                               double cp, double &im, ErrorLevel err_level)
{
  ErrorLevel err = ErrorLevel::none;
  const double Fs_khz = (1.0/dt)*1000.0;
  double one_over_cp = 1.0/cp;

  //
  //  Transmit with all - receive with all elements.
  //

  im = 0.0;
  double *y_p = Y;
  for (size_t n_t=0; n_t<num_t_elements; n_t++) {

    // Transmit
    double gx_t = gt[n_t] - xo;
    double gy_t = gt[n_t + 1*num_t_elements] - yo;
    double gz_t = gt[n_t + 2*num_t_elements] - zo;
    double t_t =  std::sqrt(gx_t*gx_t + gy_t*gy_t + gz_t*gz_t) * one_over_cp;

    for (size_t n_r=0; n_r<num_r_elements; n_r++) {

      // Recieve
      double gx_r = gr[n_r] - xo;
      double gy_r = gr[n_r + 1*num_r_elements] - yo;
      double gz_r = gr[n_r + 2*num_r_elements] - zo;
      double t_r =  std::sqrt(gx_r*gx_r + gy_r*gy_r + gz_r*gz_r) * one_over_cp;

      double t_dp = t_t + t_r; // Double-path travel time.
      auto k = dream_idx_type(t_dp*Fs_khz - delay);

      //if (n_r == n_t) { // For comparing with SAFT (for testing)
        if ((k < a_scan_len) && (k >= 0)) {
          im += y_p[k];
        } else {
          if (k >= 0) {
            err = dream_out_of_bounds_err("DAS out of bounds",k-a_scan_len+1,err_level);
          } else {
            err = dream_out_of_bounds_err("DAS out of bounds",k,err_level);
          }
        }
        //} // SAFT testing
      y_p += a_scan_len; // Jump to the next A-scan
    }
  }

  return err;
};
