/***
 *
 * Copyright (C) 2003,2006,2007,2008,2009,2014,2021,2023 Fredrik Lingvall
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

void DAS::set_running()
{
  running = true;
}

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
  dream_idx_type start;
  dream_idx_type stop;
  double *Ro;
  dream_idx_type No;
  double *Gt;                   // Transmit
  dream_idx_type num_t_elements;
  double *Gr;                   // Recieve
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
  double *Gt = D.Gt;
  dream_idx_type num_t_elements = D.num_t_elements;
  double *Gr = D.Gr;
  dream_idx_type num_r_elements  = D.num_r_elements;
  double dt=D.dt;
  dream_idx_type No=D.No;
  double *delay=D.delay, *Ro=D.Ro, cp=D.cp;
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

  for (dream_idx_type no=start; no<stop; no++) {
    double xo = Ro[no];
    double yo = Ro[no+1*No];
    double zo = Ro[no+2*No];

    double dlay = 0.0;
    if (D.delay_type == DelayType::single) {
      dlay = delay[0];
    } else { // DelayType::multiple.
      dlay = delay[no];
    }

    if (das_type == DASType::saft) {
      err = das_saft_serial(Y, a_scan_len,
                            Gt, num_t_elements, // Transmit = receive here.
                            xo, yo, zo,
                            dt, dlay,
                            cp, Im[no], tmp_lev);
    }

    if (das_type == DASType::tfm) {
      err = das_tfm_serial(Y, a_scan_len,
                           Gt, num_t_elements,
                           Gr, num_r_elements,
                           xo, yo, zo,
                           dt, dlay,
                           cp, Im[no], tmp_lev);
    }

    if (das_type == DASType::rca) {
      err = das_rca_serial(Y, a_scan_len,
                           Gt, num_t_elements,
                           Gr, num_r_elements,
                           xo, yo, zo,
                           dt, dlay,
                           cp, Im[no], tmp_lev);
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
 * Delay-and-sum (DAS) processing using SAFT or TFM ErrorLevel.
 *
 ***/

ErrorLevel DAS::das(double *Y, double *Ro, double *Gt, double *Gr,
                    double dt,
                    DelayType delay_type, double *delay,
                    double cp,
                    double *Im,
                    ErrorLevel err_level)
{
  std::thread *threads;
  dream_idx_type thread_n, nthreads;
  dream_idx_type start, stop;
  DATA *D;

  // Force SAFT if Gt is empty.
  DASType das_type = m_das_type;
  if (m_num_r_elements == 0) {
    das_type = DASType::saft;
  }

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
  if (nthreads > m_No) {
    nthreads = m_No;
  }

  // Allocate local data.
  D = (DATA*) malloc(nthreads*sizeof(DATA));

  // Allocate mem for the threads.
  threads = new std::thread[nthreads]; // Init thread data.

  for (thread_n = 0; thread_n < nthreads; thread_n++) {

    start = thread_n * m_No/nthreads;
    stop =  (thread_n+1) * m_No/nthreads;

    // Init local data.
    D[thread_n].start = start; // Local start index;
    D[thread_n].stop = stop; // Local stop index;
    D[thread_n].Y = Y;
    D[thread_n].a_scan_len = m_a_scan_len;
    D[thread_n].Im = Im;
    D[thread_n].No = m_No;
    D[thread_n].Ro = Ro;
    D[thread_n].Gt = Gt;
    D[thread_n].num_t_elements = m_num_t_elements;
    D[thread_n].Gr = Gr;
    D[thread_n].num_r_elements = m_num_r_elements;
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
  const double one_over_cp = 1.0/cp;
  const double delay_ms = delay/1000.0;

  //
  //  Transmit with one element - receive with one element.
  //

  im = 0.0;
  double *y_p = Y;
  for (dream_idx_type n_tr=0; n_tr<num_elements; n_tr++) {

    // Transmit/Receive
    double gx_tr = g[n_tr] - xo;
    double gy_tr = g[n_tr + 1*num_elements] - yo;
    double gz_tr = g[n_tr + 2*num_elements] - zo;
    double t_tr =  std::sqrt(gx_tr*gx_tr + gy_tr*gy_tr + gz_tr*gz_tr) * one_over_cp;
    t_tr += delay_ms;           // Compensate for pulse system delay.

    double t_dp = 2.0*t_tr; // Double-path travel time.
    auto k = dream_idx_type(t_dp*Fs_khz);

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
                               double *Gt, dream_idx_type num_t_elements,
                               double *Gr, dream_idx_type num_r_elements,
                               double xo, double yo, double zo,
                               double dt, double delay,
                               double cp, double &im, ErrorLevel err_level)
{
  ErrorLevel err = ErrorLevel::none;
  const double Fs_khz = (1.0/dt)*1000.0;
  const double one_over_cp = 1.0/cp;
  const double delay_ms = delay/1000.0;

  //
  //  Transmit with all - receive with all elements.
  //

  im = 0.0;
  double *y_p = Y;
  for (dream_idx_type n_t=0; n_t<num_t_elements; n_t++) {

    // Transmit
    double gx_t = Gt[n_t] - xo;
    double gy_t = Gt[n_t + 1*num_t_elements] - yo;
    double gz_t = Gt[n_t + 2*num_t_elements] - zo;
    double t_t =  std::sqrt(gx_t*gx_t + gy_t*gy_t + gz_t*gz_t) * one_over_cp;
    t_t += delay_ms;            // Compensate for pulse system delay.

    for (dream_idx_type n_r=0; n_r<num_r_elements; n_r++) {

      // Recieve
      double gx_r = Gr[n_r] - xo;
      double gy_r = Gr[n_r + 1*num_r_elements] - yo;
      double gz_r = Gr[n_r + 2*num_r_elements] - zo;
      double t_r =  std::sqrt(gx_r*gx_r + gy_r*gy_r + gz_r*gz_r) * one_over_cp;

      double t_dp = t_t + t_r; // Double-path travel time.
      auto k = dream_idx_type(t_dp*Fs_khz);

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


/***
 *
 *  das-rca - Delay-and-sum for the row-column addressed 2D array variant of TFM
 *
 ***/

ErrorLevel DAS::das_rca_serial(double *Y, // Size: a_scan_len x num_t_elements*num_r_elements (=FMC)
                               dream_idx_type a_scan_len,
                               double *Gt, dream_idx_type num_t_elements,
                               double *Gr, dream_idx_type num_r_elements,
                               double xo, double yo, double zo,
                               double dt, double delay,
                               double cp, double &im, ErrorLevel err_level)
{
  ErrorLevel err = ErrorLevel::none;
  const double Fs_khz = (1.0/dt)*1000.0;
  const double one_over_cp = 1.0/cp;
  const double delay_ms = delay/1000.0;

  // Here we have crossed striped electrodes and the length of the transmit element stripe
  // is determined by the two edge positions of the receive elements and vice versa.
  //
  // We assume:
  //
  // * the transmit elements are distributed along the x-dimension and
  //   the receive elements are distributed along the y-dimension,
  // * Gt[0] < Gt[num_t_elements-1] (x-dim edge positions) and
  //   Gr[0 + 1*num_t_elements] < Gt[num_t_elements-1 + 1*num_t_elements]
  //   (y-dim edge positions).

  // Element (stripe) lengths
  const double gt_y_min = Gr[0 + 1*num_t_elements];
  const double gt_y_max = Gr[num_t_elements-1 + 1*num_t_elements];
  const double gr_x_min = Gt[0];
  const double gr_x_max = Gt[num_r_elements-1];

  //
  //  Transmit with all - receive with all elements.
  //

  im = 0.0;
  double *y_p = Y;
  for (dream_idx_type n_t=0; n_t<num_t_elements; n_t++) {

    // Transmit
    double gx_t = Gt[n_t] - xo;
    double gy_t;
    if ( (yo >= gt_y_min) && yo <= gt_y_max) {
      gy_t = 0.0;               // We are inside the stripe aperture.
    } else {    // Here we use the distance to the edge of the stripe.
      if (yo < gt_y_min) {
        gy_t = gt_y_min - yo;
      } else { // yo > gt_y_max
        gy_t = gt_y_max - yo;
      }
    }
    double gz_t = Gt[n_t + 2*num_t_elements] - zo;
    double t_t =  std::sqrt(gx_t*gx_t + gy_t*gy_t + gz_t*gz_t) * one_over_cp;
    t_t += delay_ms;            // Compensate for system pulse delay.

    for (dream_idx_type n_r=0; n_r<num_r_elements; n_r++) {

      // Recieve
      double gx_r;
      if ( (yo >= gr_x_min) && yo <= gr_x_max) {
        gx_r = 0.0;               // We are inside the stripe aperture.
      } else {    // Here we use the distance to the edge of the stripe.
        if (yo < gr_x_min) {
          gx_r = gr_x_min - xo;
        } else { // yo > gt_y_max
          gx_r = gr_x_max - xo;
        }
      }
      double gy_r = Gr[n_r + 1*num_r_elements] - yo;
      double gz_r = Gr[n_r + 2*num_r_elements] - zo;
      double t_r =  std::sqrt(gx_r*gx_r + gy_r*gy_r + gz_r*gz_r) * one_over_cp;

      double t_dp = t_t + t_r; // Double-path travel time.
      auto k = dream_idx_type(t_dp*Fs_khz);

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
  }

  return err;
};
