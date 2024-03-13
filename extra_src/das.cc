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

#ifndef USE_OPENCL
template <class T>
DAS<T>::DAS(DASType das_type,
            dream_idx_type a_scan_len,
            dream_idx_type No,
            dream_idx_type num_t_elements,
            dream_idx_type num_r_elements, // SAFT if num_r_elements = 0;
            bool use_gpu)
  : m_out_err(ErrorLevel::none)
  , m_das_type(das_type)
  , m_a_scan_len(a_scan_len)
  , m_No(No)
  , m_num_t_elements(num_t_elements)
  , m_num_r_elements(num_r_elements)
{
  // Here we have not build with OpenCL support so switch off GPU support.
  m_use_gpu = false;
}
#endif

template <class T>
void DAS<T>::set_running()
{
  running = true;
}

template <class T>
void DAS<T>::abort(int signum)
{
  running = false;
}

template <class T>
bool DAS<T>::is_running()
{
  return running;
}

// Thread data.
template <typename T>
struct DATA
{
  const T *Y;
  dream_idx_type a_scan_len;
  dream_idx_type start;
  dream_idx_type stop;
  const T *Ro;
  dream_idx_type No;
  const T *Gt;                   // Transmit
  dream_idx_type num_t_elements;
  const T *Gr;                   // Recieve
  dream_idx_type num_r_elements;
  T dt;
  DelayType delay_type;
  const T *delay;
  T cp;
  T *Im;
  DASType das_type;
  ErrorLevel err_level;
};

/***
 *
 * Thread function.
 *
 ***/

template <class T>
void* DAS<T>::smp_das(void *arg)
{
  DATA<T> D = *(DATA<T> *)arg;
  const T *Y = D.Y;
  dream_idx_type a_scan_len=D.a_scan_len;
  const T *Gt = D.Gt;
  dream_idx_type num_t_elements = D.num_t_elements;
  const T *Gr = D.Gr;
  dream_idx_type num_r_elements  = D.num_r_elements;
  T dt=D.dt;
  dream_idx_type No=D.No;
  const T *delay=D.delay;
  const T *Ro=D.Ro;
  T cp=D.cp;
  DASType das_type=D.das_type;
  T *Im = D.Im;
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
    T xo = Ro[no];
    T yo = Ro[no+1*No];
    T zo = Ro[no+2*No];

    T dlay = 0.0;
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

    if (das_type == DASType::rca_coltx) {
      err = das_rca_serial_coltx(Y, a_scan_len,
                                 Gt, num_t_elements, // Gt = G_col
                                 Gr, num_r_elements, // Gr = G_row
                                 xo, yo, zo,
                                 dt, dlay,
                                 cp, Im[no], tmp_lev);
    }

    if (das_type == DASType::rca_rowtx) {
      err = das_rca_serial_rowtx(Y, a_scan_len,
                                 Gt, num_t_elements, // Gt = G_col
                                 Gr, num_r_elements, // Gr = G_row
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

template <class T>
ErrorLevel DAS<T>::das(const T *Y, const T *Ro, const T *Gt, const T *Gr,
                       T dt,
                       DelayType delay_type, const T *delay,
                       T cp,
                       T *Im,
                       ErrorLevel err_level)
{
  // Force SAFT if Gt is empty.
  DASType das_type = m_das_type;
  if (m_num_r_elements == 0) {
    das_type = DASType::saft;
  }

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
  if (nthreads > m_No) {
    nthreads = m_No;
  }

  // Allocate local data.
  DATA<T> *D = (DATA<T>*) malloc(nthreads*sizeof(DATA<T>));

  // Allocate mem for the threads.
  std::thread *threads = new std::thread[nthreads]; // Init thread data.

  for (dream_idx_type thread_n = 0; thread_n < nthreads; thread_n++) {

    dream_idx_type start = thread_n * m_No/nthreads;
    dream_idx_type stop =  (thread_n+1) * m_No/nthreads;

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
    for (dream_idx_type thread_n = 0; thread_n < nthreads; thread_n++) {
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
 *  das-saft - Delay-and-sum synthetic aperture focusing method (SAFT)
 *
 ***/

template <class T>
ErrorLevel DAS<T>::das_saft_serial(const T *Y, // Size: a_scan_len x num_elements
                                   dream_idx_type a_scan_len,
                                   const T *g, dream_idx_type num_elements,
                                   T xo, T yo, T zo,
                                   T dt, T delay,
                                   T cp, T &im, ErrorLevel err_level)
{
  ErrorLevel err = ErrorLevel::none;
  const T Fs_khz = (1.0/dt)*1000.0;
  const T one_over_cp = 1.0/cp;
  const T delay_ms = delay/1000.0;

  //
  //  Transmit with one element - receive with one element.
  //

  im = 0.0;
  const T *y_p = Y;
  for (dream_idx_type n_tr=0; n_tr<num_elements; n_tr++) {

    // Transmit/Receive
    T gx_tr = g[n_tr] - xo;
    T gy_tr = g[n_tr + 1*num_elements] - yo;
    T gz_tr = g[n_tr + 2*num_elements] - zo;
    T t_tr =  std::sqrt(gx_tr*gx_tr + gy_tr*gy_tr + gz_tr*gz_tr) * one_over_cp;
    t_tr += delay_ms;           // Compensate for pulse system delay.

    T t_dp = 2.0*t_tr; // Double-path travel time.
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

template <class T>
ErrorLevel DAS<T>::das_tfm_serial(const T *Y, // Size: a_scan_len x num_t_elements*num_r_elements (=FMC)
                                  dream_idx_type a_scan_len,
                                  const T *Gt, dream_idx_type num_t_elements,
                                  const T *Gr, dream_idx_type num_r_elements,
                                  T xo, T yo, T zo,
                                  T dt, T delay,
                                  T cp, T &im, ErrorLevel err_level)
{
  ErrorLevel err = ErrorLevel::none;
  const T Fs_khz = (1.0/dt)*1000.0;
  const T one_over_cp = 1.0/cp;
  const T delay_ms = delay/1000.0;

  //
  //  Transmit with all - receive with all elements.
  //

  im = 0.0;
  const T *y_p = Y;
  for (dream_idx_type n_t=0; n_t<num_t_elements; n_t++) {

    // Transmit
    T gx_t = Gt[n_t] - xo;
    T gy_t = Gt[n_t + 1*num_t_elements] - yo;
    T gz_t = Gt[n_t + 2*num_t_elements] - zo;
    T t_t =  std::sqrt(gx_t*gx_t + gy_t*gy_t + gz_t*gz_t) * one_over_cp;
    t_t += delay_ms;            // Compensate for pulse system delay.

    for (dream_idx_type n_r=0; n_r<num_r_elements; n_r++) {

      // Recieve
      T gx_r = Gr[n_r] - xo;
      T gy_r = Gr[n_r + 1*num_r_elements] - yo;
      T gz_r = Gr[n_r + 2*num_r_elements] - zo;
      T t_r =  std::sqrt(gx_r*gx_r + gy_r*gy_r + gz_r*gz_r) * one_over_cp;

      T t_dp = t_t + t_r; // Double-path travel time.
      auto k = dream_idx_type(t_dp*Fs_khz);

      //if (n_r == n_t) { // For comparing with SAFT (for testing)
      if ((k < a_scan_len) && (k >= 0)) {
        im += y_p[k];
      } else {
        if (k >= 0) {
          err = dream_out_of_bounds_err("DAS out of bounds +",k-a_scan_len+1,err_level);
        } else {
          err = dream_out_of_bounds_err("DAS out of bounds -",k,err_level);
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
 *  das-rca - Delay-and-sum for the row-column addressed 2D array variant of TFM.
 *
 ***/

// Transmit with columns - receive with rows.
template <class T>
ErrorLevel DAS<T>::das_rca_serial_coltx(const T *Y, // Size: a_scan_len x num_cols*num_rows (=FMC)
                                        dream_idx_type a_scan_len,
                                        const T *G_col, dream_idx_type num_cols,
                                        const T *G_row, dream_idx_type num_rows,
                                        T xo, T yo, T zo,
                                        T dt, T delay,
                                        T cp, T &im, ErrorLevel err_level)
{
  ErrorLevel err = ErrorLevel::none;
  const T Fs_khz = (1.0/dt)*1000.0;
  const T one_over_cp = 1.0/cp;
  const T delay_ms = delay/1000.0;

  // Here we have crossed striped electrodes and the length of the transmit element stripe
  // is determined by the two edge positions of the receive elements and vice versa.
  //
  // We assume:
  //
  // * the transmit elements are distributed along the x-dimension (transmit with columns)
  //   and the receive elements are distributed along the y-dimension (rows),
  // * G_col[0] < G_col[num_cols-1] (x-dim edge positions) and
  //   G_row[0 + 1*num_cols] < G_col[num_cols-1 + 1*num_cols]
  //   (y-dim edge positions).

  // Element (stripe) lengths
  const T gc_y_min = G_row[0 + 1*num_cols];
  const T gc_y_max = G_row[num_cols-1 + 1*num_cols];
  const T gr_x_min = G_col[0];
  const T gr_x_max = G_col[num_rows-1];

  im = 0.0;
  const T *y_p = Y;

  for (dream_idx_type n_c=0; n_c<num_cols; n_c++) {

    // Transmit
    T gx_c = G_col[n_c] - xo;
    T gy_c;
    if ( (yo >= gc_y_min) && yo <= gc_y_max) {
      gy_c = 0.0;               // We are inside the stripe aperture.
    } else {    // Here we use the distance to the edge of the stripe.
      if (yo < gc_y_min) {
        gy_c = gc_y_min - yo;
      } else { // yo > gc_y_max
        gy_c = gc_y_max - yo;
      }
    }
    T gz_c = G_col[n_c + 2*num_cols] - zo;
    T t_c =  std::sqrt(gx_c*gx_c + gy_c*gy_c + gz_c*gz_c) * one_over_cp;
    t_c += delay_ms;            // Compensate for system pulse delay.

    for (dream_idx_type n_r=0; n_r<num_rows; n_r++) {

      // Recieve
      T gx_r;
      if ( (xo >= gr_x_min) && xo <= gr_x_max) {
        gx_r = 0.0;               // We are inside the stripe aperture.
      } else {    // Here we use the distance to the edge of the stripe.
        if (xo < gr_x_min) {
          gx_r = gr_x_min - xo;
        } else { // yo > gc_y_max
          gx_r = gr_x_max - xo;
        }
      }
      T gy_r = G_row[n_r + 1*num_rows] - yo;
      T gz_r = G_row[n_r + 2*num_rows] - zo;
      T t_r =  std::sqrt(gx_r*gx_r + gy_r*gy_r + gz_r*gz_r) * one_over_cp;

      T t_dp = t_c + t_r; // Double-path travel time.
      auto k = dream_idx_type(t_dp*Fs_khz);

      if ((k < a_scan_len) && (k >= 0)) {
        im += y_p[k];
      } else {
        if (k >= 0) {
          err = dream_out_of_bounds_err("DAS out of bounds +",k-a_scan_len+1,err_level);
        } else {
          err = dream_out_of_bounds_err("DAS out of bounds -",k,err_level);
        }
      }
      y_p += a_scan_len; // Jump to the next A-scan
    }
  }

  return err;
};


// Transmit with rows - receive with cols.
template <class T>
ErrorLevel DAS<T>::das_rca_serial_rowtx(const T *Y, // Size: a_scan_len x num_cols*num_rows (=FMC)
                                        dream_idx_type a_scan_len,
                                        const T *G_col, dream_idx_type num_cols,
                                        const T *G_row, dream_idx_type num_rows,
                                        T xo, T yo, T zo,
                                        T dt, T delay,
                                        T cp, T &im, ErrorLevel err_level)
{
  ErrorLevel err = ErrorLevel::none;
  const T Fs_khz = (1.0/dt)*1000.0;
  const T one_over_cp = 1.0/cp;
  const T delay_ms = delay/1000.0;

  // Here we have crossed striped electrodes and the length of the transmit element stripe
  // is determined by the two edge positions of the receive elements and vice versa.
  //
  // We assume:
  //
  // * the transmit elements are distributed along the y-dimension (transmit with rows)
  //   and the receive elements are distributed along the x-dimension (cols),
  // * G_col[0] < G_col[num_cols-1] (x-dim edge positions) and
  //   G_row[0 + 1*num_cols] < G_col[num_cols-1 + 1*num_cols]
  //   (y-dim edge positions).

  // Element (stripe) lengths
  const T gc_y_min = G_row[0 + 1*num_cols];
  const T gc_y_max = G_row[num_cols-1 + 1*num_cols];
  const T gr_x_min = G_col[0];
  const T gr_x_max = G_col[num_rows-1];

  im = 0.0;
  const T *y_p = Y;

  for (dream_idx_type n_r=0; n_r<num_rows; n_r++) {

    // Rows
    T gx_r;
    if ( (xo >= gr_x_min) && xo <= gr_x_max) {
      gx_r = 0.0;               // We are inside the stripe aperture.
    } else {    // Here we use the distance to the edge of the stripe.
      if (xo < gr_x_min) {
        gx_r = gr_x_min - xo;
      } else { // yo > gc_y_max
        gx_r = gr_x_max - xo;
      }
    }
    T gy_r = G_row[n_r + 1*num_rows] - yo;
    T gz_r = G_row[n_r + 2*num_rows] - zo;
    T t_r =  std::sqrt(gx_r*gx_r + gy_r*gy_r + gz_r*gz_r) * one_over_cp;

    for (dream_idx_type n_c=0; n_c<num_cols; n_c++) {

      // Cols
      T gx_c = G_col[n_c] - xo;
      T gy_c;
      if ( (yo >= gc_y_min) && yo <= gc_y_max) {
        gy_c = 0.0;               // We are inside the stripe aperture.
      } else {    // Here we use the distance to the edge of the stripe.
        if (yo < gc_y_min) {
          gy_c = gc_y_min - yo;
        } else { // yo > gc_y_max
          gy_c = gc_y_max - yo;
        }
      }

      T gz_c = G_col[n_c + 2*num_cols] - zo;
      T t_c =  std::sqrt(gx_c*gx_c + gy_c*gy_c + gz_c*gz_c) * one_over_cp;
      t_c += delay_ms;            // Compensate for system pulse delay.

      T t_dp = t_c + t_r; // Double-path travel time.
      auto k = dream_idx_type(t_dp*Fs_khz);

      if ((k < a_scan_len) && (k >= 0)) {
        im += y_p[k];
      } else {
        if (k >= 0) {
          err = dream_out_of_bounds_err("DAS out of bounds +",k-a_scan_len+1,err_level);
        } else {
          err = dream_out_of_bounds_err("DAS out of bounds -",k,err_level);
        }
      }
      y_p += a_scan_len; // Jump to the next A-scan
    }
  }

  return err;
};


// https://bytefreaks.net/programming-2/c/c-undefined-reference-to-templated-class-function

// Support doubles and floats
template class DAS<double>;
template class DAS<float>;
