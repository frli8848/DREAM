/***
*
* Copyright (C) 2002,2003,2005,2006,2007,2008,2009,2014,2019,2021 Fredrik Lingvall
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

#include "dream_arr_annu.h"
#include "dreamcirc.h"
#include "affinity.h"

std::mutex err_mutex;
std::atomic<bool> running;

void ArrAnnu::abort(int signum)
{
  running = false;
}

bool ArrAnnu::is_running()
{
  return running;
}

/***
 *
 * Thread function.
 *
 ***/

void* ArrAnnu::smp_dream_arr_annu(void *arg)
{
  ErrorLevel tmp_err=ErrorLevel::none, err=ErrorLevel::none;
  DATA D = *(DATA *)arg;
  double xo, yo, zo;
  double *h = D.h;
  double dx=D.dx, dy=D.dy, dt=D.dt;
  dream_idx_type no=D.no, nt=D.nt;
  ErrorLevel tmp_lev=ErrorLevel::none, err_level=D.err_level;
  double *delay=D.delay, *ro=D.ro, v=D.v, cp=D.cp;
  Attenuation *att=D.att;
  dream_idx_type start=D.start, stop=D.stop;
  FocusMet foc_met=D.foc_met;
  bool do_apod=D.do_apod;
  ApodMet apod_met=D.apod_met;
  double *focal=D.focal, *apod=D.apod, apod_par=D.apod_par;
  dream_idx_type num_radii=D.num_radii;

  double *Gr=D.Gr;

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

  for (dream_idx_type n=start; n<stop; n++) {
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
      err = dream_arr_annu_serial(xo, yo, zo,
                                  dx, dy, dt, nt,
                                  dlay, v, cp,
                                  num_radii,Gr,
                                  foc_met, focal,
                                  apod, do_apod, apod_met, apod_par,
                                  &h[n*nt], tmp_lev);
    } else {
      err = dream_arr_annu_serial(*att, *xc_vec, *x_vec,
                                  xo, yo, zo,
                                  dx, dy, dt, nt,
                                  dlay, v, cp,
                                  num_radii, Gr,
                                  foc_met, focal,
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

/***
 *
 * Serial functions for computation of spatial impulse response of an annular array.
 *
 ***/

ErrorLevel ArrAnnu::dream_arr_annu_serial(double xo, double yo, double zo,
                                 double dx, double dy, double dt,
                                 dream_idx_type nt,
                                 double delay,
                                 double v, double cp,
                                 dream_idx_type num_radii, double *Gr,
                                 FocusMet foc_met, double *focal,
                                 double *apod, bool do_apod, ApodMet apod_met, double apod_par,
                                 double *h, ErrorLevel err_level)
{
  ErrorLevel err = ErrorLevel::none, out_err = ErrorLevel::none;

  // NB. The first ring have no inner radius (it is a normal
  // circular disc).
  dream_idx_type num_elements = (num_radii+1)/2; // The number of rings.

  // Allocate scratch data
  std::unique_ptr<double[]> h_disc = std::make_unique<double[]>(nt*num_radii); // Impulse responses for each disc.
  std::unique_ptr<double[]> h_ring = std::make_unique<double[]>(nt*num_elements); // Impulse responses for each ring.
  std::unique_ptr<double[]> ring_radii = std::make_unique<double[]>(num_elements); // Vector of ring radii.

  // FIXME: make_unique clears data

  // Clear output impulse response vector.
  for (dream_idx_type i=0; i<nt; i++) {
    h[i] = 0.0;
  }

  // Clear data
  for (dream_idx_type i=0; i<nt; i++) {
    for (dream_idx_type n=0; n<num_radii; n++) {
      h_disc[i + n*nt] = 0.0;
    }
  }

  // Compute the impulse reponses for all circular
  // (full) discs.
  for (dream_idx_type n=0; n<num_radii; n++) {

    if (n==0 && Gr[n] <= std::numeric_limits<double>::epsilon()  ) {
      for (dream_idx_type i=0; i<nt; i++) {
        h_disc[i] = 0.0; // h(i,1) = 0.0;
      }
    } else {
      err = dreamcirc(xo, yo, zo,
                      Gr[n],
                      dx, dy, dt, nt,
                      delay,
                      v, cp, &h_disc[n*nt], err_level);
      if (err != ErrorLevel::none) {
        out_err = err;
      }
    }
  }

  // Compute the "middle" radius for each array ring
  // of the transducer (which is the radius used when
  // focusing).
  annular_disc_radii(ring_radii.get(), Gr, num_radii);
  double ring_r_max = ring_radii[num_elements-1]; // Radius of the outer most (last) ring.

  for (dream_idx_type n=0; n<num_elements; n++) {

    double foc_delay = 0.0;
    if (foc_met != FocusMet::ud) {
      foc_delay = focusing_annular(foc_met, focal[0], ring_radii[n], ring_r_max, cp);
    } else {
      foc_delay = focusing_annular(foc_met, focal[n], ring_radii[n], ring_r_max, cp);
    }

    double weight = 1.0;
    if (do_apod){
      weight = apodization_annular(apod_met, n, apod, ring_radii[n], ring_r_max, apod_par);
    }

    resp_annular(h_disc.get(), h_ring.get(), nt, n); // Compute the impulse response for the n:th ring.
    superpos_annular(h_ring.get(), h, nt, weight, foc_delay, n, dt); // Add the n:th impulse response.
  }

  return out_err;
}

ErrorLevel ArrAnnu::dream_arr_annu_serial(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                                 double xo, double yo, double zo,
                                 double dx, double dy, double dt,
                                 dream_idx_type nt,
                                 double delay,
                                 double v, double cp,
                                 dream_idx_type num_radii, double *Gr,
                                 FocusMet foc_met, double *focal,
                                 double *apod, bool do_apod, ApodMet apod_met, double apod_par,
                                 double *h, ErrorLevel err_level)
{
  ErrorLevel err = ErrorLevel::none, out_err = ErrorLevel::none;

  // NB. The first ring have no inner radius (it is a normal
  // circular disc).
  dream_idx_type num_elements = (num_radii+1)/2; // The number of rings.

  // Allocate scratch data
  std::unique_ptr<double[]> h_disc = std::make_unique<double[]>(nt*num_radii); // Impulse responses for each disc.
  std::unique_ptr<double[]> h_ring = std::make_unique<double[]>(nt*num_elements); // Impulse responses for each ring.
  std::unique_ptr<double[]> ring_radii = std::make_unique<double[]>(num_elements); // Vector of ring radii.

  // Clear output impulse response vector.
  for (dream_idx_type i=0; i<nt; i++) {
    h[i] = 0.0;
  }

  // Clear data
  for (dream_idx_type i=0; i<nt; i++) {
    for (dream_idx_type n=0; n<num_radii; n++) {
      h_disc[i + n*nt] = 0.0;
    }
  }

  // Compute the impulse reponses for all circular
  // (full) discs.
  for (dream_idx_type n=0; n<num_radii; n++) {

    if (n==0 && Gr[n] <= std::numeric_limits<double>::epsilon()  ) {
      for (dream_idx_type i=0; i<nt; i++) {
        h_disc[i] = 0.0; // h(i,1) = 0.0;
      }
    } else {
      err = dreamcirc(att, xc_vec, x_vec,
                      xo, yo, zo,
                      Gr[n],
                      dx, dy, dt, nt,
                      delay,
                      v, cp, &h_disc[n*nt], err_level);
      if (err != ErrorLevel::none) {
        out_err = err;
      }
    }
  }

  // Compute the "middle" radius for each array ring
  // of the transducer (which is the radius used when
  // focusing).
  annular_disc_radii(ring_radii.get(), Gr, num_radii);
  double ring_r_max = ring_radii[num_elements-1]; // Radius of the outer most (last) ring.

  for (dream_idx_type n=0; n<num_elements; n++) {

    double foc_delay = 0.0;
    if (foc_met != FocusMet::ud) {
      foc_delay = focusing_annular(foc_met, focal[0], ring_radii[n], ring_r_max, cp);
    } else {
      foc_delay = focusing_annular(foc_met, focal[n], ring_radii[n], ring_r_max, cp);
    }

    double weight = 1.0;
    if (do_apod){
      weight = apodization_annular(apod_met, n, apod, ring_radii[n], ring_r_max, apod_par);
    }

    resp_annular(h_disc.get(), h_ring.get(), nt, n); // Compute the impulse response for the n:th ring.
    superpos_annular(h_ring.get(), h, nt, weight, foc_delay, n, dt); // Add the n:th impulse response.
  }

  return out_err;
}

/***
 *
 * Compute the radius of each ring. We use the (outer radius + inner radius)/2.
 *
 ***/

void ArrAnnu::annular_disc_radii(double *ring_r, double *Gr, dream_idx_type num_radii)
{
  //ring_r[0] = Gr[0]/2.0; // FIXME Should we use 0.0?
  ring_r[0] = 0.0;
  for (dream_idx_type n=2, ns=1; n<num_radii; n += 2, ns++) {
    ring_r[ns] = (Gr[n]+Gr[n-1]) / 2.0; // "Middle" radius.
  }
}

/***
*
* Focus delay for the annular array
*
***/

double ArrAnnu::focusing_annular(FocusMet foc_met, double focal, double ring_r, double ring_r_max, double cp)
{
  double foc_delay=0.0;

  switch  (foc_met) {

  case FocusMet::x:
  case FocusMet::y:
  case FocusMet::xy:
  case FocusMet::x_y:
    {
      double rmax = std::sqrt(ring_r_max*ring_r_max + focal*focal);
      double diff = rmax - std::sqrt(ring_r*ring_r + focal*focal);
      foc_delay = diff * 1000 / cp;
    }
    break;

  case FocusMet::ud:
    foc_delay = focal; // Here focal is the user defined time delay in [us] (not the focal depth)
    break;

  case FocusMet::none:
  default:
    foc_delay = 0.0;
    break;

  }

  return foc_delay;
}

/***
 *
 * Apodization
 *
 ***/

double ArrAnnu::apodization_annular(ApodMet apod_met, dream_idx_type n, double *apod, double ring_r,
                           double ring_r_max, double apod_par)
{
  double weight=1.0;

  switch(apod_met) {

  case ApodMet::ud:
    weight = apod[n];
    break;

  case ApodMet::triangle:
    weight = 1.0 - fabs(ring_r) / ring_r_max;
    break;

  case ApodMet::gauss:
    weight = exp(-(apod_par * ring_r*ring_r) / (ring_r_max*ring_r_max));
    break;

  case ApodMet::raised_cosine:
    weight = apod_par + std::cos(ring_r * M_PI / ring_r_max);
    break;

  case ApodMet::simply_supported:
    weight = 1.0 - ring_r*ring_r / (ring_r_max*ring_r_max);
    break;

  case ApodMet::clamped:
    weight = (1.0 - ring_r*ring_r / (ring_r_max*ring_r_max)) * (1.0 - ring_r*ring_r / (ring_r_max*ring_r_max));
    break;

  default:
    break;
  }

  return weight;
}

/***
 *
 * superpos_annular - Superposition the n:th element
 *
 * h : output response
 * h_ring : input response of actual element
 *
 ***/

void ArrAnnu::superpos_annular(double *h_ring, double *h, dream_idx_type nt,
                      double weight, double foc_delay, dream_idx_type n, double dt)
{
  double *buf = (double*) malloc(2*nt*sizeof(double));
  for (dream_idx_type i=0; i<2*nt; i++) {
    buf[i] = 0.0;
  }

  dream_idx_type delay_idx = (dream_idx_type) std::rint(foc_delay/dt) + 1;
  for (dream_idx_type i=0; i<nt; i++) { // FIXME: can delay_idx be > nt?
    if (delay_idx < nt) {
      buf[i + delay_idx] = h_ring[i + n*nt];
    }
  }

  for (dream_idx_type i=0; i<nt; i++) {
    h[i] += weight * buf[i];
  }

  free(buf);

  return;
}

/***
 *
 * Computes the response from one ring by subtracting the inner disc response from
 * the outer disc response.
 *
 ***/

void ArrAnnu::resp_annular(double *h_disc, double *h_ring, dream_idx_type nt, dream_idx_type n)
{
  dream_idx_type k = 2*n;

  if (k == 0) {
    for (dream_idx_type i=0; i<nt; i++) { // Center element.
      h_ring[i + n*nt] = h_disc[i];
    }
  } else {
    for (dream_idx_type i=0; i <nt; i++) { // n:th ring
      h_ring[i + n*nt] = h_disc[i + k*nt] - h_disc[i + (k-1)*nt]; // Outer response - inner response.
    }
  }

  return;
}

/***
 *
 * dream_arr_annu - Threaded annular array function
 *
 ***/

ErrorLevel ArrAnnu::dream_arr_annu(double alpha,
                                   double *ro, dream_idx_type no,
                                   double dx, double dy, double dt, dream_idx_type nt,
                                   DelayType delay_type, double *delay,
                                   double v, double cp,
                                   dream_idx_type num_radii, double *Gr,
                                   FocusMet foc_met, double *focal,
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

  // Check if we have attenuation
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
    D[thread_n].dx = dx;
    D[thread_n].dy = dy;
    D[thread_n].dt = dt;
    D[thread_n].nt = nt;
    D[thread_n].delay_type = delay_type;
    D[thread_n].delay = delay;
    D[thread_n].v = v;
    D[thread_n].cp = cp;
    D[thread_n].att = att_ptr;
    D[thread_n].num_radii = num_radii;
    D[thread_n].Gr = Gr;
    D[thread_n].foc_met = foc_met;
    D[thread_n].do_apod = do_apod;
    D[thread_n].apod_met = apod_met;
    D[thread_n].focal = focal;
    D[thread_n].apod = apod;
    D[thread_n].apod_par = apod_par;
    D[thread_n].h = h;
    D[thread_n].err_level = err_level;

    if (nthreads > 1) {
      // Start the threads.
      threads[thread_n] = arr_annu_thread(&D[thread_n]);
      set_dream_thread_affinity(thread_n, nthreads, threads);
    } else {
      smp_dream_arr_annu(&D[0]);
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
