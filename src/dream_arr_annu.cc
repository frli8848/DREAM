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

#include <math.h>
#include <stdlib.h>

#include "dream_arr_annu.h"
#include "dreamcirc.h"
#include "dream_error.h"

//
// Function prototypes
//

void annular_disc_radii(double *ring_r, double *gr, dream_idx_type num_radii);
double focusing_annular(FocusMet foc_met, double focal, double ring_r, double ring_r_max, double cp);
double apodization_annular(int apod_type, dream_idx_type n, double *apod, double ring_r,
                           double ring_r_max, double param);
void resp_annular(double *h_disc, double *h_ring, dream_idx_type nt, dream_idx_type n);
void superpos_annular(double *h_ring, double *h, dream_idx_type nt,
                      double weight, double foc_delay, dream_idx_type n, double dt);

/***
*
* dream_arr_annu
*
* Routine for computation of spatial impulse response of an annular array.
*
***/

int dream_arr_annu(double xo, double yo, double zo,
                   double dx, double dy, double dt,
                   dream_idx_type nt,
                   double delay,
                   double v, double cp,
                   dream_idx_type num_radii, double *gr,
                   FocusMet foc_met, double *focal,
                   double *apod, bool do_apod, int apod_type, double param,
                   double *h, int err_level)
{
  int err = NONE, out_err = NONE;

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

    if (n==0 && gr[n] <= std::numeric_limits<double>::epsilon()  ) {
      for (dream_idx_type i=0; i<nt; i++) {
        h_disc[i] = 0.0; // h(i,1) = 0.0;
      }
    } else {
      err = dreamcirc(xo, yo, zo,
                      gr[n],
                      dx, dy, dt, nt,
                      delay,
                      v, cp, &h_disc[n*nt], err_level);
      if (err != NONE) {
        out_err = err;
      }
    }
  }

  // Compute the "middle" radius for each array ring
  // of the transducer (which is the radius used when
  // focusing).
  annular_disc_radii(ring_radii.get(), gr, num_radii);
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
      weight = apodization_annular(apod_type, n, apod, ring_radii[n], ring_r_max, param);
    }

    resp_annular(h_disc.get(), h_ring.get(), nt, n); // Compute the impulse response for the n:th ring.
    superpos_annular(h_ring.get(), h, nt, weight, foc_delay, n, dt); // Add the n:th impulse response.
  }

  return out_err;
}

int dream_arr_annu(Attenuation &att, FFTCVec &xc_vec, FFTVec &x_vec,
                   double xo, double yo, double zo,
                   double dx, double dy, double dt,
                   dream_idx_type nt,
                   double delay,
                   double v, double cp,
                   dream_idx_type num_radii, double *gr,
                   FocusMet foc_met, double *focal,
                   double *apod, bool do_apod, int apod_type, double param,
                   double *h, int err_level)
{
  int err = NONE, out_err = NONE;

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

    if (n==0 && gr[n] <= std::numeric_limits<double>::epsilon()  ) {
      for (dream_idx_type i=0; i<nt; i++) {
        h_disc[i] = 0.0; // h(i,1) = 0.0;
      }
    } else {
      err = dreamcirc(att, xc_vec, x_vec,
                      xo, yo, zo,
                      gr[n],
                      dx, dy, dt, nt,
                      delay,
                      v, cp, &h_disc[n*nt], err_level);
      if (err != NONE) {
        out_err = err;
      }
    }
  }

  // Compute the "middle" radius for each array ring
  // of the transducer (which is the radius used when
  // focusing).
  annular_disc_radii(ring_radii.get(), gr, num_radii);
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
      weight = apodization_annular(apod_type, n, apod, ring_radii[n], ring_r_max, param);
    }

    resp_annular(h_disc.get(), h_ring.get(), nt, n); // Compute the impulse response for the n:th ring.
    superpos_annular(h_ring.get(), h, nt, weight, foc_delay, n, dt); // Add the n:th impulse response.
  }

  return out_err;
}

/***
 *
 * Compute the radius of each ring.  We use the (outer radius + inner radius)/2.
 *
 ***/

void annular_disc_radii(double *ring_r, double *gr, dream_idx_type num_radii)
{
  //ring_r[0] = gr[0]/2.0; // FIXME Should we use 0.0?
  ring_r[0] = 0.0;
  for (dream_idx_type n=2, ns=1; n<num_radii; n += 2, ns++) {
    ring_r[ns] = (gr[n]+gr[n-1]) / 2.0; // "Middle" radius.
  }
}

/***
*
* Focus delay for the annular array
*
***/

double focusing_annular(FocusMet foc_met, double focal, double ring_r, double ring_r_max, double cp)
{
  double foc_delay=0.0;

  switch  (foc_met) {

  case FocusMet::x:
  case FocusMet::y:
  case FocusMet::xy:
  case FocusMet::x_y:
    {
      double rmax = sqrt(ring_r_max*ring_r_max + focal*focal);
      double diff = rmax - sqrt(ring_r*ring_r + focal*focal);
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
 *  apodization
 *
 ***/

double apodization_annular(int apod_type, dream_idx_type n, double *apod, double ring_r,
                           double ring_r_max, double param)
{
  double weight=1.0;

  switch(apod_type) {

  case APOD_UD:
    weight = apod[n];
    break;

  case APOD_TRIANGLE:
    weight = 1.0 - fabs(ring_r) / ring_r_max;
    break;

  case APOD_GAUSS:
    weight = exp(-(param * ring_r*ring_r) / (ring_r_max*ring_r_max));
    break;

  case APOD_RISED_COSINE:
    weight = param + cos(ring_r * M_PI / ring_r_max);
    break;

  case APOD_SIMPLY_SUPPORTED:
    weight = 1.0 - ring_r*ring_r / (ring_r_max*ring_r_max);
    break;

  case APOD_CLAMPED:
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
 ***/

void superpos_annular(double *h_ring, double *h, dream_idx_type nt,
                      double weight, double foc_delay, dream_idx_type n, double dt)
{
  double *buf = (double*) malloc(2*nt*sizeof(double));
  for (dream_idx_type i=0; i<2*nt; i++) {
    buf[i] = 0.0;
  }

  dream_idx_type delay_idx = (dream_idx_type) rint(foc_delay/dt) + 1;
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

void resp_annular(double *h_disc, double *h_ring, dream_idx_type nt, dream_idx_type n)
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
