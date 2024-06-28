/***
 *
 * Copyright (C) 2003,2006,2007,2008,2009,2014,2015,2019,2021,2023,2024 Fredrik Lingvall
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

#include <csignal>

#include "dream_arr_cylind.h"
#include "arg_parser.h"

#include "mex.h"

/***
 *
 * Gateway function for (parallel) dream_arr_cylind.
 *
 ***/

extern void _main();

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  ArgParser ap;

  // Check for proper number of arguments

  ap.check_arg_in("dream_arr_cylind", nrhs, 13, 14);
  ap.check_arg_out("dream_arr_cylind", nlhs, 0, 2);

  //
  // Observation point.
  //

  ap.check_obs_points("dream_arr_cylind", prhs, 0);
  dream_idx_type No = mxGetM(prhs[0]); // Number of observation points.
  double *Ro = mxGetPr(prhs[0]);

  //
  // Transducer geometry
  //

  double a=0.0, b=0.0, Rcurv=0.0;
  ap.parse_geometry("dream_arr_cylind", prhs, 1, 3, a, b, Rcurv);

  //
  // Grid function (position vectors of the elements).
  //

  ap.check_array("dream_arr_cylind", prhs, 2);

  dream_idx_type num_elements = mxGetM(prhs[2]); // Number of elementents in the array.
  double *G = mxGetPr(prhs[2]);         // First column in the matrix.
  //gy    = gx + num_elements;          // Second column in the matrix.
  //gz    = gy + num_elements;          // Third column in the matrix.

  //
  // Temporal and spatial sampling parameters.
  //

  double dx=0.0, dy=0.0, dt=0.0;
  dream_idx_type nt=0;
  ap.parse_sampling("dream_arr_cylind", prhs, 3, 4, dx, dy, dt, nt);

  //
  // Start point of impulse response vector ([us]).
  //

  ap.check_delay("dream_arr_cylind", prhs, 4, No);
  double *delay = mxGetPr(prhs[4]);

  DelayType delay_type = DelayType::single;  // delay is a scalar.
  if (mxGetM(prhs[4]) * mxGetN(prhs[4]) != 1) {
    delay_type = DelayType::multiple; // delay is a vector.
  }

  //
  // Material parameters
  //

  double v=1.0, cp=1000.0, alpha=0.0;
  ap.parse_material("dream_arr_cylind", prhs, 5, v, cp, alpha);

  //
  // Focusing parameters.
  //

  // Allocate memory for the user defined focusing delays
  std::unique_ptr<double[]> focal = std::make_unique<double[]>(num_elements);

  FocusMet foc_met=FocusMet::none;

  if (nrhs >= 7) {
    ap.parse_focus_args("dream_arr_cylind", prhs, 6, foc_met, focal.get(), num_elements);
  } else {
    foc_met = FocusMet::none;
  }

  //
  // Beam steering.
  //

  SteerMet steer_met=SteerMet::none;
  double theta=0.0, phi=0.0;

  if (nrhs >= 9) {
    ap.parse_steer_args("dream_arr_cylind", prhs, 8, steer_met, theta, phi);
  } else {
    steer_met = SteerMet::none;
  }

  //
  // Apodization.
  //

  // Allocate memory for the user defined apodization weights.
  std::unique_ptr<double[]> apod = std::make_unique<double[]>(num_elements);

  double apod_par=0.0;
  bool do_apod=false;
  ApodMet apod_met=ApodMet::gauss;

  if (nrhs >= 11) {
    ap.parse_apod_args("dream_arr_cylind", prhs, 10, num_elements,
                       do_apod, apod.get(), apod_met, apod_par);
  } else {
    do_apod = false;
  }

  //
  // Error reporting.
  //

  ErrorLevel err_level=ErrorLevel::stop;

  if (nrhs == 14) {
    ap.parse_error_arg("dream_arr_cylind", prhs, 13, err_level);
  } else {
    err_level = ErrorLevel::stop; // Default.
  }

  // Create an output matrix for the impulse response
  plhs[0] = mxCreateDoubleMatrix(nt,No,mxREAL);
  double *h = mxGetPr(plhs[0]);

  SIRData hsir(h, nt, No);
  hsir.clear();

  ArrCylind arr_cylind;

  // Register signal handler.
  std::signal(SIGINT, ArrCylind::abort);

  //
  // Call the DREAM subroutine.
  //

  SIRError err = arr_cylind.dream_arr_cylind(alpha,
                                             Ro, No,
                                             a, b, Rcurv,
                                             dx, dy,  dt, nt,
                                             delay_type, delay,
                                             v, cp,
                                             num_elements, G,
                                             foc_met,  focal.get(),
                                             steer_met, theta, phi,
                                             apod.get(), do_apod, apod_met, apod_par,
                                             h, err_level);

  if (!arr_cylind.is_running()) {
    if (err != SIRError::out_of_bounds) {
      dream_err_msg("CTRL-C pressed!\n"); // Bail out.
    } else {
      dream_err_msg("SIR out-of-bounds!\n"); // Bail out.
    }
  }

  //
  // Check for Error.
  //

  if (err == SIRError::out_of_bounds) {
    dream_err_msg(""); // Bail out if error.
  }

  //
  // Return error.
  //

  if (nlhs == 2) {
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    double *err_p =  mxGetPr(plhs[1]);
    err_p[0] = (double) err;
  }

  return;
}
