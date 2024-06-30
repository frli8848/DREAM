/***
*
* Copyright (C) 2003,2005,2006,2007,2008,2009,2014,2015,2019,2021,2023,2024 Fredrik Lingvall
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
#include "dream_arr_annu.h"
#include "arg_parser.h"

#include "mex.h"

/***
 *
 * Matlab (MEX) gateway function for dream_arr_annu
 *
 ***/

extern void _main();

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  ArgParser ap;

  // Check for proper number of arguments.

  ap.check_arg_in("dream_arr_annu", nrhs, 10, 11);
  ap.check_arg_out("dream_arr_annu", nlhs, 0, 2);

  //
  // Observation point.
  //

  ap.check_obs_points("dream_arr_annu", prhs, 0);
  dream_idx_type No = mxGetM(prhs[0]); // Number of observation points.
  double *Ro = mxGetPr(prhs[0]);

  //
  // Grid function (radii vectors of the elements).
  //

  ap.check_array_annu("dream_arr_annu", prhs, 1);
  dream_idx_type num_radii = mxGetM(prhs[1])*mxGetN(prhs[1]); // Center radius + inner and outer radii for each element.
  dream_idx_type num_elements = (num_radii+1)/2; // Number of elementents in the array.
  double *Gr = mxGetPr(prhs[1]);	// Vector of annular radi,

  //
  // Temporal and spatial sampling parameters.
  //

  double dx=0.0, dy=0.0, dt=0.0;
  dream_idx_type nt=0;
  ap.parse_sampling("dream_arr_annu", prhs, 2, 4, dx, dy, dt, nt);

  //
  // Start point of impulse response vector ([us]).
  //

  ap.check_delay("dream_arr_annu", prhs, 3, No);
  double *delay = mxGetPr(prhs[3]);

  DelayType delay_type = DelayType::single;  // delay is a scalar.
  if (mxGetM(prhs[3]) * mxGetN(prhs[3]) != 1) {
    delay_type = DelayType::multiple; // delay is a vector.
  }

  //
  // Material parameters
  //

  double v=1.0, cp=1000.0, alpha=0.0;
  ap.parse_material("dream_arr_annu", prhs, 4, v, cp, alpha);

  //
  // Focusing parameters.
  //

  // Allocate space for the user defined focusing delays
  std::unique_ptr<double[]> focal = std::make_unique<double[]>(num_elements);

  FocusMet foc_met=FocusMet::none;
  if (nrhs >= 6) {
    ap.parse_focus_args("dream_arr_annu", prhs, 5, foc_met, focal.get(), num_elements);
  } else {
    foc_met = FocusMet::none;
  }

  //
  // Apodization.
  //

  // Allocate memory for the user defined apodization weights.
  std::unique_ptr<double[]> apod = std::make_unique<double[]>(num_elements);

  double apod_par=0.0;
  bool do_apod=false;
  ApodMet apod_met=ApodMet::gauss;
  if (nrhs >= 8) {
    ap.parse_apod_args("dream_arr_annu", prhs, 7, num_elements,
                       do_apod, apod.get(), apod_met, apod_par);
  } else {
    do_apod = false;
  }

  //
  // Error Reporting.
  //

  ErrorLevel err_level=ErrorLevel::stop;

  if (nrhs == 11) {
    ap.parse_error_arg("dream_arr_annu", prhs, 10, err_level);
  } else {
    err_level = ErrorLevel::stop; // Default.
  }

  // Create an output matrix for the impulse response
  plhs[0] = mxCreateDoubleMatrix(nt,No,mxREAL);
  double *h = mxGetPr(plhs[0]);

  SIRData hsir(h, nt, No);
  hsir.clear();

  ArrAnnu arr_annu;

  //
  // Call the DREAM subroutine.
  //

  SIRError err = arr_annu.dream_arr_annu(alpha,
                                Ro, No,
                                dx, dy,  dt, nt,
                                delay_type, delay,
                                v, cp,
                                num_radii, Gr,
                                foc_met, focal.get(),
                                apod.get(), do_apod, apod_met, apod_par,
                                h, err_level);

  if (!arr_annu.is_running()) {
    if (err == SIRError::out_of_bounds) {
      dream_err_msg("SIR out-of-bounds!\n"); // Bail out.
    } else {
      dream_err_msg("CTRL-C pressed!\n"); // Bail out.
    }
  }

  if (err == SIRError::warn_out_of_bounds) {
    std::cout << "Warning: SIR out-of-bounds!" << std::endl;
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
