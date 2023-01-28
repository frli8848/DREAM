/***
*
* Copyright (C) 2003,2006,2007,2008,2009,2014,2015,2016,2021,2023 Fredrik Lingvall
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

#include<csignal>

#include "das.h"
#include "arg_parser.h"

#include "mex.h"

/***
 *
 * Matlab (MEX) gateway function for das (delay-and-sum).
 *
 ***/

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  ArgParser ap;

  // Check for proper number of arguments

  if (!((nrhs == 8) || (nrhs == 9))) {
    dream_err_msg("das requires 8 or 9 input arguments!");
  }

  if (nlhs > 2) {
    dream_err_msg("Too many output arguments for das!");
  }

  //
  // Data
  //

  dream_idx_type a_scan_len = mxGetM(prhs[0]); // A-scan length
  dream_idx_type num_a_scans = mxGetN(prhs[0]);
  double *Y = mxGetPr(prhs[0]);

  //
  // Transmit array
  //

  if (mxGetN(prhs[1]) != 3) {
    dream_err_msg("Argument 2 must be a (num transmit elements) x 3 matrix!");
  }

  dream_idx_type num_t_elements = mxGetM(prhs[1]);
  double *Gt = mxGetPr(prhs[1]);

  //
  // Recieve array
  //

  dream_idx_type num_r_elements = 0;
  double *Gr = nullptr;

  if (mxGetN(prhs[2]) != 0) { // Check if we do SAFT or TFM

    if (mxGetN(prhs[2]) != 3) {
      dream_err_msg("Argument 3 must be a (num recieve elements) x 3 matrix!");
    }

    num_r_elements = mxGetM(prhs[2]);
    Gr = mxGetPr(prhs[2]);
  }

  //
  // Observation point.
  //

  // Check that arg (number of observation points) x 3 matrix
  if (!(mxGetN(prhs[3])==3)) {
    dream_err_msg("Argument 4 must be a (number of observation points) x 3 matrix!");
  }

  dream_idx_type no = mxGetM(prhs[3]); // Number of observation points.
  double *ro = mxGetPr(prhs[3]);

  //
  // Temporal and spatial sampling parameters.
  //

  if (!(mxGetM(prhs[4]) == 1 && mxGetN(prhs[4]) == 1)) {
    dream_err_msg("Argument 5 must be a scalar!");
  }

  double *s_par = mxGetPr(prhs[4]);
  double dt = s_par[0]; // Temporal discretization size (= 1/sampling freq).

  //
  // Start point of impulse response vector ([us]).
  //

  ap.check_delay("das", prhs, 5, no);
  double *delay = mxGetPr(prhs[5]);

  DelayType delay_type = DelayType::single;  // delay is a scalar.
  if (mxGetM(prhs[5]) * mxGetN(prhs[5]) != 1) {
    delay_type = DelayType::multiple; // delay is a vector.
  }

  //
  // Material parameters
  //

  // Check that arg 7 is a scalar.
  if (!(mxGetM(prhs[6]) == 1 && mxGetN(prhs[6]) == 1)) {
    dream_err_msg("Argument 7 must be a scalar!");
  }

  double *m_par = mxGetPr(prhs[6]);
  double cp = m_par[0]; // Sound speed.

  //
  // DAS method
  //

  DASType das_type;
  ap.parse_das_arg("das", prhs, 7, das_type);

  //
  // Error reporting.
  //

  ErrorLevel err=ErrorLevel::none, err_level=ErrorLevel::stop;

  if (nrhs == 9) {
    ap.parse_error_arg("das", prhs, 8, err_level);
  } else {
    err_level = ErrorLevel::stop; // Default.
  }

  // Create an output matrix for the impulse response
  plhs[0] = mxCreateDoubleMatrix(no,1,mxREAL);
  double *Im = mxGetPr(plhs[0]);

  DAS das;

  // Register signal handler.
  std::signal(SIGABRT, DAS::abort);

  err = das.das(Y, a_scan_len,
                ro,  no,
                Gt, num_t_elements,
                Gr, num_r_elements,
                dt,
                delay_type, delay,
                cp,
                das_type,
                Im,
                err_level);

  if (!das.is_running()) {
    dream_err_msg("CTRL-C pressed!\n"); // Bail out.
  }

  if (err == ErrorLevel::stop) {
    dream_err_msg("Error in DAS"); // Bail out if error.
  }

  // Return error.
  if (nlhs == 2) {
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    double *err_p =  mxGetPr(plhs[1]);
    err_p[0] = (double) err;
  }

  return;
}
