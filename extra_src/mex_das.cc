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

// Persistent smart pointer to DAS object.
// This one is only cleared when we do a
// >> clear das
std::unique_ptr<DAS> das=nullptr;

/***
 *
 * Matlab (MEX) gateway function for das (delay-and-sum).
 *
 ***/

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  ArgParser ap;

  // Check for proper number of arguments

  if ((nrhs < 8) && (nrhs > 10)) {
    dream_err_msg("das requires 8 to 10 input arguments!");
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

  dream_idx_type No = mxGetM(prhs[3]); // Number of observation points.
  double *Ro = mxGetPr(prhs[3]);

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

  ap.check_delay("das", prhs, 5, No);
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

  if (nrhs >= 9) {
    ap.parse_error_arg("das", prhs, 8, err_level);
  } else {
    err_level = ErrorLevel::stop; // Default.
  }

  //
  // Compute device
  //

  std::string device;

  if (nrhs == 10) {

    if (!mxIsChar(prhs[9])) {
      dream_err_msg("Argument 10 must be a string");
    }

    device =  ap.get_string_arg(prhs, 9);
  }

  // Create an output matrix for the impulse response
  plhs[0] = mxCreateDoubleMatrix(No,1,mxREAL);
  double *Im = mxGetPr(plhs[0]);

  bool init_das = true;
  if (das) { // das object exist - check if we can reuse previous das init
    if (!das->das_setup_has_changed(das_type, a_scan_len, No, num_t_elements, num_r_elements)) {
      init_das = false;
    }
  }

  if (init_das) {
    try {
      das = std::make_unique<DAS>(das_type, a_scan_len, No, num_t_elements, num_r_elements);
    }

    catch (std::runtime_error &err) {
      std::cout << err.what();
      return;
    }
  }

  das->set_running();

  // Register signal handler.
  std::signal(SIGABRT, DAS::abort);

#ifdef USE_OPENCL

  // Check if we should use the GPU
  if (device == "gpu" && num_r_elements > 0) { // SAFT is most likely fast enough on the CPU.

    das->cl_das(Y, Ro, Gt, Gr, dt, delay[0], cp, Im);

  } else { // Otherwise use the cpu

#endif

    if (device == "gpu") {
      std::cout << "Warning: Compute device set to 'gpu' but DREAM is build without OpenCL support!" << std::endl;
      std::cout << "Using the CPU backend!" << std::endl;
    }

    err = das->das(Y, Ro, Gt, Gr, dt, delay_type, delay, cp, Im, err_level);

#ifdef USE_OPENCL
  }
#endif

  if (!das->is_running()) {
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
