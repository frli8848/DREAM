/***
*
* Copyright (C) 2003,2006,2007,2008,2009,2014,2015,2016,2021,2023,2024 Fredrik Lingvall
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
std::unique_ptr<DAS<double>> das_d=nullptr;
std::unique_ptr<DAS<float>> das_f=nullptr;

/***
 *
 * Matlab (MEX) gateway function for das (delay-and-sum).
 *
 ***/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  ArgParser ap;

  // Check for proper number of arguments

  if ((nrhs < 8) || (nrhs > 10)) {
    dream_err_msg("das requires 8 to 10 input arguments!");
  }

  if (nlhs > 2) {
    dream_err_msg("Too many output arguments for das!");
  }

  //
  // Check if we are using single or double precision processsing
  //

  // Args 1, 2, 3, 4, and 6 must have the same datatype (float or double).
  bool use_float = mxIsSingle(prhs[0]);
  if (use_float) {
    if (mxIsSingle(prhs[1]) != use_float ||
        mxIsSingle(prhs[2]) != use_float ||
        mxIsSingle(prhs[3]) != use_float ||
        mxIsSingle(prhs[5]) != use_float) {
      dream_err_msg("First arg is single precision but one of arg 2, 3, 4, or 6 is not!");
      return;
    }
  } else {
    if (!(mxIsDouble(prhs[1]) &&
          mxIsDouble(prhs[2]) &&
          mxIsDouble(prhs[3]) &&
          mxIsDouble(prhs[5])) ) {
      dream_err_msg("First arg is double precision but one of arg 2, 3, 4, or 6 is not!");
      return;
    }
  }


  //
  // Data
  //

  dream_idx_type a_scan_len = mxGetM(prhs[0]); // A-scan length

  float *Yf = nullptr;
  double *Yd = nullptr;

  if (use_float) {
    Yf = (float*) mxGetData(prhs[0]);
  } else {
    Yd = mxGetPr(prhs[0]);
  }

  //
  // Transmit array
  //

  if (mxGetN(prhs[1]) != 3) {
    dream_err_msg("Argument 2 must be a (num transmit elements) x 3 matrix!");
  }

  dream_idx_type num_t_elements = mxGetM(prhs[1]);

  const float *Gt_f=nullptr;
  const double *Gt_d=nullptr;

  if (use_float) {
    Gt_f = (float*) mxGetData(prhs[1]);
  } else {
    Gt_d = mxGetPr(prhs[1]);
  }

  //
  // Recieve array
  //

  dream_idx_type num_r_elements = 0;

  const float *Gr_f = nullptr;
  const double *Gr_d = nullptr;


  if (mxGetN(prhs[2]) != 0) {   // SAFT is using an empty Gr.

    if (mxGetN(prhs[2]) != 3) {
      dream_err_msg("Argument 3 must be a (num recieve elements) x 3 matrix!");
    }

    num_r_elements = mxGetM(prhs[2]);

    if (use_float) {
      Gr_f = (float*) mxGetData(prhs[2]);
    } else {
      Gr_d = mxGetPr(prhs[2]);
    }
  }

  //
  // Observation point.
  //

  // Check that arg (number of observation points) x 3 matrix
  if (!(mxGetN(prhs[3])==3)) {
    dream_err_msg("Argument 4 must be a (number of observation points) x 3 matrix!");
  }

  dream_idx_type No = mxGetM(prhs[3]); // Number of observation points.

  const float *Ro_f = nullptr;
  const double *Ro_d = nullptr;

  if (use_float) {
    Ro_f = (float*) mxGetData(prhs[3]);
  } else {
    Ro_d = mxGetPr(prhs[3]);
  }

  //
  // Temporal sampling parameter.
  //

  if (!(mxGetM(prhs[4]) == 1 && mxGetN(prhs[4]) == 1)) {
    dream_err_msg("Argument 5 must be a scalar!");
  }

  float dt_f = 0.0;
  double dt_d = 0.0;

  if (use_float) {
    float *s_par_f = (float*) mxGetData(prhs[4]);
    dt_f = s_par_f[0];
  } else {
    double *s_par_d = mxGetPr(prhs[4]);
    dt_d = s_par_d[0];
  }

  //
  // Start point of impulse response vector ([us]).
  //

  ap.check_delay("das", prhs, 5, No);

  DelayType delay_type = DelayType::single;  // delay is a scalar.
  if (mxGetM(prhs[5]) * mxGetN(prhs[5]) != 1) {
    delay_type = DelayType::multiple; // delay is a vector.
  }

  const float *delay_f = nullptr;
  const double *delay_d = nullptr;

  if (use_float) {
    delay_f = (float*) mxGetData(prhs[5]);
  } else {
    delay_d = mxGetPr(prhs[5]);
  }

  //
  // Material parameter
  //

  // Check that arg 7 is a scalar.
  if (!(mxGetM(prhs[6]) == 1 && mxGetN(prhs[6]) == 1)) {
    dream_err_msg("Argument 7 must be a scalar!");
  }

  float cp_f = 0.0;
  double cp_d = 0.0;

  if (use_float) {
    float *m_par_f = (float*) mxGetData(prhs[6]);
    cp_f = m_par_f[0];
  } else {
    double *s_par_d = mxGetPr(prhs[6]);
    cp_d = s_par_d[0];
  }

  //
  // DAS method
  //

  DASType das_type;
  ap.parse_das_arg("das", prhs, 7, das_type);

  //
  // Error reporting.
  //

  ErrorLevel err_level=ErrorLevel::stop;

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

  //
  // Init DAS and output arg.
  //

  float *Im_f = nullptr;
  double *Im_d = nullptr;

  bool init_das = true;

  // We use this to make the GPU init code run
  bool use_gpu = false;
#ifdef USE_OPENCL
  if (device == "gpu") {
    use_gpu = true;
  }
#endif

  if (use_float) {

    // Create an output matrix for the impulse response.
    plhs[0] = mxCreateNumericMatrix(No, 1, mxSINGLE_CLASS , mxREAL);
    Im_f = (float*) mxGetData(plhs[0]);

    //
    // Single precision DAS
    //

    if (das_f) { // das (float) object exist - check if we can reuse previous das init
      if (!das_f->das_setup_has_changed(das_type, a_scan_len, No, num_t_elements, num_r_elements, use_gpu)) {
        init_das = false;
      }
    }

    if (init_das) {

      if (das_d) {
        das_d = nullptr; // Release the double obejct if it exist to free, in particular, GPU resources (call destructor).
      }

      try {
        das_f = std::make_unique<DAS<float>>(das_type, a_scan_len, No, num_t_elements, num_r_elements, use_gpu);
      }

      catch (std::runtime_error &err) {
        std::cout << err.what();
        return;
      }
    }

    das_d->set_running();

    // Register signal handler.
    std::signal(SIGABRT, DAS<float>::abort);

  } else {

    //
    // Double precision DAS
    //

    // Create an output matrix for the impulse response.
    plhs[0] = mxCreateDoubleMatrix(No,1,mxREAL);
    Im_d = mxGetPr(plhs[0]);

    if (das_d) { // das (double) object exist - check if we can reuse previous das init
      if (!das_d->das_setup_has_changed(das_type, a_scan_len, No, num_t_elements, num_r_elements, use_gpu)) {
        init_das = false;
      }
    }

    if (init_das) {

      if (das_f) {
        das_f = nullptr; // Release the float object if it exist to free, in particular, GPU resources (call destructor).
      }

      try {
        das_d = std::make_unique<DAS<double>>(das_type, a_scan_len, No, num_t_elements, num_r_elements, use_gpu);
      }

      catch (std::runtime_error &err) {
        std::cout << err.what();
        return;
      }
    }

    das_d->set_running();

    // Register signal handler.
    std::signal(SIGABRT, DAS<double>::abort);
  }

  SIRError err = SIRError::none;

#ifdef USE_OPENCL

  // Check if we should use the GPU
  if (device == "gpu") {

    if (use_float) { // Single precision
      das_f->cl_das(Yf, Ro_f, Gt_f, Gr_f, dt_f, delay_f[0], cp_f, Im_f);
    } else { // Double precision
      das_d->cl_das(Yd, Ro_d, Gt_d, Gr_d, dt_d, delay_d[0], cp_d, Im_d);
    }

  } else { // Otherwise use the cpu

#endif

    if (device == "gpu") {
      std::cout << "Warning: Compute device set to 'gpu' but DREAM is build without OpenCL support!" << std::endl;
      std::cout << "Using the CPU backend!" << std::endl;
    }

    if (use_float) { // Single precision

      err = das_f->das(Yf, Ro_f, Gt_f, Gr_f, dt_f, delay_type, delay_f, cp_f, Im_f, err_level);
      if (!das_f->is_running()) {
        if (err != SIRError::out_of_bounds) {
          dream_err_msg("CTRL-C pressed!\n"); // Bail out.
        } else {
          dream_err_msg("SIR out-of-bounds!\n"); // Bail out.
        }

        return;
      }

    } else { // Double precision

      err = das_d->das(Yd, Ro_d, Gt_d, Gr_d, dt_d, delay_type, delay_d, cp_d, Im_d, err_level);
      if (!das_d->is_running()) {
        if (err != SIRError::out_of_bounds) {
          dream_err_msg("CTRL-C pressed!\n"); // Bail out.
        } else {
          dream_err_msg("SIR out-of-bounds!\n"); // Bail out.
        }

        return;
      }

    }

#ifdef USE_OPENCL
  }
#endif

  if (err == SIRError::out_of_bounds) {
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
