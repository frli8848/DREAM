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

#include "das_uni.h"
#include "arg_parser.h"

#include "mex.h"

// Persistent smart pointer to DAS object.
// This one is only cleared when we do a
// >> clear das
std::unique_ptr<DAS_UNI<double>> das_uni_d=nullptr;
std::unique_ptr<DAS_UNI<float>> das_uni_f=nullptr;

/***
 *
 * Matlab (MEX) gateway function for das (delay-and-sum).
 *
 ***/

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  ArgParser ap;

  // Check for proper number of arguments

  if (nrhs != 8) {
    dream_err_msg("das requires 8 input arguments!");
  }

  if (nlhs > 1) {
    dream_err_msg("Too many output arguments for das!");
  }

  //
  // Check if we are using single or double precision processsing
  //

  // Args 1 to 7 must have the same datatype (float or double).
  bool use_float = mxIsSingle(prhs[0]);
  if (use_float) {
    if (mxIsSingle(prhs[1]) != use_float ||
        mxIsSingle(prhs[2]) != use_float ||
        mxIsSingle(prhs[3]) != use_float ||
        mxIsSingle(prhs[4]) != use_float ||
        mxIsSingle(prhs[5]) != use_float ||
        mxIsSingle(prhs[6]) != use_float) {
      dream_err_msg("First arg is single precision but one of args 2 to 7 is not!");
      return;
    }
  } else {
    if (!(mxIsDouble(prhs[1]) &&
          mxIsDouble(prhs[2]) &&
          mxIsDouble(prhs[3]) &&
          mxIsDouble(prhs[4]) &&
          mxIsDouble(prhs[5]) &&
          mxIsDouble(prhs[6])) ) {
      dream_err_msg("First arg is double precision but one of args 2 to 7 is not!");
      return;
    }
  }

  //
  // Data
  //

  dream_idx_type a_scan_len = mxGetM(prhs[0]); // A-scan length
  dream_idx_type num_a_scans = mxGetN(prhs[0]);

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

  if (mxGetM(prhs[1])*mxGetN(prhs[1]) != 3) {
    dream_err_msg("Argument 2 must be a 3 element vector!");
  }

  float min_t_f = 0.0;
  float pitch_t_f = 0.0;
  float max_t_f = 0.0;

  double min_t_d = 0.0;
  double pitch_t_d = 0.0;
  double max_t_d = 0.0;

  if (use_float) {
    min_t_f   = ((float*) mxGetData(prhs[1]))[0];
    pitch_t_f = ((float*) mxGetData(prhs[1]))[1];
    max_t_f   = ((float*) mxGetData(prhs[1]))[2];
  } else {
    min_t_d   = mxGetPr(prhs[1])[0];
    pitch_t_d = mxGetPr(prhs[1])[1];
    max_t_d   = mxGetPr(prhs[1])[2];
  }

  //
  // Recieve array
  //

  float min_r_f = 0.0;
  float pitch_r_f = 0.0;
  float max_r_f = 0.0;

  double min_r_d = 0.0;
  double pitch_r_d = 0.0;
  double max_r_d = 0.0;

  if (mxGetN(prhs[2]) != 0) {   // SAFT is using an empty Gr.

    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 3) {
      dream_err_msg("Argument 3 must be a 3 element vector!");
    }

    if (use_float) {
      min_r_f   = ((float*) mxGetData(prhs[2]))[0];
      pitch_r_f = ((float*) mxGetData(prhs[2]))[1];
      max_r_f   = ((float*) mxGetData(prhs[2]))[2];
    } else {
      min_r_d   = mxGetPr(prhs[2])[0];
      pitch_r_d = mxGetPr(prhs[2])[1];
      max_r_d   = mxGetPr(prhs[2])[2];
    }
  }

  //
  // Observation point (in a uniform grid).
  //

  // Check that arg is a 3 x 3 matrix
  if ((mxGetN(prhs[3]) != 3) || (mxGetN(prhs[3])!=3)) {
    dream_err_msg("Argument 4 must be a 3 x 3 matrix!");
  }

  dream_idx_type Nx = 0, Ny = 0, Nz = 0;

  float min_Rx_f=0.0, dx_f=0.0, max_Rx_f=0.0;
  float min_Ry_f=0.0, dy_f=0.0, max_Ry_f=0.0;
  float min_Rz_f=0.0, dz_f=0.0, max_Rz_f=0.0;

  double min_Rx_d=0.0, dx_d=0.0, max_Rx_d=0.0;
  double min_Ry_d=0.0, dy_d=0.0, max_Ry_d=0.0;
  double min_Rz_d=0.0, dz_d=0.0, max_Rz_d=0.0;

  if (use_float) {
    min_Rx_f = ((float*) mxGetData(prhs[3]))[0]; dx_f = ((float*) mxGetData(prhs[3]))[3]; max_Rx_f = ((float*) mxGetData(prhs[3]))[6];
    min_Ry_f = ((float*) mxGetData(prhs[3]))[1]; dy_f = ((float*) mxGetData(prhs[3]))[4]; max_Ry_f = ((float*) mxGetData(prhs[3]))[7];
min_Rz_f = ((float*) mxGetData(prhs[3]))[2]; dz_f = ((float*) mxGetData(prhs[3]))[5]; max_Rz_f = ((float*) mxGetData(prhs[3]))[8];
    Nx = (dream_idx_type) ((max_Rx_f - min_Rx_f)/dx_f+1.0);
    Ny = (dream_idx_type) ((max_Ry_f - min_Ry_f)/dy_f+1.0);
    Nz = (dream_idx_type) ((max_Rz_f - min_Rz_f)/dz_f+1.0);
  } else {
    min_Rx_d = mxGetPr(prhs[3])[0]; dx_d = mxGetPr(prhs[3])[3]; max_Rx_d = mxGetPr(prhs[3])[6];
    min_Ry_d = mxGetPr(prhs[3])[1]; dy_d = mxGetPr(prhs[3])[4]; max_Ry_d = mxGetPr(prhs[3])[7];
    min_Rz_d = mxGetPr(prhs[3])[2]; dz_d = mxGetPr(prhs[3])[5]; max_Rz_d = mxGetPr(prhs[3])[8];
    Nx = (dream_idx_type) ((max_Rx_d - min_Rx_d)/dx_d+1.0);
    Ny = (dream_idx_type) ((max_Ry_d - min_Ry_d)/dy_d+1.0);
    Nz = (dream_idx_type) ((max_Rz_d - min_Rz_d)/dz_d+1.0);
  }

  dream_idx_type No = Nx*Ny*Nz;

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
  ap.parse_das_arg("das_uni", prhs, 7, das_type);

//
  // Init DAS and output arg.
  //

  float *Im_f = nullptr;
  double *Im_d = nullptr;

  bool init_das = true;

  if (use_float) {

    // Create an output matrix for the impulse response.
    plhs[0] = mxCreateNumericMatrix(No, 1, mxSINGLE_CLASS , mxREAL);
    Im_f = (float*) mxGetData(plhs[0]);
    //
    // Single precision DAS_UNI
    //

    if (das_uni_f) { // das_uni (float) object exist - check if we can reuse previous das_uni init
      if (!das_uni_f->das_setup_has_changed(das_type, a_scan_len,
                                            min_t_f, pitch_t_f, max_t_f,
                                            min_r_f, pitch_r_f, max_r_f,
                                            min_Rx_f, dx_f, max_Rx_f,
                                            min_Ry_f, dy_f, max_Ry_f,
                                            min_Rz_f, dz_f, max_Rz_f)) {
        init_das = false;
      }
    }

    if (init_das) {

      if (das_uni_d) {
        das_uni_d = nullptr; // Release the double object if it exist to free, in particular, GPU resources (call destructor).
      }

      try {
        das_uni_f = std::make_unique<DAS_UNI<float>>(das_type, a_scan_len,
                                                     min_t_f, pitch_t_f, max_t_f,
                                                     min_r_f, pitch_r_f, max_r_f,
                                                     min_Rx_f, dx_f, max_Rx_f,
                                                     min_Ry_f, dy_f, max_Ry_f,
                                                     min_Rz_f, dz_f, max_Rz_f);
      }

      catch (std::runtime_error &err) {
        std::cout << err.what();
        return;
      }
    }

  } else {

    //
    // Double precision DAS_UNI
    //

    // Create an output matrix for the impulse response.
    plhs[0] = mxCreateDoubleMatrix(No,1,mxREAL);
    Im_d = mxGetPr(plhs[0]);


    if (das_uni_d) { // das_uni (double) object exist - check if we can reuse previous das_uni init
      if (!das_uni_d->das_setup_has_changed(das_type, a_scan_len,
                                            min_t_d, pitch_t_d, max_t_d,
                                            min_r_d, pitch_r_d, max_r_d,
                                            min_Rx_d, dx_d, max_Rx_d,
                                            min_Ry_d, dy_d, max_Ry_d,
                                            min_Rz_d, dz_d, max_Rz_d)) {
        init_das = false;
      }
    }

    if (init_das) {

      if (das_uni_f) {
        das_uni_f = nullptr; // Release the float object if it exist to free, in particular, GPU resources (call destructor).
      }

      try {
        das_uni_d = std::make_unique<DAS_UNI<double>>(das_type, a_scan_len,
                                            min_t_d, pitch_t_d, max_t_d,
                                            min_r_d, pitch_r_d, max_r_d,
                                            min_Rx_d, dx_d, max_Rx_d,
                                            min_Ry_d, dy_d, max_Ry_d,
                                            min_Rz_d, dz_d, max_Rz_d);
      }

      catch (std::runtime_error &err) {
        std::cout << err.what();
        return;
      }
    }
  }

  if (use_float) { // Single precision
    das_uni_f->cl_das_uni(Yf,
                          min_t_f, pitch_t_f, max_t_f,
                          min_r_f, pitch_r_f, max_r_f,
                          min_Rx_f, dx_f, max_Rx_f,
                          min_Ry_f, dy_f, max_Ry_f,
                          min_Rz_f, dz_f, max_Rz_f,
                          dt_f,
                          delay_f[0], cp_f,
                          Im_f);
  } else { // Double precision
    das_uni_d->cl_das_uni(Yd,
                          min_t_d, pitch_t_d, max_t_d,
                          min_r_d, pitch_r_d, max_r_d,
                          min_Rx_d, dx_d, max_Rx_d,
                          min_Ry_d, dy_d, max_Ry_d,
                          min_Rz_d, dz_d, max_Rz_d,
                          dt_d,
                          delay_d[0], cp_d,
                          Im_d);
  }

  return;
}
