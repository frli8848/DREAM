/***
*
* Copyright (C) 2008,2009,2016,2023 Fredrik Lingvall
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

#include <octave/oct.h>

// Persistent smart pointer to DAS object.
// This one is only cleared when we do a
// octave:1> clear das
std::unique_ptr<DAS<double>> das_d=nullptr;
std::unique_ptr<DAS<float>> das_f=nullptr;

/***
 *
 * das - Octave (oct) gateway function for DAS (delay-and-sum).
 *
 ***/

DEFUN_DLD (das, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {}  [Im] = das(Y,Gt,Gr,Ro,dt,delay,cp,method,err_level,device).\n \
\n\
DAS Computes the delay-and-sum processed reconstruction (beamformed image) for\n\
three different array geometries: SAFT, TFM, and RCA. In the synthtic aperture\n\
focusing techinque (SAFT) one use one transmitter and one reciever (the same)\n\
and moves that along the array aperture. In the total focusing technique (TFM) one\n\
first transmit with the first element and then recieve with all elements, transmit with the \n\
second element and again receive with all elements and so on until the last transmit element (ie.,\n\
transmit with all - receive with all). The last method, row-column adressed (RCA) array, is \n\
a variant of TFM where the elements are arranged in a crossed layout to form a 2D array.\n\
\n\
Data matrix:\n\
\n\
@table @code\n\
@item Y\n\
A K x N matrix where K is the A-scan length and the number of A-scans, N, depends\n\
on DAS algorithm selected (see below).\n\
@end table\n\
\n\
Transmit element grid matrix:\n\
\n\
@table @code\n\
\n\
@item Gt\n\
An Lt x 3 matrix, Gt = [x1 y1 z2; x2 y2 z2; ... xLt yoLt zLt]; where Lt is the number of tranmsit elements.\n\
@end table\n\
\n\
Receive element grid matrix:\n\
\n\
@table @code\n\
@item Gr\n\
An Lr x 3 matrix, Gr = Gt = [x1 y1 z2; x2 y2 z2; ... xLr yoLr zLr]; where Lr is the number of tranmsit elements.\n\
@end table\n\
\n\
Observation point(s) ([mm]):\n\
\n\
@table @code\n\
@item Ro\n\
An No x 3 matrix, Ro = [xo1 yo1 zo2; xo2 yo2 zo2; ... xoNo yoNo zoNo]; where No is the number of observation points.\n\
@end table\n\
\n\
Sampling parameter dt: \n\
\n\
@table @code\n\
@item dt\n\
Temporal discretization period (= 1/sampling freq) [us].\n\
@end table\n\
\n\
Data start and pulse delay compensation:\n\
\n\
@table @code\n\
@item  delay\n\
Scalar delay for all observation points or a vector with individual delays for each observation point [us].\n\
@end table\n\
\n\
Sound speed:\n\
\n\
@table @code\n\
@item cp\n\
Sound velocity of the medium [m/s].\n\
\n\
@end table\n\
DAS algorithm:\n\
das_met is a text string parameter for selecting DAS algorithm, options are:\n\
\n\
@table @code\n\
@item 'saft'\n\
When SAFT is selected data Y must be an K x Lt (SAFT is also selected if Gr=[]).\n\
@item 'tfm'\n\
When TFM is selected a linear array is assumed and data Y must be a \n\
K x Lt*Lr matrix.\n\
@item 'rca'\n\
When RCA is selected a 2D RCA array is assumed and data Y must be a \n\
K x Lt*Lr matrix.\n\
@end table\n\
\n\
Error Handling:\n\
err_level is an optional text string parameter for controlling the error behavior, options are:\n\
\n\
@table @code\n\
@item 'ignore'\n\
An error is ignored (no error message is printed and the program is not stopped) but the err output \n\
argument is negative if an error occured.\n\
@item 'warn'\n\
An error message is printed but the program in not stopped (and err is negative).\n\
@item 'stop'\n\
An error message is printed and the program is stopped.\n\
@end table\n\
\n\
Compute device:\n\
\n\
@table @code\n\
@item 'device'\n\
A string which can be one of 'cpu' or 'gpu'.\n\
@end table\n\
\n\
das is an oct-function that is a part of the DREAM Toolbox available at\n\
@url{https://github.com/frli8848/DREAM}.\n\
\n\
Copyright @copyright{} 2008-2023 Fredrik Lingvall.\n\
@end deftypefn")
{
  octave_value_list oct_retval;

  int nrhs = args.length ();

  ArgParser ap;

  //
  // Check for proper number of arguments
  //

  if ((nrhs < 8) && (nrhs > 10)) {
    error("das requires 8 to 10 input arguments!");
    return oct_retval;
  }

  if (nlhs > 2) {
    error("Too many output arguments for das!");
    return oct_retval;
  }

  //
  // Check if we are using single or double precision processsing
  //

  // Args 1, 2, 3, 4, and 6 must have the same datatype (float or double).
  bool use_float = args(0).is_single_type();
  if (use_float) {
    if (args(1).is_single_type() != use_float ||
        args(2).is_single_type() != use_float ||
        args(3).is_single_type() != use_float ||
        args(5).is_single_type() != use_float) {
      error("First arg is single precision but one of arg 2, 3, 4, or 6 is not!");
      return oct_retval;
    }
  } else {
    if (!(args(1).is_double_type() &&
          args(2).is_double_type() &&
          args(3).is_double_type() &&
          args(5).is_double_type()) ) {
      error("First arg is double precision but one of arg 2, 3, 4, or 6 is not!");
      return oct_retval;
    }
  }

  //
  // Data
  //

  dream_idx_type a_scan_len = mxGetM(0); // A-scan length
  dream_idx_type num_a_scans = mxGetN(0);

  float *Yf = nullptr;
  double *Yd = nullptr;

  if (use_float) {
    const FloatMatrix tmp0f= args(0).float_matrix_value();
    Yf = (float*) tmp0f.data();
  } else {
    const Matrix tmp0d = args(0).matrix_value();
    Yd = (double*) tmp0d.data();
  }

  //
  // Transmit array
  //

  if (mxGetN(1) != 3) {
    error("Argument 2 must be a (num transmit elements) x 3 matrix!");
    return oct_retval;
  }

  dream_idx_type num_t_elements = mxGetM(1);

  const float *Gt_f=nullptr;
  const double *Gt_d=nullptr;

  if (use_float) {
    const FloatMatrix tmp1f = args(1).float_matrix_value();
    Gt_f = tmp1f.data();
  } else {
    const Matrix tmp1d = args(1).matrix_value();
    Gt_d = tmp1d.data();
  }

  //
  // Recieve array
  //

  dream_idx_type num_r_elements = 0;
  const float *Gr_f = nullptr;
  const double *Gr_d = nullptr;

  if (mxGetM(2) != 0) { // Check if we do SAFT or TFM

    if (mxGetN(2) != 3) {
      error("Argument 3 must be a (num recieve elements) x 3 matrix!");
      return oct_retval;
    }

    num_r_elements = mxGetM(2);

    if (use_float) {
      const FloatMatrix tmp2f = args(2).float_matrix_value();
      Gr_f = tmp2f.data();
    } else {
      const Matrix tmp2d = args(2).matrix_value();
      Gr_d = tmp2d.data();
    }

  }

  //
  // Observation point.
  //

  // Check that arg (number of observation points) x 3 matrix.
  if (mxGetN(3) != 3) {
    error("Argument 4 must be a (number of observation points) x 3 matrix!");
    return oct_retval;
  }

  dream_idx_type No = mxGetM(3); // Number of observation points.

  const float *Ro_f = nullptr;
  const double *Ro_d = nullptr;

  if (use_float) {
    const FloatMatrix tmp3f= args(3).float_matrix_value();
    Ro_f = tmp3f.data();
  } else {
    const Matrix tmp3d = args(3).matrix_value();
    Ro_d = tmp3d.data();
  }

  //
  // Temporal sampling parameter.
  //

  if (!(mxGetM(4) == 1 && mxGetN(4) == 1) ) {
    error("Argument 5 must be a scalar!");
    return oct_retval;
  }

  const Matrix tmp4 = args(4).matrix_value();
  double *s_par = (double*) tmp4.data();
  double dt = s_par[0]; // Temporal discretization size (= 1/sampling freq).

  //
  // Start point of impulse response vector ([us]).
  //

  if (!ap.check_delay("das", args, 5, No)) {
    return oct_retval;
  }


  DelayType delay_type = DelayType::single;  // delay is a scalar.
  if (mxGetM(5) * mxGetN(5) != 1) {
    delay_type = DelayType::multiple; // delay is a vector.
  }

  const float *delay_f = nullptr;
  const double *delay_d = nullptr;

  if (use_float) {
    const FloatMatrix tmp5f = args(5).float_matrix_value();
    delay_f = tmp5f.data();
  } else {
    const Matrix tmp5d = args(5).matrix_value();
    delay_d =tmp5d.data();
  }

  //
  // Material parameter
  //

  // Check that arg 7 is a scalar.
  if (!(mxGetM(6)==1 && mxGetN(6)==1)) {
    error("Argument 7 must be a scalar!");
    return oct_retval;
  }

  // Sound speed.
  const Matrix tmp6 = args(6).matrix_value();
  double *m_par = (double*) tmp6.data();
  double cp = m_par[0]; // Sound speed.

  //
  // DAS method
  //

  DASType das_type;
  if (!ap.parse_das_arg("das", args, 7, das_type)) {
    return oct_retval;
  }

  //
  // Error reporting.
  //

  ErrorLevel err=ErrorLevel::none, err_level=ErrorLevel::stop;

  if (nrhs >= 9) {
    if (!ap.parse_error_arg("das", args, 8, err_level)) {
      return oct_retval;
    }
  } else {
    err_level = ErrorLevel::stop; // Default.
  }

  //
  // Compute device
  //

  std::string device;

  if (nrhs == 10) {

    if (!mxIsChar(9)) {
      error("Argument 10 must be a string");
      return oct_retval;
    }

    device = args(9).string_value();
  }

  //
  // Init DAS and output arg.
  //

  Matrix Im_mat_d;
  FloatMatrix Im_mat_f;
  float *Im_f = nullptr;
  double *Im_d = nullptr;

  bool init_das = true;

  if (use_float) {

    // Create an output matrix for the impulse response.
    Im_mat_f = FloatMatrix(No,1);
    Im_f = (float*) Im_mat_f.data();

    //
    // Single precision DAS
    //

    if (das_f) { // das (float) object exist - check if we can reuse previous das init
      if (!das_f->das_setup_has_changed(das_type, a_scan_len, No, num_t_elements, num_r_elements)) {
        init_das = false;
      }
    }

    if (init_das) {

      if (das_d) {
        das_d = nullptr; // Release the double obejct if it exist to free, in particular, GPU resources (call destructor).
      }

      try {
        das_f = std::make_unique<DAS<float>>(das_type, a_scan_len, No, num_t_elements, num_r_elements);
      }

      catch (std::runtime_error &err) {
        std::cout << err.what();
        return oct_retval;
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
    Im_mat_d = Matrix(No,1);
    Im_d = (double*) Im_mat_d.data();

    if (das_d) { // das (double) object exist - check if we can reuse previous das init
      if (!das_d->das_setup_has_changed(das_type, a_scan_len, No, num_t_elements, num_r_elements)) {
        init_das = false;
      }
    }

    if (init_das) {

      if (das_f) {
        das_f = nullptr; // Release the float object if it exist to free, in particular, GPU resources (call destructor).
      }

      try {
        das_d = std::make_unique<DAS<double>>(das_type, a_scan_len, No, num_t_elements, num_r_elements);
      }

      catch (std::runtime_error &err) {
        std::cout << err.what();
        return oct_retval;
      }
    }

    das_d->set_running();

    // Register signal handler.
    std::signal(SIGABRT, DAS<double>::abort);
  }

#ifdef USE_OPENCL

  // Check if we should use the GPU
  if (device == "gpu" && num_r_elements > 0) { // SAFT is most likely fast enough on the CPU.

    if (use_float) { // Single precision
      das_f->cl_das(Yf, Ro_f, Gt_f, Gr_f, dt, delay_f[0], cp, Im_f);
    } else { // Double precision
      das_d->cl_das(Yd, Ro_d, Gt_d, Gr_d, dt, delay_d[0], cp, Im_d);
    }

  } else { // Otherwise use the cpu

#endif

    if (device == "gpu") {
      std::cout << "Warning: Compute device set to 'gpu' but DREAM is build without OpenCL support!" << std::endl;
      std::cout << "Using the CPU backend!" << std::endl;
    }

    if (use_float) { // Single precision

      err = das_f->das(Yf, Ro_f, Gt_f, Gr_f, (float) dt, delay_type, delay_f, (float) cp, Im_f, err_level);
      if (!das_f->is_running()) {
        error("CTRL-C pressed!\n"); // Bail out.
        return oct_retval;
      }

    } else { // Double precision

      err = das_d->das(Yd, Ro_d, Gt_d, Gr_d, dt, delay_type, delay_d, cp, Im_d, err_level);
      if (!das_d->is_running()) {
        error("CTRL-C pressed!\n"); // Bail out.
        return oct_retval;
      }

    }

#ifdef USE_OPENCL
  }
#endif

  if (err == ErrorLevel::stop) {
    error("Error in DAS"); // Bail out if error.
    return oct_retval;
  }

  if (use_float) {
    oct_retval.append(Im_mat_f);
  } else {
    oct_retval.append(Im_mat_d);
  }

  // Return error.
  if (nlhs == 2) {
    Matrix err_mat(1, 1);
    double *err_p = (double*) err_mat.data();
    err_p[0] = (double) err;
    oct_retval.append(err_mat);
  }

  return oct_retval;
}
