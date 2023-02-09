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

/***
 *
 * das - Octave (oct) gateway function for DAS (delay-and-sum).
 *
 ***/

DEFUN_DLD (das, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {}  [Im] = das(Y,Gt,Gr,Ro,s_par,delay,m_par,method,err_level,device).\n \
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

  // Check for proper number of arguments

  if ((nrhs < 8) && (nrhs > 10)) {
    error("das requires 8 to 10 input arguments!");
    return oct_retval;
  }

  if (nlhs > 2) {
    error("Too many output arguments for das!");
    return oct_retval;
  }

  //
  // Data
  //

  dream_idx_type a_scan_len = mxGetM(0); // A-scan length
  dream_idx_type num_a_scans = mxGetN(0);
  const Matrix tmp0 = args(0).matrix_value();
  double *Y = (double*) tmp0.data();

  //
  // Transmit array
  //

  if (mxGetN(1) != 3) {
    error("Argument 2 must be a (num transmit elements) x 3 matrix!");
    return oct_retval;
  }

  dream_idx_type num_t_elements = mxGetM(1);
  const Matrix tmp1 = args(1).matrix_value();
  double *Gt = (double*) tmp1.data();

  //
  // Recieve array
  //

  dream_idx_type num_r_elements = 0;
  double *Gr = nullptr;

  if (mxGetM(2) != 0) { // Check if we do SAFT or TFM

    if (mxGetN(2) != 3) {
      error("Argument 3 must be a (num recieve elements) x 3 matrix!");
      return oct_retval;
    }

    num_r_elements = mxGetM(2);
    const Matrix tmp2 = args(2).matrix_value();
    Gr = (double*) tmp2.data();
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
  const Matrix tmp3 = args(3).matrix_value();
  double *Ro = (double*) tmp3.data();

  //
  // Temporal and spatial sampling parameters.
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

  const Matrix tmp5 = args(5).matrix_value();
  double *delay = (double*) tmp5.data();

  DelayType delay_type = DelayType::single;  // delay is a scalar.
  if (mxGetM(5) * mxGetN(5) != 1) {
    delay_type = DelayType::multiple; // delay is a vector.
  }

  //
  // Material parameters
  //

  // Check that arg 7 is a scalar.
  if (!(mxGetM(6)==1 && mxGetN(6)==1)) {
    error("Argument 7 must be a scalar!");
    return oct_retval;
  }

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

  // Create an output matrix for the impulse response.
  Matrix Im_mat(No,1);
  double *Im = (double*) Im_mat.data();

  std::unique_ptr<DAS> das;
  try {
    das = std::make_unique<DAS>(das_type, a_scan_len, No, num_t_elements, num_r_elements);
  }

  catch (std::runtime_error &err) {
    std::cout << err.what();
    return oct_retval;
  }

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

    if (!das->is_running()) {
      error("CTRL-C pressed!\n"); // Bail out.
      return oct_retval;
    }

#ifdef USE_OPENCL
  }
#endif

  if (err == ErrorLevel::stop) {
    error("Error in DAS"); // Bail out if error.
    return oct_retval;
  }

  oct_retval.append(Im_mat);

  // Return error.
  if (nlhs == 2) {
    Matrix err_mat(1, 1);
    double *err_p = (double*) err_mat.data();
    err_p[0] = (double) err;
    oct_retval.append(err_mat);
  }

  return oct_retval;
}
