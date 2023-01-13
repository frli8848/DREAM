/***
*
* Copyright (C) 2006,2007,2008,2009,2012,2015,2016,2021,2021,2023 Fredrik Lingvall
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

//
// Octave headers.
//

#include <octave/oct.h>

/***
 *
 *  Octave (oct) gateway function for dream_arr_annu
 *
 ***/

DEFUN_DLD (dream_arr_annu, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} [H,err] = dream_arr_annu(Ro,G,s_par,delay,m_par,foc_met,focal,...\n\
                                               apod_met,apod,win_par,err_level);\n\
\n\
 DREAM_ARR_ANNU - Computes the spatial impulse response\n\
for an annular array using parallel processing (using threads).\n\
\n\
Observation point(s) ([mm]):\n\
\n\
@table @code\n\
@item Ro\n\
An N x 3 matrix, Ro = [xo1 yo1 zo2; xo2 yo2 zo2; ... xoN yoN zoN]; where N is the number of observation points.\n\
@end table\n\
\n\
 Array grid parameter:\n\
@table @code\n\
@item G\n\
  Vector of annulus radii [mm].\n\
@end table\n\
\n\
Sampling parameters: s_par = [dx dy dt nt]; \n\
\n\
@table @code\n\
@item dx\n\
Spatial x-direction discretization size [mm].\n\
@item dy\n\
Spatial y-direction discretization size [mm].\n\
@item dt\n\
Temporal discretization period (= 1/sampling freq) [us].\n\
@item  nt\n\
Length of impulse response vector.\n\
@end table\n\
\n\
Start point of SIR:\n\
\n\
@table @code\n\
@item  delay\n\
Scalar delay for all observation points or a vector with individual delays for each observation point [us].\n\
@end table\n\
\n\
Material parameters: m_par = [v cp alpha];\n\
\n\
@table @code\n\
@item v\n\
Normal velocity [m/s].\n\
@item cp\n\
Sound velocity [m/s].\n\
@item alpha\n\
Attenuation coefficient [dB/(cm MHz)].\n\
\n\
@end table\n\
\n\
Focusing parameters: foc_met and focal:\n\
\n\
@table @code\n\
@item foc_met\n\
Focusing method, options are: 'on', 'off', and 'ud'.\n\
@item  focal\n\
Focal distance [mm]. If foc_met = 'ud' (user defined) then focal is a vector of focusing delays.\n\
@end table\n\
\n\
Apodization parameters: apod_met, apod, and win_par. The apod_met (apodization method) options are:\n\
\n\
@table @code\n\
\n\
@item 'off'\n\
No apodization.\n\
@item 'ud'\n\
User defined apodization.\n\
@item  'triangle'\n\
Triangle window.\n\
@item 'gauss'\n\
Gaussian (bell-shaped) window.\n\
@item 'raised'\n\
Raised cosine.\n\
@item 'simply'\n\
Simply supported.\n\
@item 'clamped'\n\
Clamped.\n\
@end table\n\
\n\
and the apod and win_par parameters are:\n\
\n\
@table @code\n\
@item apod\n\
Vector of apodiztion weights (used for the 'ud' option).\n\
@item win_par\n\
A scalar parameter for raised cosine and Gaussian apodization functions.\n\
@end table\n\
\n\
Error Handling: err_level;\n\
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
dream_arr_annu is an oct-function that is a part of the DREAM Toolbox available at\n\
@url{https://github.com/frli8848/DREAM}.\n\
\n\
Copyright @copyright{} 2006-2023 Fredrik Lingvall.\n\
@end deftypefn")
{
  FocusMet foc_met=FocusMet::none;
  bool    do_apod=false;
  ApodMet apod_met=ApodMet::gauss;
  double *h;

  octave_value_list oct_retval;

  int nrhs = args.length();

  ArgParser ap;

  // Check for proper number of arguments.

  if (!ap.check_arg_in("dream_arr_annu", nrhs, 10, 11)) {
    return oct_retval;
  }

  if (!ap.check_arg_out("dream_arr_annu", nlhs, 0, 2)) {
    return oct_retval;
  }

  //
  // Observation point.
  //

  if (!ap.check_obs_points("dream_arr_annu", args, 0)) {
    return oct_retval;
  }

  dream_idx_type no = mxGetM(0); // Number of observation points.
  const Matrix tmp0 = args(0).matrix_value();
  double *ro = (double*) tmp0.fortran_vec();

  //
  // Grid function (radii vectors of the elements).
  //

  if (!ap.check_array_annu("dream_arr_annu", args, 1)) {
    return oct_retval;
  }

  if ((mxGetM(1) > 1) & (mxGetN(1) > 1)) {
    error("Argument 2 must a vector (number of array elements)");
    return oct_retval;
  }

  dream_idx_type num_radii = (dream_idx_type) mxGetM(1)*mxGetN(1);
  dream_idx_type num_elements = (num_radii+1)/2;
  const Matrix tmp1 = args(1).matrix_value();
  double *gr = (double*) tmp1.fortran_vec(); // Vector of annular radi.

  //
  // Temporal and spatial sampling parameters.
  //

  double dx=0.0, dy=0.0, dt=0.0;
  dream_idx_type nt=0;
  if (!ap.parse_sampling("dream_arr_annu", args, 2, 4, dx, dy, dt, nt)) {
    return oct_retval;
  }

  //
  // Start point of impulse response vector ([us]).
  //

  if (!ap.check_delay("dream_arr_annu", args, 3, no)) {
    return oct_retval;
  }

  const Matrix tmp3 = args(3).matrix_value();
  double *delay = (double*) tmp3.fortran_vec();

  DelayType delay_type = DelayType::single;  // delay is a scalar.
  if (mxGetM(3) * mxGetN(3) != 1) {
    delay_type = DelayType::multiple; // delay is a vector.
  }

  //
  // Material parameters
  //

  double v=1.0, cp=1000.0, alpha=0.0;
  if (!ap.parse_material("dream_arr_annu", args, 4, v, cp,  alpha)) {
    return oct_retval;
  }

  //
  // Focusing parameters.
  //

  // Allocate space for the user defined focusing delays
  std::unique_ptr<double[]> focal = std::make_unique<double[]>(num_elements);

  if (nrhs >= 6) {
    if (!ap.parse_focus_args("dream_arr_annu", args, 5, foc_met, focal.get())) {
      return oct_retval;
    }
  } else {
    foc_met = FocusMet::none;
  }

  //
  // Apodization.
  //

  // Allocate space for the user defined apodization weights
  std::unique_ptr<double[]> apod = std::make_unique<double[]>(num_elements);

  double apod_par = 0.0;
  if (nrhs >= 8) {
    if (!ap.parse_apod_args("dream_arr_annu", args, 7, num_elements,
                            do_apod, apod.get(), apod_met, apod_par)) {
      return oct_retval;
    }
  } else {
    do_apod = false;
  }

  //
  // Error reporting.
  //

  ErrorLevel err=ErrorLevel::none, err_level=ErrorLevel::stop;

  if (nrhs == 11) {
    if (!ap.parse_error_arg("dream_arr_annu", args, 10, err_level)) {
      return oct_retval;
    }
  } else {
    err_level = ErrorLevel::stop; // Default.
  }

  // Create an output matrix for the impulse response
  Matrix h_mat(nt, no);
  h = h_mat.fortran_vec();

  SIRData hsir(h, nt, no);
  hsir.clear();

  ArrAnnu arr_annu;

  //
  // Call the DREAM subroutine.
  //

  err = arr_annu.dream_arr_annu(alpha,
                                ro, no,
                                dx, dy,  dt, nt,
                                delay_type, delay,
                                v, cp,
                                num_radii, gr,
                                foc_met, focal.get(),
                                apod.get(), do_apod, apod_met, apod_par,
                                h, err_level);

  if (!arr_annu.is_running()) {
    error("CTRL-C pressed!\n"); // Bail out.
    return oct_retval;
  }

  //
  // Check for Error.
  //

  if (err == ErrorLevel::stop) {
    error("Error in dream_arr_annu"); // Bail out if error.
    return oct_retval;
  }

  oct_retval.append(h_mat);

  //
  // Return error.
  //

  if (nlhs == 2) {
    Matrix err_mat(1, 1);
    double *err_p = err_mat.fortran_vec();
    err_p[0] = (double) err;
    oct_retval.append(err_mat);
  }

  return oct_retval;
}
