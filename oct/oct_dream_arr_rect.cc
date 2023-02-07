/***
*
* Copyright (C) 2006,2007,2008,2009,2012,2014,2015,2016,2021,2021,2023 Fredrik Lingvall
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

#include "dream_arr_rect.h"
#include "arg_parser.h"

//
// Octave headers.
//

#include <octave/oct.h>

/***
 *
 *  Octave (oct) gateway function for (parallel) dream_arr_rect.
 *
 ***/

DEFUN_DLD (dream_arr_rect, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} [H,err] = dream_arr_rect(Ro,geom_par,G,s_par,delay,m_par,foc_met,...\n\
                focal,steer_met,steer_par,apod_met,apod,win_par,err_level);\n\
\n\
DREAM_ARR_RECT - Computes the spatial impulse response\n\
for an array transducer with rectangular elements using parallel processing\n\
(using threads).\n\
\n\
Observation point(s) ([mm]):\n\
\n\
@table @code\n\
@item Ro\n\
An N x 3 matrix, Ro = [xo1 yo1 zo2; xo2 yo2 zo2; ... xoN yoN zoN]; where N is the number of observation points.\n\
@end table\n\
\n\
Geometrical parameters: geom_par = [a b];\n\
\n\
@table @code\n\
@item a\n\
Element width in x-direction [mm].\n\
@item b\n\
Element width in y-direction [mm].\n\
@end table\n\
\n\
Array grid parameter:\n\
\n\
@table @code\n\
@item G\n\
An L x 3 grid function matrix. Column 1 contain the x-postions \n\
of the elements, column 2 the y-positions, and column 3\n\
the z-positions, and where L is the number of elements.\n\
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
Focusing method, options are: 'off', 'x', 'y', 'xy', 'x+y', and 'ud'.\n\
@item  focal\n\
Focal distance [mm]. If foc_met = 'ud' (user defined) then focal is a vector of focusing delays.\n\
@end table\n\
\n\
Beam steering parameters: steer_met and steer_par:\n\
\n\
@table @code\n\
@item steer_met\n\
Beam steering method, options are: 'off', 'x', 'y', and 'xy'.\n\
@item  steer_par =  [theta phi];\n\
theta [deg] is the x-direction steer angle and \n\
phi [deg] is the y-direction steer angle.\n\
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
A warning message is printed but the program in not stopped.\n\
@item 'stop'\n\
An error message is printed and the program is stopped.\n\
@end table\n\
\n\
dream_arr_rect is an oct-function that is a part of the DREAM Toolbox available at\n\
@url{https://github.com/frli8848/DREAM}.\n\
\n\
Copyright @copyright{} 2006-2023 Fredrik Lingvall.\n\
@seealso {dreamrect}\n\
@end deftypefn")
{
  octave_value_list oct_retval;

  int nrhs = args.length();

  ArgParser ap;

  // Check for proper number of arguments

  if (!ap.check_arg_in("dream_arr_rect", nrhs, 13, 14)) {
    return oct_retval;
  }

  if (!ap.check_arg_out("dream_arr_rect", nlhs, 0, 2)) {
    return oct_retval;
  }

  //
  // Observation points.
  //

  if (!ap.check_obs_points("dream_arr_rect", args, 0)) {
    return oct_retval;
  }

  dream_idx_type no = mxGetM(0); // Number of observation points.
  const Matrix tmp0 = args(0).matrix_value();
  double *ro = (double*) tmp0.data();

  //
  // Transducer geometry
  //

  double a=0.0, b=0.0, dummy=0.0;
  if (!ap.parse_geometry("dream_arr_rect", args, 1, 2, a, b, dummy)) {
    return oct_retval;
  }

  //
  // Grid function (position vectors of the elements).
  //

  if (!ap.check_array("dream_arr_rect", args, 2)) {
    return oct_retval;
  }

  dream_idx_type num_elements = (dream_idx_type) mxGetM(2); // Number of elementents in the array.
  const Matrix tmp2 = args(2).matrix_value();
  double *G = (double*) tmp2.data(); // First column in the matrix.
  //gy    = gx + num_elements;		// Second column in the matrix.
  //gz    = gy + num_elements;		// Third column in the matrix.

  //
  // Temporal and spatial sampling parameters.
  //

  double dx=0.0, dy=0.0, dt=0.0;
  dream_idx_type nt=0;
  if (!ap.parse_sampling("dream_arr_rect", args, 3, 4, dx, dy, dt, nt)) {
    return oct_retval;
  }

  //
  // Start point of impulse response vector ([us]).
  //

  if (!ap.check_delay("dream_arr_rect", args, 4, no)) {
    return oct_retval;
  }

  const Matrix tmp4 = args(4).matrix_value();
  double *delay = (double*) tmp4.data();

  DelayType delay_type = DelayType::single;  // delay is a scalar.
  if (mxGetM(4) * mxGetN(4) != 1) {
    delay_type = DelayType::multiple; // delay is a vector.
  }

  //
  // Material parameters
  //

  double v=1.0, cp=1000.0, alpha=0.0;
  if (!ap.parse_material("dream_arr_rect", args, 5, v, cp, alpha)) {
    return oct_retval;
  }

  //
  // Focusing parameters.
  //

  FocusMet foc_met=FocusMet::none;

  // Allocate space for the user defined focusing delays
  std::unique_ptr<double[]> focal = std::make_unique<double[]>(num_elements);

  if (nrhs >= 7) {
    if (!ap.parse_focus_args("dream_arr_rect", args, 6, foc_met, focal.get())) {
      return oct_retval;
    }
  } else {
    foc_met = FocusMet::none;
  }

  //
  // Beam steering.
  //

  SteerMet steer_met=SteerMet::none;

  double theta=0.0, phi=0.0;
  if (nrhs >= 9) {
    if (!ap.parse_steer_args("dream_arr_rect", args, 8, steer_met, theta, phi)) {
      return oct_retval;
    }
  } else {
    steer_met = SteerMet::none;
  }

  //
  // Apodization.
  //

  bool    do_apod=false;
  ApodMet apod_met=ApodMet::gauss;

  // Allocate space for the user defined apodization weights
  std::unique_ptr<double[]> apod = std::make_unique<double[]>(num_elements);

  double apod_par=0.0;
  if (nrhs >= 11) {
    if (!ap.parse_apod_args("dream_arr_rect", args, 10, num_elements,
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

  if (nrhs == 14) {
    if (!ap.parse_error_arg("dream_arr_rect", args, 13, err_level)) {
      return oct_retval;
    }
  } else {
    err_level = ErrorLevel::stop; // Default.
  }

  // Create an output matrix for the impulse response
  Matrix h_mat(nt, no);
  double *h = (double*) h_mat.data();

  SIRData hsir(h, nt, no);
  hsir.clear();

  ArrRect arr_rect;

  // Register signal handler.
  std::signal(SIGINT, ArrRect::abort);

  //
  // Call the DREAM subroutine.
  //

  err = arr_rect.dream_arr_rect(alpha,
                                ro, no,
                                a, b,
                                dx, dy,  dt, nt,
                                delay_type, delay,
                                v, cp,
                                num_elements, G,
                                foc_met,  focal.get(),
                                steer_met, theta, phi,
                                apod.get(), do_apod, apod_met, apod_par,
                                h, err_level);

  if (!arr_rect.is_running()) {
    error("CTRL-C pressed!\n"); // Bail out.
    return oct_retval;
  }

  //
  // Check for Error.
  //

  if (err == ErrorLevel::stop) {
    error("Error in dream_arr_rect"); // Bail out if error.
    return oct_retval;
  }

  oct_retval.append(h_mat);

  //
  // Return error.
  //

  if (nlhs == 2) {
    Matrix err_mat(1, 1);
    double *err_p = (double*) err_mat.data();
    err_p[0] = (double) err;
    oct_retval.append(err_mat);
  }

  return oct_retval;
}
