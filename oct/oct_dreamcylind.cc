/***
*
* Copyright (C) 2006,2007,2008,2009,2012,2014,2015,2016,2021,2023,2024 Fredrik Lingvall
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

#include "dreamcylind.h"
#include "arg_parser.h"

//
// Octave headers.
//

#include <octave/oct.h>

/***
 *
 *  oct_dreamcylind.cc - Octave (oct) gateway function for (parallel) dreamcylind.
 *
 ***/

DEFUN_DLD (dreamcylind, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} [H,err] = dreamcylind(Ro,geom_par,s_par,delay,m_par,err_level)\n \
\n\
DREAMCYLIND - Computes the spatial impulse response\n\
for a concave (focused) cylindrical transducer using parallel \n\
processing (using threads).\n\
\n\
Observation point(s) ([mm]):\n\
\n\
@table @code\n\
@item Ro\n\
An N x 3 matrix, Ro = [xo1 yo1 zo2; xo2 yo2 zo2; ... xoN yoN zoN]; where N is the number of observation points.\n\
@end table\n\
\n\
Geometrical parameters: geom_par = [a b Rcurv];\n\
\n\
@table @code\n\
@item a\n\
x-size  of the transducer [mm].\n\
@item b\n\
y-size  of the transducer [mm].\n\
@item Rcurv\n\
Radius of the curvature [mm]. If Rcurv > 0 then it is focsued/concave and if Rcurv < 0 it is convex/defocused.\n\
\n\
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
Error Handling: err_level;\n\
err_level is an optional text string parameter for controlling the error behavior, options are:\n\
\n\
@table @code\n\
@item 'ignore'\n\
An error is ignored (no error message is printed and the program is not stopped) but the err output \n\
argument is non-zero if an error occured.\n\
@item 'warn'\n\
An error message is printed but the program in not stopped (and err is non-zero).\n\
@item 'stop'\n\
An error message is printed and the program is stopped.\n\
@end table\n\
\n\
dreamcylind is an oct-function that is a part of the DREAM Toolbox available at\n\
@url{https://github.com/frli8848/DREAM}.\n\
\n\
Copyright @copyright{} 2006-2023 Fredrik Lingvall.\n\
@seealso {dreamsphere, dreamcirc_f}\n\
@end deftypefn")
{
  octave_value_list oct_retval;

  int nrhs = args.length ();

  ArgParser ap;

  // Check for proper number of arguments

  if (!ap.check_arg_in("dreamcylind", nrhs, 5, 6)) {
    return oct_retval;
  }

  if (!ap.check_arg_out("dreamcylind", nlhs, 0, 2)) {
    return oct_retval;
  }

  //
  // Observation point.
  //

  if (!ap.check_obs_points("dreamcylind", args, 0)) {
    return oct_retval;
  }

  dream_idx_type No = mxGetM(0); // Number of observation points.
  const Matrix tmp0 = args(0).matrix_value();
  double *Ro = (double*) tmp0.data();

  //
  // Transducer geometry
  //

  double a=0.0, b=0.0, Rcurv=0.0;
  if (!ap.parse_geometry("dreamcylind", args, 1, 3, a, b, Rcurv)) {
    return oct_retval;
  }

  //
  // Temporal and spatial sampling parameters.
  //

  double dx=0.0, dy=0.0, dt=0.0;
  dream_idx_type nt=0;
  if (!ap.parse_sampling("dreamcylind", args, 2, 4, dx, dy, dt, nt)) {
    return oct_retval;
  }

  //
  // Start point of impulse response vector ([us]).
  //

  if (!ap.check_delay("dreamcylind", args, 3, No)) {
    return oct_retval;
  }

  const Matrix tmp3 = args(3).matrix_value();
  double *delay = (double*) tmp3.data();

  DelayType delay_type = DelayType::single;  // delay is a scalar.
  if (mxGetM(3) * mxGetN(3) != 1) {
    delay_type = DelayType::multiple; // delay is a vector.
  }

  //
  // Material parameters
  //

  double v=1.0, cp=1000.0, alpha=0.0;
  if (!ap.parse_material("dreamcylind", args, 4,  v, cp, alpha)) {
    return oct_retval;
  }

  //
  // Error reporting.
  //

  ErrorLevel err_level=ErrorLevel::stop;

  if (nrhs == 6) {
    if (!ap.parse_error_arg("dreamcylind", args, 5, err_level)) {
      return oct_retval;
    }
  } else {
    err_level = ErrorLevel::stop; // Default.
  }

  // Create an output matrix for the impulse response
  Matrix h_mat(nt, No);
  double *h = (double*) h_mat.data();

  SIRData hsir(h, nt, No);
  hsir.clear();

  Cylind cylind;

  // Register signal handler.
  std::signal(SIGINT, Cylind::abort);

  //
  // Call the DREAM subroutine.
  //

  SIRError err = cylind.dreamcylind(alpha,
                                    Ro, No,
                                    a,  b,  Rcurv,
                                    dx,  dy,  dt,  nt,
                                    delay_type,  delay,
                                    v,  cp,
                                    h, err_level);

  if (!cylind.is_running()) {
    if (err == SIRError::out_of_bounds) {
      error("SIR out-of-bounds!\n"); // Bail out.
    } else {
      error("CTRL-C pressed!\n"); // Bail out.
    }

    return oct_retval;
  }

  if (err == SIRError::warn_out_of_bounds) {
    std::cout << "Warning: SIR out-of-bounds!" << std::endl;
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
