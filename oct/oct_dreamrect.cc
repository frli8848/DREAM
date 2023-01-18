/***
*
* Copyright (C) 2006,2007,2008,2009,2012,2014,2015,2016,2021,2023 Fredrik Lingvall
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

#include "dreamrect.h"
#include "arg_parser.h"

#include <octave/oct.h>

/***
 *
 *  oct_dreamrect.cc - Octave (oct) gateway function for dreamrect.
 *
 ***/

DEFUN_DLD (dreamrect, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} [H,err] = dreamrect(Ro,geom_par,s_par,delay,m_par,err_level,device)\n\
\n\
DREAMRECT - Computes the spatial impulse response\n\
for a rectangular transducer using parallel processing\n\
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
x-size  of the transducer [mm].\n\
@item b\n\
y-size  of the transducer [mm].\n\
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
argument is negative if an error occured.\n\
@item 'warn'\n\
An error message is printed but the program in not stopped (and err is negative).\n\
@item 'stop'\n\
An error message is printed and the program is stopped.\n\
@item 'device'\n\
A string which can be one of 'cpu' or 'gpu'.\n\
@end table\n\
\n\
dreamrect is an oct-function that is a part of the DREAM Toolbox available at\n\
@url{https://github.com/frli8848/DREAM}.\n\
\n\
Copyright @copyright{} 2006-2023 Fredrik Lingvall.\n\
@seealso {dreamrect}\n\
@end deftypefn")
{
  std::string device;

  octave_value_list oct_retval;

  int nrhs = args.length();

  ArgParser ap;

  // Check for proper number of arguments

  if (!ap.check_arg_in("dreamrect", nrhs, 5, 7)) {
    return oct_retval;
  }

  if (!ap.check_arg_out("dreamrect", nlhs, 0, 2)) {
    return oct_retval;
  }

  //
  // Observation point.
  //

  if (!ap.check_obs_points("dreamrect", args, 0)) {
    return oct_retval;
  }

  dream_idx_type no = mxGetM(0); // Number of observation points.
  const Matrix tmp0 = args(0).matrix_value();
  double *ro = (double*) tmp0.data();

  //
  // Transducer geometry
  //

  double a=0.0, b=0.0, dummy=0.0;
  if (!ap.parse_geometry("dreamrect", args, 1, 2, a, b, dummy)) {
    return oct_retval;
  }

  //
  // Temporal and spatial sampling parameters.
  //

  double dx=0.0, dy=0.0, dt=0.0;
  dream_idx_type nt=0;
  if (!ap.parse_sampling("dreamrect", args, 2, 4, dx, dy, dt, nt)) {
    return oct_retval;
  }

  //
  // Start point of impulse response vector ([us]).
  //

  if (!ap.check_delay("dreamrect", args, 3, no)) {
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
  if (!ap.parse_material("dreamrect", args, 4, v, cp, alpha)) {
    return oct_retval;
  }

  //
  // Error reporting.
  //

  ErrorLevel err=ErrorLevel::none, err_level=ErrorLevel::stop;

  if (nrhs >= 6) {
    if (!ap.parse_error_arg("dreamrect", args, 5, err_level)) {
      return oct_retval;
    }
  } else {
    err_level = ErrorLevel::stop; // Default.
  }

  //
  // Compute device
  //

  if (nrhs == 7) {

    if (!mxIsChar(6)) {
      error("Argument 7 must be a string");
      return oct_retval;
    }

    device = args(6).string_value();
  }

  // Create an output matrix for the impulse response
  Matrix h_mat(nt, no);
  double *h = (double*) h_mat.data();

  SIRData hsir(h, nt, no);
  hsir.clear();

  Rect rect;

  // Register signal handler.
  std::signal(SIGABRT, Rect::abort);

  //
  // Call the DREAM subroutine.
  //

#ifdef USE_OPENCL

  // Check if we should use the GPU
  if (device == "gpu") {
    rect.cl_dreamrect(ro, no, a, b,  dx, dy, dt, nt,  delay[0], v, cp, h);

  } else { // Otherwise use the cpu
#endif

    err = rect.dreamrect(alpha,
                         ro,  no,
                         a,  b,
                         dx,  dy,  dt, nt,
                         delay_type, delay,
                         v,cp,
                         h, err_level);

    if (!rect.is_running()) {
      error("CTRL-C pressed!\n"); // Bail out.
      return oct_retval;
    }

#ifdef USE_OPENCL
  }
#endif

  //
  // Check for Error.
  //

  // NB. The GPU code do not return any error codes (yet).
  if (err == ErrorLevel::stop) {
    error("Error in dreamrect"); // Bail out if error.
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
