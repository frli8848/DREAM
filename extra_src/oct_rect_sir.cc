/***
*
* Copyright (C) 2008,2009,2012,2015,2016,2021,2020,2023 Fredrik Lingvall
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

#include "rect_sir.h"
#include "arg_parser.h"

#include <octave/oct.h>

/***
 *
 *  rect_sir - Octave (oct) gateway function for RECT_SIR.
 *
 ***/

DEFUN_DLD (rect_sir, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {}  [Y] = rect_sir(Ro,geom_par,s_par,delay,m_par,device).\n \
\n\
RECT_SIR Computes the time-continous (analytic) spatial impulse response(s) for a\n\
rectangular transducer.\n\
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
@item b\n\
Dimensions of the transducer [mm].\n\
@end table\n\
\n\
Sampling parameters: s_par = [dt nt]; \n\
\n\
@table @code\n\
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
Material parameters: m_par = [v cp];\n\
\n\
@table @code\n\
@item v\n\
Normal velocity [m/s].\n\
@item cp\n\
Sound velocity [m/s].\n\
@end table\n\
\n\
Compute device:\n\
\n\
@table @code\n\
@item 'device'\n\
A string which can be one of 'cpu' or 'gpu'.\n\
@end table\n\
\n\
rect_sir is an oct-function that is a part of the DREAM Toolbox available at\n\
@url{https://github.com/frli8848/DREAM}.\n\
\n\
Copyright @copyright{} 2008-2023 Fredrik Lingvall.\n\
@seealso {dreamrect,circ_sir}\n\
@end deftypefn")
{
  std::string device;
  octave_value_list oct_retval;

  int nrhs = args.length ();

  ArgParser ap;

  // Check for proper number of arguments

  if (!ap.check_arg_in("rect_sir", nrhs, 5, 6)) {
    return oct_retval;
  }

  if (!ap.check_arg_out("rect_sir", nlhs, 0, 1)) {
    return oct_retval;
  }

  //
  // Observation point.
  //

  if (!ap.check_obs_points("rect_sir", args, 0)) {
    return oct_retval;
  }

  dream_idx_type No = mxGetM(0); // Number of observation points.
  const Matrix tmp0 = args(0).matrix_value();
  double *Ro = (double*) tmp0.data();

  //
  // Transducer geometry
  //

  double a=0.0, b=0.0, dummy=0.0;
  if (!ap.parse_geometry("rect_sir", args, 1, 2, a, b, dummy)) {
    return oct_retval;
  }

  //
  // Temporal and spatial sampling parameters.
  //

  // Check that arg 3 is a 2 element vector
  if (!((mxGetM(2) == 2 && mxGetN(2) == 1) ||
        (mxGetM(2) == 1 && mxGetN(2) == 2))) {
    error("Argument 3 must be a vector of length 2!");
    return oct_retval;
  }

  const Matrix tmp2 = args(2).matrix_value();
  double *s_par = (double*) tmp2.data();
  double dt = s_par[0]; // Temporal discretization size (= 1/sampling freq).
  dream_idx_type nt = (dream_idx_type) s_par[1];	// Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //

  if (!ap.check_delay("rect_sir", args, 3, No)) {
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

  // Check that arg 5 is a 2 element vector.
  if (!((mxGetM(4) == 2 && mxGetN(4) == 1) ||
        (mxGetM(4) == 1 && mxGetN(4) == 2))) {
    error("Argument 5 must be a vector of length 2!");
    return oct_retval;
  }

  const Matrix tmp4 = args(4).matrix_value();
  double *m_par = (double*) tmp4.data();
  double v  = m_par[0];      // Normal velocity of transducer surface.
  double cp = m_par[1];      // Sound speed.

  /*
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
  */

  //
  // Compute device
  //

  if (nrhs == 6) {

    if (!mxIsChar(5)) {
      error("Argument 6 must be a string");
      return oct_retval;
    }

    device = args(5).string_value();
  }

  //
  // Create an output matrix for the impulse response(s).
  //

  Matrix h_mat(nt, No);
  double *h = (double*) h_mat.data();

  SIRData hsir(h, nt, No);
  hsir.clear();

  RectSir rect_sir;

  // Register signal handler.
  std::signal(SIGABRT, RectSir::abort);

  //
  // Call the rect_sir subroutine.
  //

  // Check if we should use the GPU

#ifdef USE_OPENCL
  if (device == "gpu") {
    rect_sir.cl_rect_sir(Ro, No, a, b, dt, nt, delay[0], v, cp, h);

  } else { // Otherwise use the cpu
#endif

    ErrorLevel err = rect_sir.rect_sir(Ro, No,
                                       a, b,
                                       dt, nt,
                                       delay_type, delay,
                                       v, cp,
                                       h);

    if (!rect_sir.is_running()) {
      dream_err_msg("CTRL-C pressed!\n"); // Bail out.
    }

    //
    // Check for Error. FIXME: Do we need to here?
    //

    // NB. The GPU and analytical SIR code do not return any error codes.
    if (err == ErrorLevel::stop) {
      error("Error in rect_sir"); // Bail out if error.
      return oct_retval;
    }

#ifdef USE_OPENCL
  }
#endif

  oct_retval.append(h_mat);

  /*
  if (nlhs == 2) {
    Matrix err_mat(1, 1);
    double *err_p = (double*) err_mat.data();
    err_p[0] = (double) err;
    oct_retval.append(err_mat);
  }
  */

  return oct_retval;
}
