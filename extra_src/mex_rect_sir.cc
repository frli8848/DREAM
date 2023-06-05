/***
*
* Copyright (C) 2004,2006,2007,2008,2009,2014,2015,2016,2019,2021,2023 Fredrik Lingvall
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

#include "mex.h"

/***
 *
 * Matlab (MEX) gateway function for RECT_SIR.
 *
 ***/

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  ArgParser ap;

  // Check for proper number of arguments
  ap.check_arg_in("rect_sir", nrhs, 5, 6); // NB. 6:th arg is for 'gpu' (OpenCL) which is not implemented here yet.
  ap.check_arg_out("rect_sir", nlhs, 0, 1);

  //
  // Observation point.
  //

  ap.check_obs_points("rect_sir", prhs, 0);
  dream_idx_type No = mxGetM(prhs[0]); // Number of observation points.
  double *Ro = mxGetPr(prhs[0]);

  //
  // Transducer geometry
  //

  double a=0.0, b=0.0, dummy=0.0;
  ap.parse_geometry("rect_sir", prhs, 1, 2, a, b, dummy);

  //
  // Temporal and spatial sampling parameters.
  //

  // Check that arg 3 is a 2 element vector
  if (!((mxGetM(prhs[2]) == 2 && mxGetN(prhs[2]) == 1) ||
        (mxGetM(prhs[2]) == 1 && mxGetN(prhs[2]) == 2))) {
    dream_err_msg("Argument 3 must be a vector of length 2!");
  }

  double *s_par = mxGetPr(prhs[2]);
  double dt = s_par[0];		// Temporal discretization size (= 1/sampling freq).
  double nt = (size_t) s_par[1];	// Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //

  ap.check_delay("rect_sir", prhs, 3, No);
  double *delay = mxGetPr(prhs[3]);

  DelayType delay_type = DelayType::single;  // delay is a scalar.
  if (mxGetM(prhs[3]) * mxGetN(prhs[3]) != 1) {
    delay_type = DelayType::multiple; // delay is a vector.
  }

  //
  // Material parameters
  //

  // Check that arg 5 is a 2 element vector.
  if (!((mxGetM(prhs[4]) == 2 && mxGetN(prhs[4]) == 1) ||
        (mxGetM(prhs[4]) == 1 && mxGetN(prhs[4]) == 2))) {
    dream_err_msg("Argument 5 must be a vector of length 2!");
  }

  double *m_par = mxGetPr(prhs[4]);
  double v = m_par[0];       // Normal velocity of transducer surface.
  double cp = m_par[1];      // Sound speed.

  //
  // Create an output matrix for the impulse response(s).
  //

  plhs[0] = mxCreateDoubleMatrix(nt,No,mxREAL);
  double *h = mxGetPr(plhs[0]);

  SIRData hsir(h, nt, No);
  hsir.clear();

  RectSir rect_sir;

  // Register signal handler.
  std::signal(SIGABRT, RectSir::abort);

  //
  // Call the analytic rect_sir subroutine.
  //

  ErrorLevel err = rect_sir.rect_sir(Ro, No,
                                     a, b,
                                     dt, nt,
                                     delay_type, delay,
                                     v, cp,
                                     h);

  if (!rect_sir.is_running()) {
    dream_err_msg("CTRL-C pressed!\n"); // Bail out.
  }

  // FIXME. The analytical SIR code do not return any error codes.
  if (err == ErrorLevel::stop) {
    dream_err_msg("Error in rect_sir"); // Bail out if error.
  }

  return;
}
