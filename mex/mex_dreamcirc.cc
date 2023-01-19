/***
*
* Copyright (C) 2003,2005,2006,2007,2008,2009,2014,2015,2019,2021,2023 Fredrik Lingvall
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

#include "dreamcirc.h"
#include "arg_parser.h"

#include "mex.h"

/***
 *
 * Matlab (MEX) gateway function for dreamcirc.
 *
 ***/

extern void _main();

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  ArgParser ap;

  // Check for proper number of arguments

  ap.check_arg_in("dreamcirc", nrhs, 5, 7); // NB. 7:th arg is for 'gpu' (OpenCL) which is not implemented here yet.
  ap.check_arg_out("dreamcirc", nlhs, 0, 2);

  //
  // Observation point.
  //

  ap.check_obs_points("dreamcirc", prhs, 0);
  dream_idx_type no = mxGetM(prhs[0]); // Number of observation points.
  double *ro = mxGetPr(prhs[0]);

  //
  // Transducer geometry
  //

  double R=0.0, dummy1=0.0, dummy2=0.0;
  ap.parse_geometry("dreamcirc", prhs, 1, 1, R, dummy1, dummy2);

  //
  // Temporal and spatial sampling parameters.
  //

  double dx=0.0, dy=0.0, dt=0.0;
  dream_idx_type nt=0;
  ap.parse_sampling("dreamcirc", prhs, 2, 4, dx, dy, dt, nt);

  //
  // Start point of impulse response vector ([us]).
  //

  ap.check_delay("dreamcirc", prhs, 3, no);
  double *delay = mxGetPr(prhs[3]);

  DelayType delay_type = DelayType::single;  // delay is a scalar.
  if (mxGetM(prhs[3]) * mxGetN(prhs[3]) != 1) {
    delay_type = DelayType::multiple; // delay is a vector.
  }

  //
  // Material parameters
  //

  double v=1.0, cp=1000.0, alpha=0.0;
  ap.parse_material("dreamcirc", prhs, 4,  v, cp, alpha);

  //
  // Error reporting.
  //

  ErrorLevel err=ErrorLevel::none, err_level=ErrorLevel::stop;

  if (nrhs == 6) {
    ap.parse_error_arg("dreamcirc", prhs, 5, err_level);
  } else {
    err_level = ErrorLevel::stop; // Default.
  }

  // Create an output matrix for the impulse response
  plhs[0] = mxCreateDoubleMatrix(nt,no,mxREAL);
  double *h = mxGetPr(plhs[0]);

  SIRData hsir(h, nt, no);
  hsir.clear();

  Circ circ;

  // Register signal handler.
  std::signal(SIGABRT, Circ::abort);

  //
  // Call the DREAM subroutine.
  //

  err = circ.dreamcirc(alpha,
                       ro, no,
                       R,
                       dx, dy, dt, nt,
                       delay_type, delay,
                       v,cp,
                       h, err_level);

  if (!circ.is_running()) {
    dream_err_msg("CTRL-C pressed!\n"); // Bail out.
  }

  //
  // Check for Error.
  //

  if (err == ErrorLevel::stop) {
    dream_err_msg(""); // Bail out if error.
  }

  //
  // Return error.
  //

  if (nlhs == 2) {
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    double *err_p =  mxGetPr(plhs[1]);
    err_p[0] = (double) err;
  }

  return;
}
