/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2014,2016,2021,2023 Fredrik Lingvall
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

#include <cmath>

#include "attenuation.h"
#include "arg_parser.h"

#include "mex.h"

/***
 *
 * Matlab (MEX) gateway function for dreamrect.
 *
 ***/

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  ArgParser ap;

  // Check for proper number of arguments

  ap.check_arg_in("dream_att", nrhs, 4, 4);
  ap.check_arg_out("dream_att", nlhs, 0, 1);

  //
  // Observation point.
  //

  ap.check_obs_points("dream_att", prhs, 0);
  dream_idx_type No = mxGetM(prhs[0]); // Number of observation points.
  double *Ro = mxGetPr(prhs[0]);

  //
  // Temporal and spatial sampling parameters.
  //

  // Check that arg 2 is a 4 element vector
  if (!((mxGetM(prhs[1])==2 && mxGetN(prhs[1])==1) || (mxGetM(prhs[1])==1 && mxGetN(prhs[1])==2))) {
    dream_err_msg("Argument 2 must be a vector of length 2!");
  }

  double *s_par = mxGetPr(prhs[1]);
  double dt = s_par[0];		// Temporal discretization size (= 1/sampling freq).
  dream_idx_type nt = (dream_idx_type) s_par[1];	// Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //

  ap.check_delay("dream_att", prhs, 2, No);
  double *delay = mxGetPr(prhs[2]);

  //
  // Material parameters
  //

  // Check that arg 4 is a 2 element vector.
  if (!((mxGetM(prhs[3])==2 && mxGetN(prhs[3])==1) || (mxGetM(prhs[3])==1 && mxGetN(prhs[3])==2))) {
    dream_err_msg("Argument 4 must be a vector of length 2!");
  }

  double *m_par = mxGetPr(prhs[3]);
  double cp    = m_par[0]; // Sound speed.
  double alpha = m_par[1]; // Attenuation coefficient [dB/(cm MHz)].

  //
  // Create an output matrix for the impulse response
  //

  plhs[0] = mxCreateDoubleMatrix(nt,No,mxREAL);
  double *h = mxGetPr(plhs[0]);

  SIRData hsir(h, nt, No);
  hsir.clear();

  //
  // Call the attenuation subroutine.
  //

  Attenuation att(nt, dt, alpha);
  FFTCVec xc(nt);
  FFTVec x(nt);

  if (mxGetM(prhs[2]) * mxGetN(prhs[2]) == 1) {
    for (dream_idx_type n=0; n<No; n++) {
      double xo = Ro[n];
      double yo = Ro[n+1*No];
      double zo = Ro[n+2*No];

      double r = std::sqrt(xo*xo + yo*yo + zo*zo);
      dream_idx_type it = (dream_idx_type) ( (r * 1.0e3/cp - delay[0])/dt + 1);
      att.att(xc, x, r, it, &h[n*nt], 1.0);
    }
  } else {
    for (dream_idx_type n=0; n<No; n++) {
      double xo = Ro[n];
      double yo = Ro[n+1*No];
      double zo = Ro[n+2*No];

      double r = std::sqrt(xo*xo + yo*yo + zo*zo);
      dream_idx_type it = (dream_idx_type) ( (r * 1.0e3/cp - delay[n])/dt + 1);
      att.att(xc, x, r, it, &h[n*nt], 1.0);
    }
  }

  return;
}
