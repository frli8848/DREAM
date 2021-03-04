/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2014,2016,2021 Fredrik Lingvall
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
#include "dream_error.h"

#include "mex.h"

/***
 *
 * Matlab (MEX) gateway function for dreamrect.
 *
 ***/

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *ro, *s_par, *m_par;
  int    it, nt, no, n;
  double xo, yo, zo, dt;
  double *delay, cp, alpha, r;
  double *h;

  // Check for proper number of arguments

  if (nrhs != 4) {
    dream_err_msg("dream_att requires 4 input arguments!");
  } else {
    if (nlhs > 1) {
      dream_err_msg("dream_att requires one output argument!");
    }
  }

  //
  // Observation point.
  //

  // Check that arg (number of observation points) x 3 matrix
  if (!mxGetN(prhs[0])==3) {
    dream_err_msg("Argument 1 must be a (number of observation points) x 3 matrix!");
  }

  no = mxGetM(prhs[0]); // Number of observation points.
  ro = mxGetPr(prhs[0]);

  //
  // Temporal and spatial sampling parameters.
  //

  // Check that arg 2 is a 4 element vector
  if (!((mxGetM(prhs[1])==2 && mxGetN(prhs[1])==1) || (mxGetM(prhs[1])==1 && mxGetN(prhs[1])==2))) {
    dream_err_msg("Argument 2 must be a vector of length 2!");
  }

  s_par = mxGetPr(prhs[1]);
  dt    = s_par[0];		// Temporal discretization size (= 1/sampling freq).
  nt    = (dream_idx_type) s_par[1];	// Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //

  // Check that arg 3 is a scalar (or vector).
  if ( (mxGetM(prhs[2]) * mxGetN(prhs[2]) !=1) && ((mxGetM(prhs[2]) * mxGetN(prhs[2])) != no)) {
    dream_err_msg("Argument 3 must be a scalar or a vector with a length equal to the number of observation points!");
  }

  delay = mxGetPr(prhs[2]);

  //
  // Material parameters
  //

  // Check that arg 4 is a 2 element vector.
  if (!((mxGetM(prhs[3])==2 && mxGetN(prhs[3])==1) || (mxGetM(prhs[3])==1 && mxGetN(prhs[3])==2))) {
    dream_err_msg("Argument 4 must be a vector of length 2!");
  }

  m_par = mxGetPr(prhs[3]);
  cp    = m_par[0]; // Sound speed.
  alpha  = m_par[1]; // Attenuation coefficient [dB/(cm MHz)].

  //
  // Create an output matrix for the impulse response
  //
  plhs[0] = mxCreateDoubleMatrix(nt,no,mxREAL);
  h = mxGetPr(plhs[0]);

  //
  // Call the attenuation subroutine.
  //

  Attenuation att(nt, dt, alpha);
  FFTCVec xc(nt);
  FFTVec x(nt);

  if (mxGetM(prhs[2]) * mxGetN(prhs[2]) == 1) {
    for (n=0; n<no; n++) {
      xo = ro[n];
      yo = ro[n+1*no];
      zo = ro[n+2*no];

      r = std::sqrt(xo*xo + yo*yo + zo*zo);
      it = (dream_idx_type) ( (r * 1.0e3/cp - delay[0])/dt + 1);
      att.att(xc, x, r, it, &h[n*nt], 1.0);
    }
  } else {
    for (n=0; n<no; n++) {
      xo = ro[n];
      yo = ro[n+1*no];
      zo = ro[n+2*no];

      r = std::sqrt(xo*xo + yo*yo + zo*zo);
      it = (dream_idx_type) ( (r * 1.0e3/cp - delay[n])/dt + 1);
      att.att(xc, x, r, it, &h[n*nt], 1.0);
    }
  }

  return;
}
