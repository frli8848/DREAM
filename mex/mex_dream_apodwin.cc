/***
*
* Copyright (C) 2002,2003,2006,2008,2009,2014,2016,2021 Fredrik Lingvall
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

#include <string.h>
#include <uchar.h>

#include "arr_functions.h"
#include "dream_error.h"

#include "mex.h"

/***
 *
 * Matlab (mex) gateway function for dream_apodwin (for the apodization functions in
 * arr_functions.c)
 *
 ***/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  bool do_apod = false;			// default off.
  char apod_str[40];
  int buflen;
  ApodMet apod_met=ApodMet::gauss;
  dream_idx_type i, num_elements=0;
  bool is_set = false;
  double *apod=nullptr, weight, xs, ys, ramax, apod_par;
  double *h;

  // Check for proper number of arguments
  if (nrhs != 3) {
    dream_err_msg("dream_apodwin requires three input arguments!");
  } else {
    if (nlhs > 1) {
      dream_err_msg("dream_apodwin requires one output argument!");
    }
  }

  buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0]) * sizeof(mxChar)) + 1;
  mxGetString(prhs[0],apod_str,buflen);

  is_set = false;
  if (!strcmp(apod_str,"off")) {
    do_apod = false;
    is_set = true;
  }

  if (!strcmp(apod_str,"ud")) {
    do_apod = true;
    apod_met = ApodMet::ud;
    dream_err_msg(" 'ud'- (user defined) meaningless for this function!");
  }

  if (!strcmp(apod_str,"triangle")) {
    do_apod = true;
    apod_met = ApodMet::triangle;
    is_set = true;
  }

  if (!strcmp(apod_str,"gauss")) {
    do_apod = true;
    apod_met = ApodMet::gauss;
    is_set = true;
  }

  if (!strcmp(apod_str,"raised")) {
    do_apod = true;
    apod_met = ApodMet::raised_cosine;
    is_set = true;
  }

  if (!strcmp(apod_str,"simply")) {
    do_apod = true;
    apod_met = ApodMet::simply_supported;
    is_set = true;
  }

  if (!strcmp(apod_str,"clamped")) {
    do_apod = true;
    apod_met = ApodMet::clamped;
    is_set = true;
  }

  if (is_set == false) {
    dream_err_msg("Unknown apodization level!");
  }

  //
  // Apodization.
  //

  if (!mxIsChar(prhs[0])) {
    dream_err_msg("Argument 1 must be a string");
  }

  if (mxGetM(prhs[1]) * mxGetN(prhs[1]) !=1){
    dream_err_msg("Argument 2 must be a scalar!");
  }

  num_elements = (int) mxGetScalar(prhs[1]);

  //if (!mxIsInt32(prhs[1]) || num_elements < 0)
  if (num_elements < 0)
    dream_err_msg("Argument 2 must be a positive integer!");

  //
  // apod_par - Parameter used for raised cos and Gaussian apodization functions.
  //

  if (mxGetM(prhs[2]) * mxGetN(prhs[2]) !=1)
    dream_err_msg("Argument 3 must be a scalar!");


  apod_par = mxGetScalar(prhs[2]);

  // Create a matrix for return arguments
  plhs[0] = mxCreateDoubleMatrix(num_elements,1,mxREAL);
  h = mxGetPr(plhs[0]);

  ramax = 1;
  ys = 0;
  if (do_apod) {
    for (i=0; i<num_elements; i++) {
      xs = 2*ramax * (0.5 - ((double) i / (double) num_elements));
      apodization(apod_met, i, apod, &weight, xs, ys, ramax, apod_par);
      h[i] = weight;
    }
  }
}
