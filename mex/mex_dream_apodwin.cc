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
  char apod_met[40];
  int buflen;
  int apod_type=0, i;
  mwSize num_elements=0;
  bool is_set = false;
  double *apod=nullptr, weight, xs, ys, ramax, param;
  double *h;

  // Check for proper number of arguments
  if (nrhs != 3) {
    dream_err_msg("dream_apodwin requires three input arguments!");
  }
  else
    if (nlhs > 1) {
      dream_err_msg("dream_apodwin requires one output argument!");
    }


  buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0]) * sizeof(mxChar)) + 1;
  mxGetString(prhs[0],apod_met,buflen);

  is_set = false;
  if (!strcmp(apod_met,"off")) {
    do_apod = false;
    is_set = true;
  }

  if (!strcmp(apod_met,"ud")) {
    do_apod = true;
    apod_type = 0;
    dream_err_msg(" 'ud'- (user defined) meaningless for this function!");
  }

  if (!strcmp(apod_met,"triangle")) {
    do_apod = true;
    apod_type = 1;
    is_set = true;
  }

  if (!strcmp(apod_met,"gauss")) {
    do_apod = true;
    apod_type = 2;
    is_set = true;
  }

  if (!strcmp(apod_met,"raised")) {
    do_apod = true;
    apod_type = 3;
    is_set = true;
  }

  if (!strcmp(apod_met,"simply")) {
    do_apod = true;
    apod_type = 4;
    is_set = true;
  }

  if (!strcmp(apod_met,"clamped")) {
    do_apod = true;
    apod_type = 5;
    is_set = true;
  }

  if (is_set == false)
    dream_err_msg("Unknown apodization level!");


  //
  // Apodization.
  //

  // apod_type = 0 - user defined, 1 traingle, 2 Gauss, 3 raised cosine, 4 simply supported, 5 clamped.

  if (!mxIsChar(prhs[0]))
    dream_err_msg("Argument 1 must be a string");

  if (mxGetM(prhs[1]) * mxGetN(prhs[1]) !=1)
    dream_err_msg("Argument 2 must be a scalar!");

  num_elements = (int) mxGetScalar(prhs[1]);

  //if (!mxIsInt32(prhs[1]) || num_elements < 0)
  if (num_elements < 0)
    dream_err_msg("Argument 2 must be a positive integer!");

  //
  // param - Parameter used for raised cos and Gaussian apodization functions.
  //

  if (mxGetM(prhs[2]) * mxGetN(prhs[2]) !=1)
    dream_err_msg("Argument 3 must be a scalar!");


  param = mxGetScalar(prhs[2]);

  // Create a matrix for return arguments
  plhs[0] = mxCreateDoubleMatrix(num_elements,1,mxREAL);
  h = mxGetPr(plhs[0]);

  ramax = 1;
  ys = 0;
  if (do_apod) {
    for (i=0; i<num_elements; i++) {
      xs = 2*ramax * (0.5 - ((double) i / (double) num_elements));
      apodization(apod_type, i, apod, &weight, xs, ys, ramax, param);
      h[i] = weight;
    }
  }
}
