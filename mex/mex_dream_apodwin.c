/***
*
* Copyright (C) 2002,2003,2006,2008,2009,2014,2016 Fredrik Lingvall
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
#include "mex.h"
#include "arr_functions.h"
#include "dream_error.h"

/***
 *
 * Matlab (mex) gateway function for dream_apodwin (for the apodization functions in
 * arr_functions.c)
 *
 ***/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char   apod_met[40];
  int    buflen;
  mwSize iweight=0,iapo=0,i,isize=0;
  int set = FALSE;
  double *RESTRICT apod=NULL, weight, xs, ys, ramax, param;
  double *RESTRICT yr;


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

  do_apod = false;			// default off.

  set = FALSE;
  if (!strcmp(apod_met,"off")) {
    do_apod = false;
    set = TRUE;
  }

  if (!strcmp(apod_met,"ud")) {
    do_apod = true;
    iapo = 0;
    dream_err_msg(" 'ud'- (user defined) meaningless for this function!");
  }

  if (!strcmp(apod_met,"triangle")) {
    do_apod = true;
    iapo = 1;
    set = TRUE;
  }

  if (!strcmp(apod_met,"gauss")) {
    do_apod = true;
    iapo = 2;
    set = TRUE;
  }

  if (!strcmp(apod_met,"raised")) {
    do_apod = true;
    iapo = 3;
    set = TRUE;
  }

  if (!strcmp(apod_met,"simply")) {
    do_apod = true;
    iapo = 4;
    set = TRUE;
  }

  if (!strcmp(apod_met,"clamped")) {
    do_apod = true;
    iapo = 5;
    set = TRUE;
  }

  if (set == FALSE)
    dream_err_msg("Unknown apodization level!");


  //
  // Apodization.
  //

  // iapo = 0 - user defined, 1 traingle, 2 Gauss, 3 raised cosine, 4 simply supported, 5 clamped.

  if (!mxIsChar(prhs[0]))
    dream_err_msg("Argument 1 must be a string");

  if (mxGetM(prhs[1]) * mxGetN(prhs[1]) !=1)
    dream_err_msg("Argument 2 must be a scalar!");

  isize = (int) mxGetScalar(prhs[1]);

  //if (!mxIsInt32(prhs[1]) || isize < 0)
  if (isize < 0)
    dream_err_msg("Argument 2 must be a positive integer!");

  //
  // param - Parameter used for raised cos and Gaussian apodization functions.
  //

  if (mxGetM(prhs[2]) * mxGetN(prhs[2]) !=1)
    dream_err_msg("Argument 3 must be a scalar!");


  param = mxGetScalar(prhs[2]);

  // Create a matrix for return arguments
  plhs[0] = mxCreateDoubleMatrix(isize,1,mxREAL);
  yr = mxGetPr(plhs[0]);
  //yi = mxGetPi(plhs[0]);

  ramax = 1;

  ys = 0;
  if (iweight != 1) {
    for (i=0; i<isize; i++) {
      xs = 2*ramax * (0.5 - ((double) i / (double) isize));
      apodization(iweight,iapo,i,apod,&weight,xs,ys,ramax,param,isize);
      yr[i] = weight;
    }
  }
}
