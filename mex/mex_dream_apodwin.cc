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

#include "arr_functions.h"
#include "arg_parser.h"

#include "mex.h"

/***
 *
 * Matlab (mex) gateway function for dream_apodwin (for the apodization functions in
 * arr_functions.c)
 *
 ***/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  bool do_apod=false;
  ApodMet apod_met=ApodMet::gauss;
  double *h;

  ArgParser ap;

  ap.check_arg_in("dream_apodwin", nrhs, 3, 3);
  ap.check_arg_out("dream_apodwin", nlhs, 0, 1);

  dream_idx_type num_elements = (dream_idx_type) mxGetScalar(prhs[1]);

  // Allocate space for the user defined apodization weights
  std::unique_ptr<double[]> apod = std::make_unique<double[]>(num_elements);

  double apod_par;
  ap.parse_apod_args("dream_apodwin", prhs, 0, num_elements,
                     do_apod, apod.get(), apod_met, apod_par);

  // Create an output matrix for the impulse response
  plhs[0] = mxCreateDoubleMatrix(num_elements,1,mxREAL);
  h = mxGetPr(plhs[0]);

  double weight=1.0;
  double ramax = 1.0;
  double ys = 0.0;
  if (do_apod) {
    for (dream_idx_type n=0; n<num_elements; n++) {
      double xs = 2*ramax * (0.5 - ((double) n / (double) num_elements));
      apodization(apod_met, n, apod.get(), &weight, xs, ys, ramax, apod_par);
      h[n] = weight;
    }
  }
}
