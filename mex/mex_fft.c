/***
*
* Copyright (C) 2002,2004,2006,2007,2008,2009,2014,2016 Fredrik Lingvall
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

// $Revision: 890 $ $Date: 2016-09-27 16:56:06 +0200 (Tue, 27 Sep 2016) $ $LastChangedBy: frli8848 $

#include <uchar.h>
#include "mex.h"
#include "fft.h"

/***
 * Gateway function for fft.c
 *
 * This is  a MEX file for MATLAB.
 *
 ***/

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  dream_idx_type nt;
  double *RESTRICT xr, *RESTRICT xi, *RESTRICT yr, *RESTRICT yi;

  // Check for proper number of arguments
  if (nrhs != 1) {
    mexErrMsgTxt("my_ifft requires one input arguments!");
  }
  else
    if (nlhs > 1) {
      mexErrMsgTxt("my_ifft requires one output argument!");
    }
  
  //  Get the input parameters.
  xr    = mxGetPr(prhs[0]);
  xi    = mxGetPi(prhs[0]);
  nt =  (dream_idx_type) mxGetM(prhs[0]) * mxGetN(prhs[0]);
  
  // Create a matrix for return arguments
  plhs[0] = mxCreateDoubleMatrix(nt,1,mxCOMPLEX);
  yr = mxGetPr(plhs[0]);
  yi = mxGetPi(plhs[0]);
  
  cc_ifft(xr, xi, yr, yi, nt);
}
      
