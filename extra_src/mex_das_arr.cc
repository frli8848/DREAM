/***
*
* Copyright (C) 2003,2006,2007,2008,2009,2014,2015,2016,2021 Fredrik Lingvall
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
#include <signal.h>
#include <uchar.h>
#include "mex.h"
#include "das_arr.h"
#include "dream_error.h"

//
// Globals
//

int running;

//
// Function prototypes.
//

void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

//
// typedef:s
//

typedef void (*sighandler_t)(int);

/***
 *
 * Signal handlers.
 *
 ***/

void sighandler(int signum) {
  //mexPrintf("Caught signal SIGTERM.\n");
  running = FALSE;
}

void sig_abrt_handler(int signum) {
  //mexPrintf("Caught signal SIGABRT.\n");
}

void sig_keyint_handler(int signum) {
  //mexPrintf("Caught signal SIGINT.\n");
}


/***
 *
 * Gateway function for das_arr (delay-and-sum for arrays).
 *
 * This is  a MEX file for MATLAB.
 *
 ***/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *ro,*s_par,*m_par;
  double *steer_par;
  char   apod_met[50],foc_met[50],steer_met[50];
  int    buflen;
  double xo,yo,zo,dt;
  size_t nt,no,n;
  double param=0,*delay,cp;
  int    isize=0;
  double *gx, *gy, *gz;
  int    ifoc=0;
  double focal=0, *ud_focal=NULL;
  int    ister=0;
  double theta=0,phi=0,*apod=NULL;
  int    iweight=0, iapo=0;
  double *h, *err_p;
  int    err_level=STOP, err=NONE, out_err = NONE, set = FALSE;
  char   err_str[50];
  sighandler_t old_handler, old_handler_abrt, old_handler_keyint;

  // Check for proper number of arguments

  if (!((nrhs == 12) || (nrhs == 13))) {
    dream_err_msg("das_arr requires 12 or 13 input arguments!");
  }
  else
    if (nlhs > 2) {
      dream_err_msg("Too many output arguments for das_arr !");
    }

  //
  // Observation point.
  //

  // Check that arg (number of observation points) x 3 matrix
  if (!mxGetN(prhs[0])==3)
    dream_err_msg("Argument 1 must be a (number of observation points) x 3 matrix!");

  no = mxGetM(prhs[0]); // Number of observation points.
  ro = mxGetPr(prhs[0]);

  //
  // Grid function (position vectors of the elements).
  //

  isize = (int) mxGetM(prhs[1]); // Number of elementents in the array.
  if (mxGetN(prhs[1]) !=3 )
    dream_err_msg("Argument 2 must a (number of array elements) x 3 matrix!");

  gx    = mxGetPr(prhs[1]);	// First column in the matrix.
  gy    = gx + isize;		// Second column in the matrix.
  gz    = gy + isize;		// Third column in the matrix.

  //
  // Temporal and spatial sampling parameters.
  //

  // Check that arg 3 is a 2 element vector
  if (!((mxGetM(prhs[2])==2 && mxGetN(prhs[2])==1) || (mxGetM(prhs[2])==1 && mxGetN(prhs[2])==2)))
    dream_err_msg("Argument 3 must be a vector of length 2!");

  s_par = mxGetPr(prhs[2]);
  dt    = s_par[0];		// Temporal discretization size (= 1/sampling freq).
  nt    = (size_t) s_par[1];	// Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //

  // Check that arg 4 is a scalar.
  if ( (mxGetM(prhs[3]) * mxGetN(prhs[3]) !=1) && ((mxGetM(prhs[3]) * mxGetN(prhs[3])) != no))
    dream_err_msg("Argument 4 must be a scalar or a vector with a length equal to the number of observation points!");

  delay = mxGetPr(prhs[3]);

  //
  // Material parameters
  //

  // Check that arg 5 is a scalar
  if (!(mxGetM(prhs[4])==1 && mxGetN(prhs[4])==1))
    dream_err_msg("Argument 5 must be a scalar!");

  m_par = mxGetPr(prhs[4]);
  cp    = m_par[0]; // Sound speed.

  //
  // Focusing parameters.
  //

  //  ifoc = 1 - no foc, 2 foc x ,3 foc y, 4 foc xy (del=fsqrt(x*x+y*y)), 5 focx+focy.

  if (nrhs >= 6) {

   if (!mxIsChar(prhs[5]))
      dream_err_msg("Argument 6 must be a string");

    buflen = (mxGetM(prhs[5]) * mxGetN(prhs[5]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[5],foc_met,buflen);

    set = FALSE;

    if (!strcmp(foc_met,"off")) {
      ifoc = 1;
      set = TRUE;
    }

    if (!strcmp(foc_met,"x")) {
      ifoc = 2;
      set = TRUE;
    }

    if (!strcmp(foc_met,"y")) {
      ifoc = 3;
      set = TRUE;
    }

    if (!strcmp(foc_met,"xy")) {
      ifoc = 4;
      set = TRUE;
    }

    if (!strcmp(foc_met,"x+y")) {
      ifoc = 5;
      set = TRUE;
    }

    if (!strcmp(foc_met,"ud")) {
      ifoc = 6;
      set = TRUE;

      if (mxGetM(prhs[6]) * mxGetN(prhs[6]) != isize ) {
        dream_err_msg("The time delay vector (argument 7) for user defined ('ud') focusing\n delays must have the same length as the number of array elements.!");
      }
      ud_focal = mxGetPr(prhs[6]);
    }
    else {

      // Check that arg 7 is a scalar.
      if (mxGetM(prhs[6]) * mxGetN(prhs[6]) !=1 )
        dream_err_msg("Argument 7 to must be a scalar!");

      // Focal point (in mm).
      focal = mxGetScalar(prhs[6]);
    }

    if (set == FALSE)
      dream_err_msg("Unknown focusing method!");

  } else
    ifoc = 1;


  //
  // Beam steering.
  //

  // Beam steering: ister = 1 - no steering, 2 steer ph=ax ,3 steer y ph=by, 4 steer xy ph=ax+by.

  if (nrhs >= 8) {

   if (!mxIsChar(prhs[7]))
      dream_err_msg("Argument 8 must be a string");

    buflen = (mxGetM(prhs[7]) * mxGetN(prhs[7]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[7],steer_met,buflen);

    ister = 1;			// Default no steering
    set = FALSE;

    if (!strcmp(steer_met,"off")) {
      ister = 1;
      set = TRUE;
    }

    if (!strcmp(steer_met,"x")) {
      ister = 2;
      set = TRUE;
    }

    if (!strcmp(steer_met,"y")) {
      ister = 3;
      set = TRUE;
    }

    if (!strcmp(steer_met,"xy")) {
      ister = 4;
      set = TRUE;
    }

    if (set == FALSE)
      dream_err_msg("Unknown beamsteering method!");

    // Check that arg 9 is a 2 element vector
    if (!((mxGetM(prhs[8])==2 && mxGetN(prhs[8])==1) || (mxGetM(prhs[8])==1 && mxGetN(prhs[8])==2)))
      dream_err_msg("Argument 9 must be a vector of length 2!");

    steer_par = mxGetPr(prhs[8]);
    theta  = steer_par[0];		// Angle in x-direction.
    phi    = steer_par[1];		// Angle in y-direction.

  } else
    ister = 1;

  //
  // Apodization.
  //

  // iweight = 1 - no apodization, 2  apodization.
  // iapo = 0 - user defined, 1 traingle, 2 Gauss, 3 raised cosine, 4 simply supported, 5 clamped.

  if (nrhs >= 10) {

    if (!mxIsChar(prhs[9]))
      dream_err_msg("Argument 10 must be a string");


    buflen = (mxGetM(prhs[9]) * mxGetN(prhs[9]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[9],apod_met,buflen);

    iweight = 1;			// default off.
    set = FALSE;

    if (!strcmp(apod_met,"off")) {
      iweight = 1;
      set = TRUE;
    }

    if (!strcmp(apod_met,"ud")) {
      iweight = 2;
      iapo = 0;
      set = TRUE;

      // Vector of apodization weights.
      if (mxGetM(prhs[10]) * mxGetN(prhs[10]) != isize)
        dream_err_msg("The length of argument 11 (apodization vector) must be the same as the number of array elements!");

      apod = mxGetPr(prhs[10]);
    }

    if (!strcmp(apod_met,"triangle")) {
      iweight = 2;
      iapo = 1;
      set = TRUE;
    }

    if (!strcmp(apod_met,"gauss")) {
      iweight = 2;
      iapo = 2;
      set = TRUE;
    }

    if (!strcmp(apod_met,"raised")) {
      iweight = 2;
      iapo = 3;
      set = TRUE;
    }

    if (!strcmp(apod_met,"simply")) {
      iweight = 2;
      iapo = 4;
      set = TRUE;
    }

    if (!strcmp(apod_met,"clamped")) {
      iweight = 2;
      iapo = 5;
      set = TRUE;
    }

   if (set == FALSE)
      dream_err_msg("Unknown apodization method!");

    // Parameter for raised cos and Gaussian apodization functions.
    if (mxGetM(prhs[11]) * mxGetN(prhs[11]) !=1 )
      dream_err_msg("Argument 12 must be a scalar");

    param = mxGetScalar(prhs[11]);
  }
  else
    iweight = 1;

  //
  // Error reporting.
  //
  if (nrhs == 13) {

   if (!mxIsChar(prhs[12]))
      dream_err_msg("Argument 13 must be a string");

    buflen = (mxGetM(prhs[12]) * mxGetN(prhs[12]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[12],err_str,buflen);

    set = FALSE;

    if (!strcmp(err_str,"ignore")) {
      err_level = IGNORE;
      set = TRUE;
    }

    if (!strcmp(err_str,"warn")) {
      err_level = WARN;
      set = TRUE;
    }

    if (!strcmp(err_str,"stop")) {
      err_level = STOP;
      set = TRUE;
    }

    if (set == FALSE)
      dream_err_msg("Unknown error level!");
  }
  else
    err_level = STOP; // Default.

  // Create an output matrix for the impulse response
  plhs[0] = mxCreateDoubleMatrix(nt,no,mxREAL);
  h = mxGetPr(plhs[0]);

  //
  // Register signal handlers.
  //

  if ((old_handler = signal(SIGTERM, &sighandler)) == SIG_ERR) {
    mexPrintf("Couldn't register SIGTERM  signal handler.\n");
  }

  if (( old_handler_abrt=signal(SIGABRT, &sighandler)) == SIG_ERR) {
    mexPrintf("Couldn't register SIGABRT signal handler.\n");
  }

  if (( old_handler_keyint=signal(SIGINT, &sighandler)) == SIG_ERR) {
    mexPrintf("Couldn't register SIGINT signal handler.\n");
  }

  //
  // Call the DAS subroutine.
  //

  running = TRUE;

  if (ifoc != 6) {

    if (mxGetM(prhs[3]) * mxGetN(prhs[3]) == 1) {
      for (n=0; n<no; n++) {
        xo = ro[n];
        yo = ro[n+1*no];
        zo = ro[n+2*no];

        err = das_arr(xo,yo,zo,dt,nt,delay[0],cp,isize,gx,gy,gz,
                             ifoc,focal,ister,theta,phi,apod,iweight,iapo,param,&h[n*nt],err_level);
        if (err != NONE)
          out_err = err;
      }
    } else {
      for (n=0; n<no; n++) {
        xo = ro[n];
        yo = ro[n+1*no];
        zo = ro[n+2*no];

        err = das_arr(xo,yo,zo,dt,nt,delay[n],cp,isize,gx,gy,gz,
                             ifoc,focal,ister,theta,phi,apod,iweight,iapo,param,&h[n*nt],err_level);
        if (err != NONE)
          out_err = err;
      }
    }
  }
  else { // User defined focusing.

    if (mxGetM(prhs[3]) * mxGetN(prhs[3]) == 1) {
      for (n=0; n<no; n++) {
        xo = ro[n];
        yo = ro[n+1*no];
        zo = ro[n+2*no];

        err = das_arr_ud(xo,yo,zo,dt,nt,delay[0],cp,isize,gx,gy,gz,
                         ifoc,ud_focal,ister,theta,phi,apod,iweight,iapo,param,&h[n*nt],err_level);
        if (err != NONE)
          out_err = err;

        if (!running) {
          break; // CTRL-C pressed.
        }

      }
    } else {
      for (n=0; n<no; n++) {
        xo = ro[n];
        yo = ro[n+1*no];
        zo = ro[n+2*no];

        err = das_arr_ud(xo,yo,zo,dt,nt,delay[n],cp,isize,gx,gy,gz,
                         ifoc,ud_focal,ister,theta,phi,apod,iweight,iapo,param,&h[n*nt],err_level);
        if (err != NONE)
          out_err = err;

        if (!running) {
          break; // CTRL-C pressed.
        }

      }
    }
  }

  //
  // Restore old signal handlers.
  //

  if (signal(SIGTERM, old_handler) == SIG_ERR) {
    mexPrintf("Couldn't register old SIGTERM signal handler.\n");
  }

  if (signal(SIGABRT,  old_handler_abrt) == SIG_ERR) {
    mexPrintf("Couldn't register old SIGABRT signal handler.\n");
  }

  if (signal(SIGINT, old_handler_keyint) == SIG_ERR) {
    mexPrintf("Couldn't register old SIGINT signal handler.\n");
  }

  if (!running) {
    mexErrMsgTxt("CTRL-C pressed!\n"); // Bail out.
  }

  // Return error.
  if (nlhs == 2) {
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    err_p =  mxGetPr(plhs[1]);
    err_p[0] = (double) out_err;
  }

  return;
}