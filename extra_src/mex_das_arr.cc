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

#include <iostream>
#include <csignal>
#include <cstring>

#include "das_arr.h"
#include "dream_error.h"

#include "mex.h"

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
  running = false;
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
  char   apod_str[50],foc_str[50],steer_str[50];
  int    buflen;
  double xo,yo,zo,dt;
  size_t nt,no,n;
  double param=0,*delay,cp;
  int    isize=0;
  double *gx, *gy, *gz;
  FocusMet foc_met=FocusMet::none;
  double *focal=nullptr;
  SteerMet steer_met=SteerMet::none;
  double theta=0,phi=0,*apod=NULL;
  bool    do_apod=false;
  ApodMet apod_met=ApodMet::gauss;
  double *h, *err_p;
  ErrorLevel err_level=ErrorLevel::stop, err=ErrorLevel::none, out_err = ErrorLevel::none;
  bool is_set = false;
  char err_str[50];
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

  DelayType delay_type=DelayType::single;

  if ( (mxGetM(prhs[3]) * mxGetN(prhs[3]) !=1) && ((mxGetM(prhs[3]) * mxGetN(prhs[3])) != no)) {
    dream_err_msg("Argument 4 must be a scalar or a vector with a length equal to the number of observation points!");
  } else {
    if (mxGetM(prhs[3]) * mxGetN(prhs[3]) ==1) {
      delay_type=DelayType::single;
    } else {
      delay_type=DelayType::multiple;
    }
  }

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

  if (nrhs >= 6) {

   if (!mxIsChar(prhs[5]))
      dream_err_msg("Argument 6 must be a string");

    buflen = (mxGetM(prhs[5]) * mxGetN(prhs[5]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[5],foc_str,buflen);

    is_set = false;

    if (!strcmp(foc_str,"off")) {
      foc_met = FocusMet::none;
      is_set = true;
    }

    if (!strcmp(foc_str,"x")) {
      foc_met = FocusMet::x;
      is_set = true;
    }

    if (!strcmp(foc_str,"y")) {
      foc_met = FocusMet::y;
      is_set = true;
    }

    if (!strcmp(foc_str,"xy")) {
      foc_met = FocusMet::xy;
      is_set = true;
    }

    if (!strcmp(foc_str,"x+y")) {
      foc_met = FocusMet::x_y;
      is_set = true;
    }

    if (!strcmp(foc_str,"ud")) {
      foc_met = FocusMet::ud;
      is_set = true;

      if (mxGetM(prhs[6]) * mxGetN(prhs[6]) != isize ) {
        dream_err_msg("The time delay vector (argument 7) for user defined ('ud') focusing\n delays must have the same length as the number of array elements.!");
      }
      focal = mxGetPr(prhs[6]);
    } else {

      // Check that arg 7 is a scalar.
      if (mxGetM(prhs[6]) * mxGetN(prhs[6]) !=1 )
        dream_err_msg("Argument 7 to must be a scalar!");

      // Focal point (in mm).
      focal = mxGetPr(prhs[6]);
    }

    if (is_set == false) {
      dream_err_msg("Unknown focusing method!");
    }

  } else {
    foc_met = FocusMet::none;
  }

  //
  // Beam steering.
  //

  if (nrhs >= 8) {

   if (!mxIsChar(prhs[7]))
      dream_err_msg("Argument 8 must be a string");

    buflen = (mxGetM(prhs[7]) * mxGetN(prhs[7]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[7],steer_str,buflen);

    steer_met = SteerMet::none; // Default no steering
    is_set = false;

    if (!strcmp(steer_str,"off")) {
      steer_met = SteerMet::none;
      is_set = true;
    }

    if (!strcmp(steer_str,"x")) {
      steer_met = SteerMet::x;
      is_set = true;
    }

    if (!strcmp(steer_str,"y")) {
      steer_met = SteerMet::y;
      is_set = true;
    }

    if (!strcmp(steer_str,"xy")) {
      steer_met = SteerMet::xy;
      is_set = true;
    }

    if (is_set == false) {
      dream_err_msg("Unknown beamsteering method!");
    }

    // Check that arg 9 is a 2 element vector
    if (!((mxGetM(prhs[8])==2 && mxGetN(prhs[8])==1) || (mxGetM(prhs[8])==1 && mxGetN(prhs[8])==2))) {
      dream_err_msg("Argument 9 must be a vector of length 2!");
    }

    steer_par = mxGetPr(prhs[8]);
    theta  = steer_par[0];		// Angle in x-direction.
    phi    = steer_par[1];		// Angle in y-direction.

  } else {
    steer_met = SteerMet::none;
  }

  //
  // Apodization.
  //

  if (nrhs >= 10) {

    if (!mxIsChar(prhs[9]))
      dream_err_msg("Argument 10 must be a string");


    buflen = (mxGetM(prhs[9]) * mxGetN(prhs[9]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[9],apod_str,buflen);

    do_apod = false;			// default off.
    is_set = false;

    if (!strcmp(apod_str,"off")) {
      do_apod = false;
      is_set = true;
    }

    if (!strcmp(apod_str,"ud")) {
      do_apod = true;
      apod_met = ApodMet::ud;
      is_set = true;

      // Vector of apodization weights.
      if (mxGetM(prhs[10]) * mxGetN(prhs[10]) != isize)
        dream_err_msg("The length of argument 11 (apodization vector) must be the same as the number of array elements!");

      apod = mxGetPr(prhs[10]);
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
      apod_met = ApodMet::simply_supported;
      is_set = true;
    }

    if (is_set == false) {
      dream_err_msg("Unknown apodization method!");
    }

    // Parameter for raised cos and Gaussian apodization functions.
    if (mxGetM(prhs[11]) * mxGetN(prhs[11]) !=1 ) {
      dream_err_msg("Argument 12 must be a scalar");
    }

    param = mxGetScalar(prhs[11]);
  } else {
    do_apod = false;
  }

  //
  // Error reporting.
  //

  if (nrhs == 13) {

   if (!mxIsChar(prhs[12]))
      dream_err_msg("Argument 13 must be a string");

    buflen = (mxGetM(prhs[12]) * mxGetN(prhs[12]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[12],err_str,buflen);

    is_set = false;

    if (!strcmp(err_str,"ignore")) {
      err_level = ErrorLevel::ignore;
      is_set = true;
    }

    if (!strcmp(err_str,"warn")) {
      err_level = ErrorLevel::warn;
      is_set = true;
    }

    if (!strcmp(err_str,"stop")) {
      err_level = ErrorLevel::stop;
      is_set = true;
    }

    if (is_set == false) {
      dream_err_msg("Unknown error level!");
    }
  } else {
    err_level = ErrorLevel::stop; // Default.
  }

  // Create an output matrix for the impulse response
  plhs[0] = mxCreateDoubleMatrix(nt,no,mxREAL);
  h = mxGetPr(plhs[0]);

  //
  // Register signal handlers.
  //

  if ((old_handler = std::signal(SIGTERM, &sighandler)) == SIG_ERR) {
    std::cerr << "Couldn't register SIGTERM signal handler!" << std::endl;
  }

  if ((old_handler_abrt = std::signal(SIGABRT, &sighandler)) == SIG_ERR) {
    std::cerr << "Couldn't register SIGABRT signal handler!" << std::endl;
  }

  if ((old_handler_keyint = std::signal(SIGINT, &sighandler)) == SIG_ERR) {
    std::cerr << "Couldn't register SIGINT signal handler!" << std::endl;
  }

  //
  // Call the DAS subroutine.
  //

  running = true;

  for (n=0; n<no; n++) {
    xo = ro[n];
    yo = ro[n+1*no];
    zo = ro[n+2*no];

    double dlay = 0.0;
    if (delay_type == DelayType::single) {
      dlay = delay[0];
    } else { // DelayType::multiple.
      dlay = delay[n];
    }

    err = das_arr(xo, yo, zo,
                  dt, nt,
                  dlay,
                  cp,
                  isize,
                  gx, gy, gz,
                  foc_met, focal,
                  steer_met, theta, phi,
                  apod, do_apod, apod_met, param,
                  &h[n*nt],err_level);

    if (err != ErrorLevel::none) {
      out_err = err;
    }
  }

  //
  // Restore old signal handlers.
  //

  if (std::signal(SIGTERM, old_handler) == SIG_ERR) {
    std::cerr << "Couldn't register old SIGTERM signal handler!" << std::endl;
  }

  if (std::signal(SIGABRT, old_handler_abrt) == SIG_ERR) {
    std::cerr << "Couldn't register old SIGABRT signal handler!" << std::endl;
  }

  if (std::signal(SIGINT, old_handler_keyint) == SIG_ERR) {
    std::cerr << "Couldn't register old SIGINT signal handler!" << std::endl;
  }

  if (!running) {
    dream_err_msg("CTRL-C pressed!\n"); // Bail out.
  }

  // Return error.
  if (nlhs == 2) {
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    err_p =  mxGetPr(plhs[1]);
    err_p[0] = (double) out_err;
  }

  return;
}
