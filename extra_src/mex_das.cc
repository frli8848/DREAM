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

#include "das.h"
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
 * Matlab (MEX) gateway function for das (delay-and-sum).
 *
 ***/

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *ro, *s_par, *m_par;
  size_t nt, no, n;
  double xo, yo, zo, dt;
  double *delay, cp;
  double *h, *err_p;
  ErrorLevel err_level=ErrorLevel::stop, err=ErrorLevel::none, out_err = ErrorLevel::none;
  bool is_set = false;
  char err_str[50];
  int  buflen;
  sighandler_t old_handler, old_handler_abrt, old_handler_keyint;

  // Check for proper number of arguments

  if (!((nrhs == 5) || (nrhs == 6))) {
    dream_err_msg("das requires 5 or 6 input arguments!");
  }
  else
    if (nlhs > 2) {
      dream_err_msg("Too many output arguments for das!");
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
  // Temporal and spatial sampling parameters.
  //

  // Check that arg 2 is a 2 element vector
  if (!((mxGetM(prhs[1])==2 && mxGetN(prhs[1])==1) || (mxGetM(prhs[1])==1 && mxGetN(prhs[1])==2)))
    dream_err_msg("Argument 2 must be a vector of length 2!");

  s_par = mxGetPr(prhs[1]);
  dt    = s_par[0]; // Temporal discretization size (= 1/sampling freq).
  nt    = (size_t) s_par[1];	// Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //

  // Check that arg 3 is a scalar (or vector).
  if ( (mxGetM(prhs[2]) * mxGetN(prhs[2]) !=1) && ((mxGetM(prhs[2]) * mxGetN(prhs[2])) != no))
    dream_err_msg("Argument 3 must be a scalar or a vector with a length equal to the number of observation points!");

  delay = mxGetPr(prhs[2]);

  //
  // Material parameters
  //

  // Check that arg 4 is a scalar.
  if (!(mxGetM(prhs[3])==1 && mxGetN(prhs[3])==1))
    dream_err_msg("Argument 5 must be a scalar!");

  m_par = mxGetPr(prhs[3]);
  cp    = m_par[0]; // Sound speed.

  // Error reporting.
  if (nrhs == 5) {

    if (!mxIsChar(prhs[4]))
      dream_err_msg("Argument 5 must be a string");

    buflen = (mxGetM(prhs[4]) * mxGetN(prhs[4]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[4],err_str,buflen);

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

    if (is_set == false)
      dream_err_msg("Unknown error level!");

  }
  else
    err_level = ErrorLevel::stop; // Default.


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

  if (mxGetM(prhs[2]) * mxGetN(prhs[2]) == 1) {
    for (n=0; n<no; n++) {
      xo = ro[n];
      yo = ro[n+1*no];
      zo = ro[n+2*no];

      err = das(xo,yo,zo,dt,nt,delay[0],cp,&h[n*nt],err_level);

      if (err != ErrorLevel::none)
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

      err = das(xo,yo,zo,dt,nt,delay[n],cp,&h[n*nt],err_level);

      if (err != ErrorLevel::none)
        out_err = err;

      if (!running) {
        break; // CTRL-C pressed.
      }

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
