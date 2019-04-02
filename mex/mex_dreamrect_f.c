/***
*
* Copyright (C) 2003,2006,2007,2008,2009,2014,2016 Fredrik Lingvall
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
#include "dreamrect_f.h"
#include "dream_error.h"

#ifdef USE_FFTW
#include "att.h"
#endif

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
  //printf("Caught signal SIGTERM.\n");
  running = FALSE;
}

void sig_abrt_handler(int signum) {
  //printf("Caught signal SIGABRT.\n");
}

void sig_keyint_handler(int signum) {
  //printf("Caught signal SIGINT.\n");
}



/***
 *
 * Matlab (MEX) gateway function for dreamrect_f.
 *
 ***/

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *RESTRICT ro, *RESTRICT geom_par, *RESTRICT s_par, *RESTRICT m_par;
  dream_idx_type  nt,no,n;
  int    ifoc=0;
  char   foc_met[50];
  int    buflen;
  double xo,yo,zo,a,b, dx, dy, dt;
  double *RESTRICT delay, v, cp, alpha, focal=0;
  double *RESTRICT h, *err_p;
  int    err_level=STOP, err=NONE, out_err = NONE, set = FALSE;
  char   err_str[50];
  sighandler_t old_handler, old_handler_abrt, old_handler_keyint;

  // Check for proper number of arguments

  if (!((nrhs == 7) || (nrhs == 8))) {
    dream_err_msg("dreamrect_f requires 7 or 8 input arguments!");
  }
  else
    if (nlhs > 2) {
      dream_err_msg("Too many output arguments for dreamrect_f!");
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
  // Transducer geometry
  //

  // Check that arg 2 is a 2 element vector
  if (!((mxGetM(prhs[1])==2 && mxGetN(prhs[1])==1) || (mxGetM(prhs[1])==1 && mxGetN(prhs[1])==2)))
    dream_err_msg("Argument 2 must be a 2 element vector!");

  geom_par = mxGetPr(prhs[1]);
  a = geom_par[0];		// x-size.
  b = geom_par[1];		// y-size.

  //
  // Temporal and spatial sampling parameters.
  //

  // Check that arg 3 is a 4 element vector
  if (!((mxGetM(prhs[2])==4 && mxGetN(prhs[2])==1) || (mxGetM(prhs[2])==1 && mxGetN(prhs[2])==4)))
    dream_err_msg("Argument 3 must be a vector of length 4!");

  s_par = mxGetPr(prhs[2]);
  dx    = s_par[0];		// Spatial x discretization size.
  dy    = s_par[1];		// Spatial dy iscretization size.
  dt    = s_par[2];		// Temporal discretization size (= 1/sampling freq).
  nt    = (dream_idx_type) s_par[3];	// Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //

  // Check that arg 4 is a scalar (or vector).
  if ( (mxGetM(prhs[3]) * mxGetN(prhs[3]) !=1) && ((mxGetM(prhs[3]) * mxGetN(prhs[3])) != no))
    dream_err_msg("Argument 4 must be a scalar or a vector with a length equal to the number of observation points!");

  delay = mxGetPr(prhs[3]);

  //
  // Material parameters
  //

  // Check that arg 5 is a 3 element vectora
  if (!((mxGetM(prhs[4])==3 && mxGetN(prhs[4])==1) || (mxGetM(prhs[4])==1 && mxGetN(prhs[4])==3)))
    dream_err_msg("Argument 5 must be a vector of length 3!");

  m_par = mxGetPr(prhs[4]);
  v     = m_par[0]; // Normal velocity of transducer surface.
  cp    = m_par[1]; // Sound speed.
  alpha  = m_par[2]; // Attenuation coefficient [dB/(cm MHz)].


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

  if (set == FALSE)
      dream_err_msg("Unknown focusing method!");

    // Check that arg 7 is a scalar.
    if (mxGetM(prhs[6]) * mxGetN(prhs[6]) !=1 )
      dream_err_msg("Argument 7 must be a scalar!");

    // Focal point (in mm).
    focal = mxGetScalar(prhs[6]);

  } else
    ifoc = 1;

  // Error reporting.
  if (nrhs == 8) {

    if (!mxIsChar(prhs[7]))
      dream_err_msg("Argument 8 must be a string");

    buflen = (mxGetM(prhs[7]) * mxGetN(prhs[7]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[7],err_str,buflen);

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
    printf("Couldn't register SIGTERM  signal handler.\n");
  }

  if (( old_handler_abrt=signal(SIGABRT, &sighandler)) == SIG_ERR) {
    printf("Couldn't register SIGABRT signal handler.\n");
  }

  if (( old_handler_keyint=signal(SIGINT, &sighandler)) == SIG_ERR) {
    printf("Couldn't register SIGINT signal handler.\n");
  }

  //
  // Call the DREAM subroutine.
  //

  running = TRUE;

#ifdef USE_FFTW
  if (alpha != (double) 0.0)
    att_init(nt,1);
#endif

  if (mxGetM(prhs[3]) * mxGetN(prhs[3]) == 1) {
    for (n=0; n<no; n++) {
      xo = ro[n];
      yo = ro[n+1*no];
      zo = ro[n+2*no];
      err = dreamrect_f(xo,yo,zo,a,b,dx,dy,dt,nt,delay[0],v,cp,alpha,ifoc,focal,&h[n*nt],err_level);
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
      err = dreamrect_f(xo,yo,zo,a,b,dx,dy,dt,nt,delay[n],v,cp,alpha,ifoc,focal,&h[n*nt],err_level);
      if (err != NONE)
        out_err = err;

      if (!running) {
        break; // CTRL-C pressed.
      }

    }
  }

  //
  // Restore old signal handlers.
  //

  if (signal(SIGTERM, old_handler) == SIG_ERR) {
    printf("Couldn't register old SIGTERM signal handler.\n");
  }

  if (signal(SIGABRT,  old_handler_abrt) == SIG_ERR) {
    printf("Couldn't register old SIGABRT signal handler.\n");
  }

  if (signal(SIGINT, old_handler_keyint) == SIG_ERR) {
    printf("Couldn't register old SIGINT signal handler.\n");
  }

#ifdef USE_FFTW
  if (alpha != (double) 0.0)
    att_close();
#endif

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
