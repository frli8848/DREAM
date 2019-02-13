/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2012,2014,2016 Fredrik Lingvall
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

#include <string.h>
#include <signal.h>
#include <uchar.h>
#include "mex.h"
#include "dreamline.h"
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
 * Matlab (MEX) gateway function for dreamline.
 *
 ***/

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *RESTRICT ro, *RESTRICT geom_par, *RESTRICT s_par, *RESTRICT m_par;
  dream_idx_type nt,no,n;
  double xo,yo,zo,a,dx,dy,dt;
  double *RESTRICT delay,v,cp,alfa;
  double *RESTRICT h=NULL, *err_p;
  float  *RESTRICT sh=NULL;			// Single precision.
  int    err_level=STOP, err=NONE, out_err = NONE, set = FALSE;
  char   err_str[50];
  int    buflen;
  sighandler_t old_handler, old_handler_abrt, old_handler_keyint;
  int data_format = DOUBLE_PRECISION;

  // Check for proper number of arguments

  if (!((nrhs == 5) || (nrhs == 6))) {
    dream_err_msg("dreamline requires 5 or 6 input arguments!");
  }
  else
    if (nlhs > 2) {
      dream_err_msg("Too many output arguments for dreamline!");
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

  // Check that arg 2 is a scalar.
  if (!(mxGetM(prhs[1])==1 && mxGetN(prhs[1])==1))
    dream_err_msg("Argument 2 must be a scalar!");

  geom_par = mxGetPr(prhs[1]);
  a = geom_par[0];		// Width of strip.

  //
  // Temporal and spatial sampling parameters.
  //

  // Check that arg 3 is a 4 or 5 element vector
  if (!((mxGetM(prhs[2])==4 && mxGetN(prhs[2])==1) || (mxGetM(prhs[2])==1 && mxGetN(prhs[2])==4) ||
        (mxGetM(prhs[2])==5 && mxGetN(prhs[2])==1) || (mxGetM(prhs[2])==1 && mxGetN(prhs[2])==5)) )
    dream_err_msg("Argument 3 must be a 4 or 5 element vector!");

  s_par = mxGetPr(prhs[2]);
  dx    = s_par[0];		// Spatial discretization size (x-direction).
  dy    = s_par[1];		// Spatial discretization size (y-direction = width of strip).
  dt    = s_par[2];		// Temporal discretization size (= 1/sampling freq).
  nt    = (dream_idx_type) s_par[3];	// Length of SIR.

  if ((mxGetM(prhs[2])==5 && mxGetN(prhs[2])==1) || (mxGetM(prhs[2])==1 && mxGetN(prhs[2])==5)) {
    data_format = (int) s_par[4];	// Single or double precision.

    if ( (data_format != DOUBLE_PRECISION) && (data_format != SINGLE_PRECISION) )
    dream_err_msg("Unknown data format detected in 5 element in the 3rd arguement!");
  }

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

  // Check that arg 5 is a 3 element vector.
  if (!((mxGetM(prhs[4])==3 && mxGetN(prhs[4])==1) || (mxGetM(prhs[4])==1 && mxGetN(prhs[4])==3)))
    dream_err_msg("Argument 5 must be a vector of length 3!");

  m_par = mxGetPr(prhs[4]);
  v     = m_par[0]; // Normal velocity of transducer surface.
  cp    = m_par[1]; // Sound speed.
  alfa  = m_par[2]; // Attenuation coefficient [dB/(cm MHz)].

  // Error reporting.
  if (nrhs == 6) {

    if (!mxIsChar(prhs[5]))
      dream_err_msg("Argument 6 must be a string");

    buflen = (mxGetM(prhs[5]) * mxGetN(prhs[5]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[5],err_str,buflen);

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

  if (data_format == DOUBLE_PRECISION) {
    // Create a double precision output matrix for the impulse response(s).
    plhs[0] = mxCreateDoubleMatrix(nt,no,mxREAL);
    h = mxGetPr(plhs[0]);
  }

  if (data_format == SINGLE_PRECISION) {
    // Create a single precision output matrix for the impulse response(s).
    plhs[0] = mxCreateNumericMatrix(nt,no,mxSINGLE_CLASS,0);
    sh = mxGetPr(plhs[0]);
  }

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

  running=TRUE;

#ifdef USE_FFTW
  if (alfa != (double) 0.0)
    att_init(nt,1);
#endif

  if (mxGetM(prhs[3]) * mxGetN(prhs[3]) == 1) {

    if (data_format == DOUBLE_PRECISION) {
      for (n=0; n<no; n++) {
        xo = ro[n];
        yo = ro[n+1*no];
        zo = ro[n+2*no];

        err = dreamline(xo,yo,zo,a,dx,dy,dt,nt,delay[0],v,cp,alfa,&h[n*nt],err_level);
        if (err != NONE)
          out_err = err;

        if (!running) {
          break; // CTRL-C pressed.
        }
      }
    }

    if (data_format == SINGLE_PRECISION) {
      for (n=0; n<no; n++) {
        xo = ro[n];
        yo = ro[n+1*no];
        zo = ro[n+2*no];

        //err = sdreamline(xo,yo,zo,a,dx,dy,dt,nt,delay[0],v,cp,alfa,&sh[n*nt],err_level);
        if (err != NONE)
          out_err = err;

        if (!running) {
          break; // CTRL-C pressed.
        }
      }
    }

  } else {

    if (data_format == DOUBLE_PRECISION) {

      for (n=0; n<no; n++) {
        xo = ro[n];
        yo = ro[n+1*no];
        zo = ro[n+2*no];

        //err = dreamline(xo,yo,zo,a,dx,dy,dt,nt,delay[n],v,cp,alfa,&h[n*nt],err_level);
        if (err != NONE)
          out_err = err;

        if (!running) {
          break; // CTRL-C pressed.
        }
      }
    }

    if (data_format == SINGLE_PRECISION) {

      for (n=0; n<no; n++) {
        xo = ro[n];
        yo = ro[n+1*no];
        zo = ro[n+2*no];

        //err = sdreamline(xo,yo,zo,a,dx,dy,dt,nt,delay[n],v,cp,alfa,&sh[n*nt],err_level);
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
    printf("Couldn't register old SIGTERM signal handler.\n");
  }

  if (signal(SIGABRT,  old_handler_abrt) == SIG_ERR) {
    printf("Couldn't register old SIGABRT signal handler.\n");
  }

  if (signal(SIGINT, old_handler_keyint) == SIG_ERR) {
    printf("Couldn't register old SIGINT signal handler.\n");
  }

#ifdef USE_FFTW
  if (alfa != (double) 0.0)
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
