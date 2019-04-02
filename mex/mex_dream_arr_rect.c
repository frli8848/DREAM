/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2014,2016 Fredrik Lingvall
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
#include "dream_arr_rect.h"
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
 * Gateway function for dream_arr_rect.
 *
 * This is  a MEX file for MATLAB.
 *
 ***/

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *RESTRICT ro, *RESTRICT geom_par, *RESTRICT s_par, *RESTRICT m_par;
  double *RESTRICT steer_par;
  char   apod_met[50],foc_met[50],steer_met[50];
  int    buflen;
  double xo,yo,zo,a,b,dx,dy,dt;
  dream_idx_type    nt,no,n;
  double param=0, *RESTRICT delay, v, cp, alfa;
  int    isize;
  double *RESTRICT gx, *RESTRICT gy, *RESTRICT gz;
  int    ifoc=0;
  double focal=0, *RESTRICT ud_focal=NULL;
  int    ister=0;
  double theta=0, phi=0, *RESTRICT apod=NULL;
  int    iweight=0,iapo=0;
  double *RESTRICT h, *err_p;
  int    err_level=STOP, err=NONE, out_err = NONE , set = FALSE;
  char   err_str[50];
  sighandler_t old_handler, old_handler_abrt, old_handler_keyint;

  // Check for proper number of arguments

  if (!((nrhs == 13) || (nrhs == 14))) {
    dream_err_msg("dream_arr_rect requires 13 or 14 input arguments!");
  }
  else
    if (nlhs > 2) {
      dream_err_msg("Too many output arguments for dream_arr_rect!");
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
    dream_err_msg("Argument 2 must be a vector of length 2!");

  geom_par = mxGetPr(prhs[1]);
  a = geom_par[0];		// min x-pos.
  b = geom_par[1];		// max x-pos.

  //
  // Grid function (position vectors of the elements).
  //

  isize = (int) mxGetM(prhs[2]); // Number of elementents in the array.
  if (mxGetN(prhs[2]) !=3 )
    dream_err_msg("Argument 3  must a (number of array elements) x 3 matrix!");

  gx    = mxGetPr(prhs[2]);	// First column in the matrix.
  gy    = gx + isize;		// Second column in the matrix.
  gz    = gy + isize;		// Third column in the matrix.

  //
  // Temporal and spatial sampling parameters.
  //


  // Check that arg 4 is a 4 element vector
  if (!((mxGetM(prhs[3])==4 && mxGetN(prhs[3])==1) || (mxGetM(prhs[3])==1 && mxGetN(prhs[3])==4)))
    dream_err_msg("Argument 4 must be a vector of length 4!");

  s_par = mxGetPr(prhs[3]);
  dx    = s_par[0];		// Spatial x discretization size.
  dy    = s_par[1];		// Spatial dy iscretization size.
  dt    = s_par[2];		// Temporal discretization size (= 1/sampling freq).
  nt    = (dream_idx_type) s_par[3];	// Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //

  // Check that arg 5 is a scalar.
  if ( (mxGetM(prhs[4]) * mxGetN(prhs[4]) !=1) && ((mxGetM(prhs[4]) * mxGetN(prhs[4])) != no))
    dream_err_msg("Argument 5 must be a scalar or a vector with a length equal to the number of observation points!");

  delay = mxGetPr(prhs[4]);

  //
  // Material parameters
  //

  // Check that arg 6 is a 3 element vector.
  if (!((mxGetM(prhs[5])==3 && mxGetN(prhs[5])==1) || (mxGetM(prhs[5])==1 && mxGetN(prhs[5])==3)))
    dream_err_msg("Argument 6 must be a vector of length 3!");

  m_par = mxGetPr(prhs[5]);
  v     = m_par[0]; // Normal velocity of transducer surface.
  cp    = m_par[1]; // Sound speed.
  alfa  = m_par[2]; // Attenuation coefficient [dB/(cm MHz)].

  //
  // Focusing parameters.
  //

  //  ifoc = 1 - no foc, 2 foc x ,3 foc y, 4 foc xy (del=fsqrt(x*x+y*y)), 5 focx+focy.

  if (nrhs >= 7) {

    if (!mxIsChar(prhs[6]))
      dream_err_msg("Argument 7 must be a string");

    buflen = (mxGetM(prhs[6]) * mxGetN(prhs[6]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[6],foc_met,buflen);

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

      if (mxGetM(prhs[7]) * mxGetN(prhs[7]) != isize ) {
        mexPrintf("The time delay vector (argument 8) for user defined ('ud') focusing\n") ;
        dream_err_msg("delays must have the same length as the number of array elements.!");
      }
      ud_focal = mxGetPr(prhs[7]);
    }
    else {

      // Check that arg 8 is a scalar.
      if (mxGetM(prhs[7]) * mxGetN(prhs[7]) !=1 )
        dream_err_msg("Argument 8  must be a scalar!");

      // Focal point (in mm).
      focal = mxGetScalar(prhs[7]);
    }

    if (set == FALSE)
      dream_err_msg("Unknown focusing method!");

  } else
    ifoc = 1;


  //
  // Beam steering.
  //

  // Beam steering: ister = 1 - no steering, 2 steer ph=ax ,3 steer y ph=by, 4 steer xy ph=ax+by.

  if (nrhs >= 9) {

    if (!mxIsChar(prhs[8]))
      dream_err_msg("Argument 9 must be a string");

    buflen = (mxGetM(prhs[8]) * mxGetN(prhs[8]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[8],steer_met,buflen);

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

    // Check that arg 10 is a 2 element vector
    if (!((mxGetM(prhs[9])==2 && mxGetN(prhs[9])==1) || (mxGetM(prhs[9])==1 && mxGetN(prhs[9])==2)))
      dream_err_msg("Argument 10 must be a vector of length 2!");

    steer_par = mxGetPr(prhs[9]);
    theta  = steer_par[0];		// Angle in x-direction.
    phi    = steer_par[1];		// Angle in y-direction.

  } else
    ister = 1;

  //
  // Apodization.
  //

  // iweight = 1 - no apodization, 2  apodization.
  // iapo = 0 - user defined, 1 traingle, 2 Gauss, 3 raised cosine, 4 simply supported, 5 clamped.

  if (nrhs >= 11) {

    if (!mxIsChar(prhs[10]))
      dream_err_msg("Argument 11 must be a string");

    buflen = (mxGetM(prhs[10]) * mxGetN(prhs[10]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[10],apod_met,buflen);

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
      if (mxGetM(prhs[11]) * mxGetN(prhs[11]) != isize)
        dream_err_msg("The length of argument 12 (apodization vector) must be the same as the number of array elements!");

      apod = mxGetPr(prhs[11]);
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
      dream_err_msg("Unknown apodization!");

    // Parameter for raised cos and Gaussian apodization functions.
    if (mxGetM(prhs[12]) * mxGetN(prhs[12]) !=1 )
      dream_err_msg("Argument 13 must be a scalar");

    param = mxGetScalar(prhs[12]);
  }
  else
    iweight = 1;


  // Error reporting.
  if (nrhs == 14) {

    if (!mxIsChar(prhs[13]))
      dream_err_msg("Argument 14 must be a string");

    buflen = (mxGetM(prhs[13]) * mxGetN(prhs[13]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[13],err_str,buflen);

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
  if (alfa != (double) 0.0)
    att_init(nt,1);
#endif

  if (ifoc != 6) {

    if (mxGetM(prhs[4]) * mxGetN(prhs[4]) == 1) {
      for (n=0; n<no; n++) {
        xo = ro[n];
        yo = ro[n+1*no];
        zo = ro[n+2*no];

        err = dream_arr_rect(xo,yo,zo,a,b,dx,dy,dt,nt,delay[0],v,cp,alfa,isize,gx,gy,gz,
                             ifoc,focal,ister,theta,phi,apod,iweight,iapo,param,&h[n*nt],err_level);
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

        err = dream_arr_rect(xo,yo,zo,a,b,dx,dy,dt,nt,delay[n],v,cp,alfa,isize,gx,gy,gz,
                             ifoc,focal,ister,theta,phi,apod,iweight,iapo,param,&h[n*nt],err_level);
        if (err != NONE)
          out_err = err;

        if (!running) {
          break; // CTRL-C pressed.
        }

      }
    }
  }
  else { // User defined focusing.

    if (mxGetM(prhs[4]) * mxGetN(prhs[4]) == 1) {
      for (n=0; n<no; n++) {
        xo = ro[n];
        yo = ro[n+1*no];
        zo = ro[n+2*no];

        err = dream_arr_rect_ud(xo,yo,zo,a,b,dx,dy,dt,nt,delay[0],v,cp,alfa,isize,gx,gy,gz,
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

        err = dream_arr_rect_ud(xo,yo,zo,a,b,dx,dy,dt,nt,delay[n],v,cp,alfa,isize,gx,gy,gz,
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
