/***
*
* Copyright (C) 2003,2004,2006,2007,2008,2009,2014,2015,2016,2019 Fredrik Lingvall
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
#include <stdlib.h>
#include <math.h>
#include <thread>
#include <signal.h>
#include <uchar.h>
#include "mex.h"
#include "dream.h"
#include "affinity.h"
#include "dream_error.h"

#define SINGLE 0
#define MULTIPLE 1

#if defined(_MSC_VER) || defined(__LCC__)
#include "msvc_rint.h"
#endif

/***
 *
 * (Monostatic) Synthetic Aperture Focusing Technique - Parallel version.
 *
 ***/

//
// Globals
//

int running;

//
// typedef:s
//

typedef struct
{
  int    no;
  int    start;
  int    stop;
  double *ro;
  double dt;
  int delay_method;
  double *delay;
  double cp;
  double a;
  int    K;
  int    L;
  double *r_trans;
  double *B;
  double *Bsaft;
} DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//
void* smp_process(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

/***
 *
 * Thread function.
 *
 ***/

void* smp_process(void *arg)
{
  size_t n, l;
  int  k_shift;			// Can be negative!
  DATA D = *(DATA *)arg;
  double *RESTRICT Bsaft = D.Bsaft;
  double *RESTRICT B = D.B;
  double dt=D.dt;
  int  no=D.no, K=D.K, L=D.L;
  double *RESTRICT delay=D.delay, *ro=D.ro, cp=D.cp;
  double xo, yo, zo, r_xy;
  size_t start=D.start, stop=D.stop;
  double *RESTRICT r_trans = D.r_trans, a = D.a;
  double x_trans, y_trans, z_trans, d, z, z_prim, tmp;

  if (D.delay_method == SINGLE) {
    for (n=start; n<stop; n++) {
      xo = ro[n];
      yo = ro[n+1*no];
      zo = ro[n+2*no];

      // The SAFT loop.
      for (l=0; l< (size_t) L; l++) {

        x_trans = r_trans[l];
        y_trans = r_trans[l+1*L];
        z_trans = r_trans[l+2*L];

        // Horizontal distance between transducer and observation point.
        r_xy = sqrt( (x_trans - xo)*(x_trans - xo) + (y_trans - yo)*(y_trans - yo) );

        // if r_xy is inside the synthetic aperture.
        if (r_xy <= a/2) {

          // Vertical distance between transducer and observation point.
          z = zo - z_trans;

          z_prim = sqrt(z*z + r_xy*r_xy);
          tmp = (2*z_prim/cp * 1e3 - 2*delay[0])/dt;

          // Better to round just one time!
          k_shift = (int) rint(tmp);

          // Rounding err.
          d =  tmp - ((double) k_shift);

          // Linear interpolation.
          if ((k_shift+1 < K) && (k_shift-1 > 0) ) {
            if (d >=0) {
              Bsaft[n] += (1.0 - d) * B[k_shift     + l*K];
              Bsaft[n] += d * B[(k_shift+1) + l*K];
            }
            else {
              Bsaft[n] -= d * B[k_shift     + l*K];
              Bsaft[n] += (1.0 + d) * B[(k_shift-1) + l*K];
            }
          } // if
        }
      }

      if (!running) {
        mexPrintf("Thread for observation points %d -> %d bailing out!\n",start+1,stop);
        return(NULL);
      }

    }

  } else { // MULTIPLE delays.

    for (n=start; n<stop; n++) {
      xo = ro[n];
      yo = ro[n+1*no];
      zo = ro[n+2*no];

      // The SAFT loop.
      for (l=0; l < (size_t) L; l++) {

        x_trans = r_trans[l];
        y_trans = r_trans[l+1*L];
        z_trans = r_trans[l+2*L];

        // Horizontal distance between transducer and observation point.
        r_xy = sqrt( (x_trans - xo)*(x_trans - xo) + (y_trans - yo)*(y_trans - yo) );

        // if r_xy is inside the synthetic aperture.
        if (r_xy <= a/2) {

          // Vertical distance between transducer and observation point.
          z = zo - z_trans;

          z_prim = sqrt(z*z + r_xy*r_xy);
          tmp = (2*z_prim/cp * 1e3 - 2*delay[n])/dt;

          // Better to round just one time!
          k_shift = (int) rint(tmp);

          // Rounding err.
          d =  tmp - ((double) k_shift);

          // Linear interpolation.
          if ((k_shift+1 < K) && (k_shift-1 > 0) ) {
            if (d >=0) {
              Bsaft[n] += (1.0 - d) * B[k_shift     + l*K];
              Bsaft[n] += d * B[(k_shift+1) + l*K];
            }
            else {
              Bsaft[n] -= d * B[k_shift     + l*K];
              Bsaft[n] += (1.0 + d) * B[(k_shift-1) + l*K];
            }
          } // if
        }
      }

      if (!running) {
        mexPrintf("Thread for observation points %d -> %d bailing out!\n",start+1,stop);
        return(NULL);
      }
    }
  }

  return(NULL);
}

/***
 *
 * Signal handlers.
 *
 ***/

void sighandler(int signum) {
  //printf("Caught signal SIGTERM.\n");
  running = false;
}

void sig_abrt_handler(int signum) {
  //printf("Caught signal SIGABRT.\n");
}

void sig_keyint_handler(int signum) {
  //printf("Caught signal SIGINT.\n");
}

/***
 *
 * saft_p  - Matlab (MEX) gateway function for SAFT_P.
 *
 ***/

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *RESTRICT s_par,*m_par;
  size_t  K, L, no;
  double *RESTRICT B, *RESTRICT Bsaft, cp , a, dt, *RESTRICT ro;
  double *RESTRICT r_trans, *RESTRICT delay;
  sighandler_t   old_handler, old_handler_abrt, old_handler_keyint;
  std::thread *threads;
  unsigned int thread_n, nthreads;
  size_t start, stop;
  DATA   *D;

  // Check for proper number of arguments

  if (nrhs != 7) {
    mexErrMsgTxt("saft requires 7 input arguments!");
  }
  else
    if (nlhs > 1) {
      mexErrMsgTxt("Too many output arguments for saft!");
    }

  //
  // B-scan.
  //

  K = mxGetM(prhs[0]); // Number of temporal samples.
  L = mxGetN(prhs[0]); // Number of transducer positions.
  B = mxGetPr(prhs[0]);

  //
  // Transducer position matrix.
  //

  // Check that arg 2 is a L x 3 is a matrix.
  if (!(mxGetM(prhs[1]) == L && mxGetN(prhs[1]) == 3))
    mexErrMsgTxt("Arg 2 must be a (number of transducer positions) x 3 matrix!");

  r_trans = mxGetPr(prhs[1]);

  //
  // Start point of impulse response vector ([us]).
  //

  // Check that arg 3 is a scalar (or vector).
  if ( (mxGetM(prhs[2]) * mxGetN(prhs[2]) !=1) && ((mxGetM(prhs[2]) * mxGetN(prhs[2])) != L))
    dream_err_msg("Argument 3 (delay(s)) must be a scalar or a vector with a length equal to the number of A-scans!");

  delay = mxGetPr(prhs[2]);

  //
  // Sampling parameter.
  //

  // Check that arg 4 is a scalar.
  if (!((mxGetM(prhs[3])==1 && mxGetN(prhs[3])==1)))
    dream_err_msg("Argument 4 (temporal sampling period) must be a scalar");

  s_par = mxGetPr(prhs[3]);
  dt    = s_par[0]; // Temporal discretization size (= 1/sampling freq).

  //
  // Material parameter.
  //

  // Check that arg 5 is a scalar.
  if (!((mxGetM(prhs[4])==1 && mxGetN(prhs[4])==1)))
    dream_err_msg("Argument 5 (sound speed) must be a scalar!");

  m_par = mxGetPr(prhs[4]);
  cp    = m_par[0]; // Sound speed.

  //
  // Observation point matrix.
  //

  // Check that arg 6 is a (number of observation points) x 3 matrix.
  if (!mxGetN(prhs[5])==3)
    dream_err_msg("Argument 6 must be a (number of observation points) x 3 matrix!");

  no = mxGetM(prhs[5]); // Number of observation points.
  ro = mxGetPr(prhs[5]);


  //
  // Synthetic Aperture
  //

  // Check that arg 7 is scalar.
  if ((mxGetM(prhs[6])!=1) && (mxGetN(prhs[6])!=1))
    mexErrMsgTxt("Argument 7 must be a scalar!");

  a = mxGetScalar(prhs[6]);

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11)
  nthreads = std::thread::hardware_concurrency();

  // Read OMP_NUM_THREADS env var
  if(const char* env_p = std::getenv("OMP_NUM_THREADS")) {
    unsigned int omp_threads = std::stoul(env_p);
    if (omp_threads < nthreads) {
      nthreads = omp_threads;
    }
  }

  // nthreads can't be larger then the number of observation points.
  if (nthreads > (unsigned int) no) {
    nthreads = no;
  }

  //
  // Create an output matrix for the processed image.
  //

  plhs[0] = mxCreateDoubleMatrix(no,1,mxREAL);
  Bsaft = mxGetPr(plhs[0]);

  //
  // Register signal handlers.
  //

  if ((old_handler = signal(SIGTERM, &sighandler)) == SIG_ERR) {
    printf("Couldn't register signal handler.\n");
  }

  if ((old_handler_abrt=signal(SIGABRT, &sighandler)) == SIG_ERR) {
    printf("Couldn't register signal handler.\n");
  }

  if ((old_handler_keyint=signal(SIGINT, &sighandler)) == SIG_ERR) {
    printf("Couldn't register signal handler.\n");
  }

  //
  // Call the SAFT subroutine.
  //

  running = true;

  // Allocate local data.
  D = (DATA*) malloc(nthreads*sizeof(DATA));

  // Allocate mem for the threads.
  threads = new std::thread[nthreads]; // Init thread data.

  for (thread_n = 0; thread_n < nthreads; thread_n++) {

    start = thread_n * no/nthreads;
    stop =  (thread_n+1) * no/nthreads;

    // Init local data.
    D[thread_n].start = start; // Local start index;
    D[thread_n].stop = stop; // Local stop index;
    D[thread_n].no = no;
    D[thread_n].ro = ro;
    D[thread_n].cp = cp;
    D[thread_n].dt = dt;

    if (mxGetM(prhs[3]) * mxGetN(prhs[3]) == 1)
      D[thread_n].delay_method = SINGLE; // delay is a scalar.
    else
      D[thread_n].delay_method = MULTIPLE; // delay is a vector.

    D[thread_n].delay = delay;
    D[thread_n].B = B;
    D[thread_n].Bsaft = Bsaft;
    D[thread_n].a = a;
    D[thread_n].K = K;
    D[thread_n].L = L;
    D[thread_n].r_trans = r_trans;

    // Start the threads.
    threads[thread_n] = std::thread(smp_process, &D[thread_n]);
    set_dream_thread_affinity(thread_n, nthreads, threads);
  }

  // Wait for all threads to finish.
  for (thread_n = 0; thread_n < nthreads; thread_n++)
    threads[thread_n].join();

  // Free memory.
  free((void*) D);

  //
  // Restore old signal handlers.
  //

  if (signal(SIGTERM, old_handler) == SIG_ERR) {
    printf("Couldn't register old signal handler.\n");
  }

  if (signal(SIGABRT,  old_handler_abrt) == SIG_ERR) {
    printf("Couldn't register signal handler.\n");
  }

  if (signal(SIGINT, old_handler_keyint) == SIG_ERR) {
    printf("Couldn't register signal handler.\n");
  }

  if (!running) {
    mexErrMsgTxt("CTRL-C pressed!\n"); // Bail out.
  }

  return;
}
