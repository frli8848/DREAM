/***
*
* Copyright (C) 2003,2004,2006,2007,2008,2009,2015,2016,2021 Fredrik Lingvall
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
#include <string>
#include <stdlib.h>
#include <math.h>
#include <thread>
#include <signal.h>
#include <uchar.h>
#include "mex.h"

#include "dream.h"
#include "dream_error.h"

#define printf mexPrintf

#define EQU 0
#define SUM 1
#define NEG 2

/***
 *
 *  Parallel (threaded) convolution.
 *
 ***/

//
// Globals
//

volatile int running;
volatile int in_place;
int mode = EQU;

//
// typedef:s
//

typedef struct
{
  size_t line_start;
  size_t line_stop;
  double *A;
  size_t A_M;
  size_t A_N;
  double *B;
  size_t B_M;
  size_t B_N;
  double *Y;
} DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//

void* smp_dream_conv_p(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

#ifdef __cplusplus
extern "C"
#endif
void conv(double *xr, size_t nx, double *yr, size_t ny, double *zr,
          int in_place, int mode);

/***
 *
 * Thread function.
 *
 ***/

void* smp_dream_conv_p(void *arg)
{
  DATA D = *(DATA *)arg;
  size_t    line_start=D.line_start, line_stop=D.line_stop, n;
  double *A = D.A, *B = D.B, *Y = D.Y;
  size_t A_M = D.A_M, B_M = D.B_M, B_N = D.B_N;

  // Do the convolution.

  for (n=line_start; n<line_stop; n++) {

    if (B_N > 1) // B is a matrix.
      conv( &A[0+n*A_M], A_M, &B[0+n*B_M], B_M, &Y[0+n*(A_M+B_M-1)],in_place,mode);
    else // B is a vector.
      conv( &A[0+n*A_M], A_M, B, B_M, &Y[0+n*(A_M+B_M-1)],in_place,mode);

    if (running == false) {
      mexPrintf("conv_p: thread for column %d -> %d bailing out!\n",line_start+1,line_stop);
      break;
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
 *  - Matlab (MEX) gateway function for CONV_P.
 *
 ***/

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *A,*B, *Y;
  sighandler_t   old_handler, old_handler_abrt, old_handler_keyint;
  size_t line_start, line_stop, A_M, A_N, B_M, B_N, n;
  int    buflen, is_set = false;
  DATA   *D;
  std::thread *threads;
  size_t  thread_n, nthreads;
  char   *the_str = NULL;

  in_place = false;

  // Check for proper inputs arguments.

  switch (nrhs) {

  case 0:
  case 1:
    mexErrMsgTxt("conv_p requires 2 to 4 input arguments!");
    break;

  case 2:
    if (nlhs > 1) {
      mexErrMsgTxt("Too many output arguments for conv_p!");
    }
    break;

  case 3:
    if (nlhs > 0) {
      mexErrMsgTxt("No output arguments required for conv_p in in-place operating mode!");
    }
    if (mxIsChar(prhs[2])) {
      mexErrMsgTxt("Arg 3 must be a matrix (not a string)");
    }
    break;

  case 4:
    if (mxIsChar(prhs[3])) { // 4th arg is a mode string.
      //the_str = (char*) mxGetChars(prhs[4]);
      buflen = mxGetM(prhs[3])*mxGetN(prhs[3])+1;
      the_str = (char*) malloc(buflen * sizeof(char));
      mxGetString(prhs[3], the_str, buflen); // Obsolete in Matlab 7.x ?

      // Valid strings are:
      //  '='  : In-place replace mode.
      //  '+=' : In-place add mode.
      //  '-=' : In-place sub mode.

      is_set = false;

      if (strcmp(the_str,"=") == 0) {
        mode = EQU;
        is_set = true;
      }

      if (strcmp(the_str,"+=") == 0) {
        mode = SUM;
        is_set = true;
      }

      if (strcmp(the_str,"-=") == 0) {
        mode = NEG;
        is_set = true;
      }

      if (is_set == false)
        mexErrMsgTxt("Non-valid string in arg 4!");

    }
    free(the_str);
    break;

  default:
    mexErrMsgTxt("conv_p requires 2 to 4 input arguments!");
    break;
  }

  A_M = mxGetM(prhs[0]);
  A_N = mxGetN(prhs[0]);
  A = mxGetPr(prhs[0]);

  B_M = mxGetM(prhs[1]);
  B_N = mxGetN(prhs[1]);
  B = mxGetPr(prhs[1]);

  // Check that arg 2.
  if ( B_M != 1 && B_N !=1 && B_N != A_N)
    mexErrMsgTxt("Argument 2 must be a vector or a matrix with the same number of rows as arg 1!");

  if (  B_M == 1 || B_N == 1 ) { // B is a vector.
    B_M = B_M*B_N;
    B_N = 1;
  }

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

  // nthreads can't be larger then the number of columns in the A matrix.
  if (nthreads > A_N) {
    nthreads = A_N;
  }

  if (nrhs == 2) { // Test for in-place/normal mode.

    //
    // Normal mode.
    //

    in_place = false;
    plhs[0] = mxCreateDoubleMatrix(A_M+B_M-1, A_N, mxREAL);
    Y = mxGetPr(plhs[0]);

  } else {

    //
    // in-place mode.
    //

    in_place = true;

    if ( mxGetM(prhs[2]) != A_M+B_M-1)
      mexErrMsgTxt("Wrong number of rows in argument 3!");

    if ( mxGetN(prhs[2]) != A_N)
      mexErrMsgTxt("Wrong number of columns in argument 3!");

    Y = mxGetPr(prhs[2]);
  }

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
  // Call the CONV subroutine.
  //

  running = true;

  if (nthreads>1) { // Use threads

    // Allocate local data.
    D = (DATA*) malloc(nthreads*sizeof(DATA));
    if (!D)
      mexErrMsgTxt("Failed to allocate memory for thread data!");

    // Allocate mem for the threads.
    threads = new std::thread[nthreads]; // Init thread data.
    if (!threads)
      mexErrMsgTxt("Failed to allocate memory for threads!");

    for (thread_n = 0; thread_n < nthreads; thread_n++) {

      line_start = thread_n * A_N/nthreads;
      line_stop =  (thread_n+1) * A_N/nthreads;

      // Init local data.
      D[thread_n].line_start = line_start; // Local start index;
      D[thread_n].line_stop = line_stop; // Local stop index;
      D[thread_n].A = A;
      D[thread_n].A_M = A_M;
      D[thread_n].A_N = A_N;
      D[thread_n].B = B;
      D[thread_n].B_M = B_M;
      D[thread_n].B_N = B_N;
      D[thread_n].Y = Y;

      // Start the threads.
      threads[thread_n] = std::thread(smp_dream_conv_p, &D[thread_n]);

    } // for (thread_n = 0; thread_n < nthreads; thread_n++)

    // Wait for all threads to finish.
    for (thread_n = 0; thread_n < nthreads; thread_n++)
        threads[thread_n].join();

    // Free memory.
    if (D)
      free((void*) D);

  } else {			// Do not use threads

    if (B_N > 1) {		// B is a matrix.
      for (n=0; n<A_N; n++) {

        conv( &A[0+n*A_M], A_M, &B[0+n*B_M], B_M, &Y[0+n*(A_M+B_M-1)],in_place,mode);

        if (running==false) {
          printf("conv_p: bailing out!\n");
          break;
        }
      }
    } else {			// B is a vector.
      for (n=0; n<A_N; n++) {

        conv( &A[0+n*A_M], A_M, B, B_M, &Y[0+n*(A_M+B_M-1)],in_place,mode);

        if (running==false) {
          printf("conv_p: bailing out!\n");
          break;
        }
      }
    }
  }

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
    // This may kill Matlab 7 release 14 (bug in Matlab).
    mexErrMsgTxt("CTRL-C pressed!\n"); // Bail out.
  }

  return;
}
