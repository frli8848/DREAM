/***
*
* Copyright (C) 2006,2007,2008,2009,2014,2015,2016,2019 Fredrik Lingvall
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
#include "affinity.h"
#include "dream_error.h"

#include "dream.h"

/***
 *
 *  Parallel (threaded) matrix and sub-matrix copy.
 *
 ***/

//
// Globals
//

volatile int running;

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
  size_t r;
  size_t k;
  size_t len;
} DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//
void* smp_copy_p(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

/***
 *
 * Thread function.
 *
 ***/

void* smp_copy_p(void *arg)
{
  DATA D = *(DATA *)arg;
  size_t    line_start=D.line_start, line_stop=D.line_stop, n;
  double *RESTRICT A = (double*) D.A, *RESTRICT B = D.B;
  size_t A_M = D.A_M, B_M = D.B_M, r = D.r, k = D.k, len = D.len;

  for (n=line_start; n<line_stop; n++) {

    memcpy( &A[r+(k+n)*A_M], &B[0+n*B_M], len*sizeof(double));

    if (running==false) {
      mexPrintf("copy_p: thread for column %d -> %d bailing out!\n",line_start+1,line_stop);
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
 * Matlab (MEX) gateway function for copy_p.
 *
 ***/
extern void _main();

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double         *A,*B;
  sighandler_t   old_handler, old_handler_abrt, old_handler_keyint;
  size_t line_start, line_stop, A_M, A_N, B_M, B_N;
  size_t r_M, r_N, k_M, k_N;
  double         *r, *k;
  std::thread     *threads;
  size_t  thread_n, nthreads;
  DATA           *D;

  // Check for proper number of arguments


  if (nrhs != 4)
    mexErrMsgTxt("copy_p requires 4 input arguments!");

  if (nlhs > 0)
    mexErrMsgTxt("Too many output arguments for copy_p, none required!");

  A_M = mxGetM(prhs[0]);
  A_N = mxGetN(prhs[0]);
  A   = mxGetPr(prhs[0]);

  r_M = mxGetM(prhs[1]);
  r_N = mxGetN(prhs[1]);
  r   = mxGetPr(prhs[1]);

  k_M = mxGetM(prhs[2]);
  k_N = mxGetN(prhs[2]);
  k   = mxGetPr(prhs[2]);

  B_M = mxGetM(prhs[3]);
  B_N = mxGetN(prhs[3]);
  B   = mxGetPr(prhs[3]);

  // Check that arg 2.
  if ( r_M * r_N !=2 )
    mexErrMsgTxt("Argument 2 must be a 2 element vector!");

  if ( (size_t) r[1] < 1)
    mexErrMsgTxt("1st element of argument 2 must be > 1!");

  if (( (size_t) r[0] > A_M) || ( (size_t) r[1] > A_M))
    mexErrMsgTxt("One element of argument 2 exceeds row dimension of arg 1!");

  // Check that arg 3.
  if ( k_M * k_N !=2 )
    mexErrMsgTxt("Argument 3 must be a 2 element vector!");

  if ( (size_t) k[1] < 1)
    mexErrMsgTxt("1st element of argument 3 must be > 1!");

  if (( (size_t) k[0] > A_N) || ( (size_t) k[1] > A_N))
    mexErrMsgTxt("One element of argument 3 exceeds colomn dimension of arg 1!");

  // Check that arg 4.
  if (B_M < (size_t) (r[1]-r[0])+1 )
    mexErrMsgTxt("Argument 4 has not the number of rows as given by arg 2!");

  if (B_N != (size_t) (k[1]-k[0])+1 )
    mexErrMsgTxt("Argument 4 has not the number of columns as given by arg 3!");


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
  if (nthreads > A_N)
    nthreads = A_N;

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
  // Call the mem_copy subroutine.
  //

  running = true;

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
    D[thread_n].r = (size_t) (r[0]-1);
    D[thread_n].k = (size_t) (k[0]-1);
    D[thread_n].len = (size_t) (r[1]-r[0]+1);

    // Start the threads.
    threads[thread_n] = std::thread(smp_copy_p, &D[thread_n]);
    set_dream_thread_affinity(thread_n, nthreads, threads);
  }

  // Wait for all threads to finish.
  for (thread_n = 0; thread_n < nthreads; thread_n++)
    threads[thread_n].join();

  // Free memory.
  if (D) {
    free((void*) D);
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
