/***
*
* Copyright (C) 2003,2004,2006,2007,2008,2009,2015,2016,2021,2022,2025 Fredrik Lingvall
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

#include <csignal>
#include <thread>

#include "mex.h"

#include "dream.h"
#include "affinity.h"
#include "conv.h"

/***
 *
 *  Parallel (threaded) convolution.
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
  size_t col_start;
  size_t col_stop;
  double *A;
  size_t A_M;
  size_t A_N;
  double *B;
  size_t B_M;
  size_t B_N;
  double *Y;
  ConvMode conv_mode;
} DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//

void* smp_dream_conv(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

/***
 *
 * Thread function.
 *
 ***/

void* smp_dream_conv(void *arg)
{
  DATA D = *(DATA *)arg;
  size_t    col_start=D.col_start, col_stop=D.col_stop;
  double *A = D.A, *B = D.B, *Y = D.Y;
  size_t A_M = D.A_M, B_M = D.B_M, B_N = D.B_N;
  ConvMode conv_mode = D.conv_mode;

  // Do the convolution.

  for (dream_idx_type n=col_start; n<col_stop; n++) {

    if (B_N > 1) // B is a matrix.
      conv( &A[0+n*A_M], A_M, &B[0+n*B_M], B_M, &Y[0+n*(A_M+B_M-1)], conv_mode);
    else // B is a vector.
      conv( &A[0+n*A_M], A_M, B, B_M, &Y[0+n*(A_M+B_M-1)], conv_mode);

    if (running == false) {
      std::cout << "conv_p: thread for column " << col_start+1 << " -> " << col_stop << "bailing out!\n";
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
  double *A = nullptr,*B = nullptr, *Y = nullptr;
  sighandler_t   old_handler, old_handler_abrt, old_handler_keyint;
  size_t col_start, col_stop, A_M, A_N, B_M, B_N;
  bool is_set = false;
  DATA *D = nullptr;
  std::thread *threads;
  size_t  thread_n, nthreads;
  ConvMode conv_mode = ConvMode::equ;

  // Check for proper inputs arguments.

  switch (nrhs) {

  case 0:
  case 1:
    dream_err_msg("conv_p requires 2 to 4 input arguments!");
    break;

  case 2:
    if (nlhs > 1) {
      dream_err_msg("Too many output arguments for conv_p!");
    }
    break;

  case 3:
    if (nlhs > 0) {
      dream_err_msg("No output arguments required for conv_p in in-place operating mode!");
    }
    if (mxIsChar(prhs[2])) {
      dream_err_msg("Arg 3 must be a matrix (not a string)");
    }
    break;

  case 4:
    if (mxIsChar(prhs[3])) { // 4th arg is a mode string.

      char *str = mxArrayToString(prhs[3]);
      std::string ip_mode(str);
      mxFree(str);

      // Valid strings are:
      //  '='  : In-place replace mode.
      //  '+=' : In-place add mode.
      //  '-=' : In-place sub mode.

      is_set = false;

      if (ip_mode.compare("=") == 0) {
        conv_mode=ConvMode::equ;
        is_set = true;
      }

      if (ip_mode.compare("+=") == 0) {
        conv_mode=ConvMode::sum;
        is_set = true;
      }

      if (ip_mode.compare("-=") == 0) {
        conv_mode=ConvMode::neg;
        is_set = true;
      }

      if (is_set == false) {
        dream_err_msg("Non-valid string in arg 4!");
      }
    }
    break;

  default:
    dream_err_msg("conv_p requires 2 to 4 input arguments!");
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
    dream_err_msg("Argument 2 must be a vector or a matrix with the same number of rows as arg 1!");

  if (  B_M == 1 || B_N == 1 ) { // B is a vector.
    B_M = B_M*B_N;
    B_N = 1;
  }

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11)
  nthreads = std::thread::hardware_concurrency();

  // Read DREAM_NUM_THREADS env var
  if(const char* env_p = std::getenv("DREAM_NUM_THREADS")) {
    unsigned int dream_threads = std::stoul(env_p);
    if (dream_threads < nthreads) {
      nthreads = dream_threads;
    }
  }

  // nthreads can't be larger then the number of columns in the A matrix.
  if (nthreads > A_N) {
    nthreads = A_N;
  }

  if (nthreads < 1) {
    nthreads = 1;
  }

  if (nrhs == 2) { // Test for in-place/normal mode.

    //
    // Normal mode.
    //

    plhs[0] = mxCreateDoubleMatrix(A_M+B_M-1, A_N, mxREAL);
    Y = mxGetPr(plhs[0]);

  } else {

    //
    // in-place mode.
    //

    if ( mxGetM(prhs[2]) != A_M+B_M-1)
      dream_err_msg("Wrong number of rows in argument 3!");

    if ( mxGetN(prhs[2]) != A_N)
      dream_err_msg("Wrong number of columns in argument 3!");

    Y = mxGetPr(prhs[2]);
  }

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
  // Call the CONV subroutine.
  //

  running = true;

  // Allocate local data.
  D = (DATA*) malloc(nthreads*sizeof(DATA));
  if (!D) {
    dream_err_msg("Failed to allocate memory for thread data!");
  }

  // Allocate mem for the threads.
  threads = new std::thread[nthreads]; // Init thread data.
  if (!threads)
    dream_err_msg("Failed to allocate memory for threads!");

  for (thread_n = 0; thread_n < nthreads; thread_n++) {

    col_start = thread_n * A_N/nthreads;
    col_stop =  (thread_n+1) * A_N/nthreads;

    // Init local data.
    D[thread_n].col_start = col_start; // Local start index;
    D[thread_n].col_stop = col_stop; // Local stop index;
    D[thread_n].A = A;
    D[thread_n].A_M = A_M;
    D[thread_n].A_N = A_N;
    D[thread_n].B = B;
    D[thread_n].B_M = B_M;
    D[thread_n].B_N = B_N;
    D[thread_n].Y = Y;
    D[thread_n].conv_mode = conv_mode;

    if (nthreads > 1) {
      // Start the threads.
      threads[thread_n] = std::thread(smp_dream_conv, &D[thread_n]);
      set_dream_thread_affinity(thread_n, nthreads, threads);
    } else {
      smp_dream_conv(&D[0]);
    }
  }

  if (nthreads > 1) {
    // Wait for all threads to finish.
    for (thread_n = 0; thread_n < nthreads; thread_n++) {
      threads[thread_n].join();
    }
  }

  // Free memory.
  if (D) {
    free((void*) D);
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

  return;
}
