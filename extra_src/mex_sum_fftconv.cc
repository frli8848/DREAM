/***
*
* Copyright (C) 2006,2007,2008,2009,2012,2016,2019,2021,2022 Fredrik Lingvall
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

// FIXME: We cannot use FFTW for mex-files due to lib conflicts!

#include <iostream>
#include <csignal>
#include <cstring>
#include <thread>
#include <complex> // C++

#include <fftw3.h>

#include "dream.h"
#include "affinity.h"
#include "fftconv.h"

#include "mex.h"

/***
 *
 *  Parallel (threaded) FFTW based convolution and summation.
 *
 ***/

//
// Globals
//
volatile int running;
int plan_method = 4; // Default to ESTIMATE method.

//
// typedef:s
//

typedef struct
{
  size_t line_start;
  size_t line_stop;
  int L;
  double **H;
  size_t H_M;
  size_t H_N;
  double *U;
  size_t U_M;
  size_t U_N;
  double *Y;
  FFT *fft;
  ConvMode conv_mode;
} DATA;


typedef void (*sighandler_t)(int);

//
// Function prototypes.
//
void* smp_dream_sum_fftconv(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

/***
 *
 * Thread function.
 *
 ***/

void* smp_dream_sum_fftconv(void *arg)
{
  DATA D = *(DATA *)arg;
  size_t line_start=D.line_start, line_stop=D.line_stop, n;
  double **H = D.H, *U = D.U, *Y = D.Y;
  size_t H_M = D.H_M, U_M = D.U_M; //, U_N = D.U_N;
  size_t L = D.L;
  FFT fft = *D.fft;
  ConvMode conv_mode = D.conv_mode;

  dream_idx_type fft_len = H_M + U_M - 1;

  // Input vectors.
  FFTVec a_v(fft_len);
  FFTVec b_v(fft_len);
  FFTVec c_v(fft_len);
  double *a = a_v.get(), *b = b_v.get(), *c  = c_v.get();

  // Fourier Coefficients.
  FFTCVec af_v(fft_len);
  FFTCVec bf_v(fft_len);
  FFTCVec cf_v(fft_len);
  std::complex<double> *af  = af_v.get(), *bf  = bf_v.get(), *cf  = cf_v.get();

  //
  // Do the convolution.
  //

  // NB. U is always a matrix here

  for (n=line_start; n<line_stop; n++) {

    add_fftconv(fft,
                H, L, H_M,
                n,
                U, U_M, // U must be U_M x L
                &Y[0+n*(H_M+U_M-1)],
                a, b, c, af, bf, cf,
                conv_mode);

    if (running==false) {
      std::cout << "sum_fftconv: thread for column " << line_start+1 << " -> " << line_stop << " bailing out!\n";
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
 * Matlab (mex) gateway function for SUM_FFTCONV.
 *
 ***/

extern void _main();

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *U, *Y = NULL;
  sighandler_t  old_handler, old_handler_abrt, old_handler_keyint;
  size_t line_start, line_stop, H_M, H_N, H_L, U_M, U_N, n, k;
  std::thread *threads;
  size_t      thread_n, nthreads;
  DATA   *D;
  // Input vectors (only used for creating fftw plans).
  double *a  = NULL, *b  = NULL;
  double *c  = NULL;
  // Fourier Coefficients (only used for creating fftw plans).
  std::complex<double> *af  = NULL, *bf  = NULL, *cf  = NULL;
  size_t fft_len, return_wisdom = false, load_wisdom = false;
  char *wisdom_str = NULL;
  int buflen;
  const mwSize *dv; // This don't compile in XP with MinGW (not defined in matrix.h?!).
                    // Needs Matlab R2006b (eg. R2006b's matlab/extern/include/matrix.h)
  //const long *dv; // This gives a runtime error.
  int dims;

  //
  // Set the method which fftw computes plans
  //

  // If we want to save a plan (in the second output arg)
  // then use the more time-consuming MEAUSURE method.
  if (nlhs == 2) {
    plan_method = 3; // 3 = MEASURE.
  }

  // Check for proper inputs arguments.

  switch (nrhs) {

  case 0:
  case 1:
    mexErrMsgTxt("sum_fftconv requires 2 to 3 input arguments!");
    break;

  case 2:
    if (nlhs > 1) {
      mexErrMsgTxt("Too many output arguments for sum_fftconv!");
    }
    if (nlhs == 2)
      return_wisdom = true;

    break;

  case 3:
    if (mxIsChar(prhs[2])) { // 3rd arg is a string.
      //wisdom_str = (char*) mxGetChars(prhs[3]);
      buflen = mxGetM(prhs[2])*mxGetN(prhs[3])+1;
      wisdom_str = (char*) fftw_malloc(buflen * sizeof(char));
      mxGetString(prhs[2], wisdom_str, buflen); // Obsolete in Matlab 7.x

      //
      // If 3:th arg is a string then only a wisdom string is valid.
      //

      if (strcmp("wisdom",wisdom_str) < 0) {
        mexErrMsgTxt("The string in arg 3 do not seem to be in fftw wisdom format!");
      }
      else {
        load_wisdom = true;
      }
    } else { // 3d arg not a string then assume in-place mode.
      fftw_forget_wisdom(); // Clear wisdom history (a new wisdom will be created below).
      if (nlhs > 0) {
        mexErrMsgTxt("No output arguments required for fftconv in in-place operating mode!");
      }
    }
    break;

  case 4:
    if (mxIsChar(prhs[3])) { // 4th arg is a string.
      //wisdom_str = (char*) mxGetChars(prhs[3]);
      buflen = mxGetM(prhs[3])*mxGetN(prhs[3])+1;
      wisdom_str = (char*) fftw_malloc(buflen * sizeof(char));
      mxGetString(prhs[3], wisdom_str, buflen); // Obsolete in Matlab 7.x ?

      if (!strcmp("wisdom",wisdom_str)) {
        mexErrMsgTxt("The string in arg 4 do not seem to be in fftw wisdom format!");

      }
      else {
        load_wisdom = true;
      }
    }
    else { // 4th arg not a string
      mexErrMsgTxt("Argument 4 do not seem to be in fftw wisdom string format!");

    }
    break;

  default:
    mexErrMsgTxt("sum_fftconv requires 2 to 4 input arguments!");

    break;
  }

  dims = mxGetNumberOfDimensions(prhs[0]);
  if (dims != 3)
    mexErrMsgTxt("Argument 1 should be a 3D Matrix\n");

  dv = mxGetDimensions(prhs[0]);
  H_M = dv[0];
  H_N = dv[1];
  H_L = dv[2];

  // Store pointers to the L A-matrices in a vector.
  double **H = (double**) malloc(H_L*sizeof(double*));
  double* tmp_p = mxGetPr(prhs[0]);
  for (size_t k=0; k<H_L; k++) {
    H[k] = &(tmp_p[H_M*H_N*k]);
  }

  U_M = mxGetM(prhs[1]);
  U_N = mxGetN(prhs[1]);
  U = mxGetPr(prhs[1]);

  if (H_L != U_N) {
    mexErrMsgTxt("3rd dimension of arg 1 must match the number of columns in arg 2\n");
  }

  if (U_M == 1 || U_N == 1 ) { // U is a vector.
    U_M = U_M*U_N;
    U_N = 1;
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

  // We cannot have more threads than the number of observation points.
  if (nthreads > H_N)
    nthreads = H_N;

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
  // Init the FFTW plans.
  //

  if(load_wisdom) {
    if (!fftw_import_wisdom_from_string(wisdom_str)) {
      mexErrMsgTxt("Failed to load fftw wisdom!");

    } else
      fftw_free(wisdom_str); // Clean up.
  }

  // NB We call the FFTW planners only from the main thread once
  // so no need to lock it with a mutex
  FFT fft(fft_len, nullptr, plan_method);

  //
  // Normal (non in-place) mode.
  //

  if (nrhs == 2 || (nrhs == 3 && load_wisdom)) { // Normal mode.

    plhs[0] = mxCreateDoubleMatrix(H_M+U_M-1, H_N, mxREAL);
    Y = mxGetPr(plhs[0]);
  }

  //
  // In-place mode.
  //

  if ( (nrhs == 3 && !load_wisdom) || nrhs == 4 ) {

    if ( mxGetM(prhs[2]) != H_M+U_M-1)
      mexErrMsgTxt("Wrong number of rows in argument 4!");

    if ( mxGetN(prhs[2]) != H_N)
      mexErrMsgTxt("Wrong number of columns in argument 4!");

    Y = mxGetPr(prhs[2]);
  }

  //
  // Call the CONV subroutine.
  //

  running = true;

  // Allocate local data.
  D = (DATA*) malloc(nthreads*sizeof(DATA));
  if (!D) {
    mexErrMsgTxt("Failed to allocate memory for thread data!");
  }

  // Allocate mem for the threads.
  threads = new std::thread[nthreads]; // Init thread data.
  if (!threads)
    mexErrMsgTxt("Failed to allocate memory for threads!");

  for (thread_n = 0; thread_n < nthreads; thread_n++) {

    line_start = thread_n * H_N/nthreads;
    line_stop =  (thread_n+1) * H_N/nthreads;

    // Init local data.
    D[thread_n].line_start = line_start; // Local start index;
    D[thread_n].line_stop  = line_stop;  // Local stop index;
    D[thread_n].H = H;
    D[thread_n].H_M = H_M;
    D[thread_n].H_N = H_N;
    D[thread_n].U = U;
    D[thread_n].U_M = U_M;
    D[thread_n].U_N = U_N;
    D[thread_n].Y = Y;
    D[thread_n].L = H_L;

    if (nthreads > 1) {
      // Start the threads.
      threads[thread_n] = std::thread(smp_dream_sum_fftconv, &D[thread_n]);
    } else {
      smp_dream_sum_fftconv(&D[0]);
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

  if (H) {
    free(H);
  }

  //
  // Return the FFTW Wisdom so that the plans can be re-used.
  //

  if (return_wisdom) {
    wisdom_str = fftw_export_wisdom_to_string();
    plhs[1] = mxCreateString(wisdom_str);
    fftw_free(wisdom_str);
  }

  return;
}
