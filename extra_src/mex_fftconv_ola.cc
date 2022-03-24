/***
*
* Copyright (C) 2010,2011,2012,2014,2015,2016,2021,2022 Fredrik Lingvall
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
#include <complex>

#include "mex.h"

#include "dream.h"
#include "affinity.h"
#include "fftconv.h"

/***
 *
 *  Parallel (threaded) FFTW overlap-and-add based convolution.
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
  dream_idx_type col_start;
  dream_idx_type col_stop;
  double *A;
  dream_idx_type A_M;
  dream_idx_type A_N;
  double *B;
  dream_idx_type B_M;
  dream_idx_type B_N;
  dream_idx_type block_len;
  double *Y;
  FFT *fft;
  ConvMode conv_mode;
} DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//

void* smp_dream_fftconv_ola(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

/***
 *
 * Thread function
 *
 ***/

void* smp_dream_fftconv_ola(void *arg)
{
  DATA D = *(DATA *)arg;
  dream_idx_type col_start=D.col_start, col_stop=D.col_stop;
  double *A = D.A, *B = D.B;
  double *Y = D.Y;
  FFT fft = *D.fft;
  dream_idx_type A_M = D.A_M, B_M = D.B_M, B_N = D.B_N;
  dream_idx_type block_len = D.block_len;
  ConvMode conv_mode = D.conv_mode;
  dream_idx_type fft_len = block_len + B_M - 1;

  // Input vectors.
  FFTVec a_v(fft_len);
  FFTVec b_v(fft_len);
  FFTVec c_v(fft_len);
  double *a = a_v.get(), *b = b_v.get(), *c = c_v.get();

  FFTVec y_v(fft_len);
  double *y = y_v.get();

  // Fourier Coefficients.
  FFTCVec af_v(fft_len);
  FFTCVec bf_v(fft_len);
  FFTCVec cf_v(fft_len);
  std::complex<double> *af = af_v.get(), *bf = bf_v.get(), *cf = cf_v.get();

  FFTVec a_block_v(fft_len);
  double *a_block = a_block_v.get();

  //
  // Do the convolution.
  //

  if (B_N > 1) {// B is a matrix.

    dream_idx_type n=col_start;
    double *Ap = &A[0+n*A_M];
    double *Bp = &B[0+n*B_M];

    for (n=col_start; n<col_stop; n++) {

      //
      // The overlap-and-add algorithm.
      //

      dream_idx_type k = 0;
      dream_idx_type len = 0;
      while (k < A_M ) {

        //
        // Convolve one block
        //

        if (k+block_len < A_M) {

          fftconv(fft, Ap, block_len, Bp, B_M, y,
                  a, b, c, af, bf, cf, conv_mode);

          len = block_len;

        } else {

          // We must do the convolution (eg. the FFT) of the same lenght here (=fft_len) but the first
          // arg is now shorter than above so copy and zero-pad the first arg to fftconv.

          len = A_M-k;

          // Copy.
          memcpy(a_block, &A[k + n*A_M], len*sizeof(double));

          // Zero-pad.
          memset(&a_block[A_M-k], 0, (block_len - len)*sizeof(double));

          fftconv(fft, a_block, block_len, Bp, B_M, y,
                  a, b, c, af, bf, cf, conv_mode);
        }

        Ap += len;

        //
        // Add the overlap.
        //

        if (k+fft_len <  A_M+B_M-1) {
          for (dream_idx_type l=0; l<fft_len; l++) {
            Y[k+l + n*(A_M+B_M-1)] += y[l];
          }
        } else {
          for (dream_idx_type l=0; l<(A_M+B_M-1)-k; l++) {
            Y[k+l + n*(A_M+B_M-1)] += y[l];
          }
        }

        k += block_len;
      }

      if (running==false) {
        std::cout << "fftconv_ola: thread for column " << col_start+1 << " -> " <<
          col_stop << " bailing out at column = " << n <<"!\n";
        break;
      }

      Bp += B_M;

    } // end-for

  } else { // B is a vector.

    dream_idx_type n=col_start;
    double *Ap = &A[0+n*A_M];

    for (n=col_start; n<col_stop; n++) {

      //
      // The overlap-and-add algorithm.
      //

      dream_idx_type k = 0;
      dream_idx_type len = 0;
      while (k < A_M ) {

        // Convolve one block
        if (k+block_len < A_M) {

          fftconv(fft, Ap, block_len, B, B_M, y,
                  a, b, c, af, bf, cf, conv_mode);

          len = block_len;

        } else {

          // We must do the convolution (eg. the FFT) of the same lenght here (=fft_len) but the first
          // arg is now shorter than above so copy and zero-pad the first arg to fftconv.

          len = A_M-k;

          // Copy.
          memcpy(a_block, &A[k + n*A_M], len*sizeof(double));

          // Zero-pad.
          memset(&a_block[A_M-k], 0, (block_len - len)*sizeof(double));

          fftconv(fft, a_block, block_len, B, B_M, y,
                  a, b, c, af, bf, cf, conv_mode);

        }

        Ap += len;

        // Add the overlap.
        if (k+fft_len < A_M+B_M-1) {
          for (dream_idx_type l=0; l<fft_len; l++) {
            Y[k+l + n*(A_M+B_M-1)] += y[l];
          }
        } else {
          for (dream_idx_type l=0; l<(A_M+B_M-1)-k; l++) {
            Y[k+l + n*(A_M+B_M-1)] += y[l];
          }
        }

        k += block_len;
      } // while

      if (running==false) {
        std::cout << "fftconv_ola: thread for column " << col_start+1 << " -> " <<
          col_stop << " bailing out at column = " << n <<"!\n";
        break;
      }

      Ap += A_M;

    } // for
  } // if

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
 * Matlab (mex) gateway function for FFTCONV_OLA.
 *
 ***/

extern void _main();

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  sighandler_t old_handler, old_handler_abrt, old_handler_keyint;
  std::thread *threads;
  double *Y = nullptr;
  DATA *D = nullptr;
  int plan_method = 4; // Default to FFTW_ESTIMATE
  dream_idx_type fft_len;
  bool return_wisdom = false, load_wisdom = false;
  bool is_set = false;
  dream_idx_type block_len;
  ConvMode conv_mode=ConvMode::equ;

  //
  // Set the method which fftw computes plans
  //

  // If we want to save a plan (in the second output arg)
  // then use the more time-consuming MEAUSURE method.
  if (nlhs == 2) {
    plan_method = 3; // 3 = MEASURE. This takes too long on long fft:s
    //plan_method = 4; // 4 = ESTIMATE.
  }

  //
  // Check for proper inputs arguments.
  //

  // Num inputs
  if ( (nrhs < 3) ||  (nrhs > 6) ) {
    dream_err_msg("fftconv_ola requires 3 to 6 input arguments!");
  }

  // Num outputs
  if (nlhs > 2) {
    dream_err_msg("Too many output arguments for fftconv_ola!");
  }

  if (nlhs == 2) {
    return_wisdom = true;
  }

  dream_idx_type A_M = mxGetM(prhs[0]);
  dream_idx_type A_N = mxGetN(prhs[0]);
  double *A = mxGetPr(prhs[0]);

  dream_idx_type B_M = mxGetM(prhs[1]);
  dream_idx_type B_N = mxGetN(prhs[1]);
  double *B = mxGetPr(prhs[1]);

  // Check dims of arg 2.
  if ( B_M != 1 && B_N !=1 && B_N != A_N) {
    dream_err_msg("Argument 2 must be a vector or a matrix with the same number of rows as arg 1!");
  }

  if ( B_M == 1 || B_N == 1 ) { // B is a vector.
    B_M = B_M*B_N;
    B_N = 1;
  }

  // The segment length (block length).

  if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 1) {
    dream_err_msg("Argument 3 must be a scalar!");
  }
  block_len = (dream_idx_type) mxGetPr(prhs[2])[0]; // Get scalar

  if (block_len < 0 || block_len > A_M) {
    dream_err_msg("Argument 3 is out of bounds! Must be > 0 and less than number of rows of arg 1!");
  }

  fft_len = block_len+B_M-1;

  // NB We call the FFTW planners only from the main thread once
  // so no need to lock it with a mutex
  FFT fft(fft_len, nullptr, plan_method);

  std::string wisdom_str ="";

  // Check for proper inputs arguments.
  switch (nrhs) {

  case 3:
    break;

  case 4: // Output arg in in-place mode or a wisdom string in normal mode.
    {
      if ( mxIsChar(prhs[3]) ) { // 4th arg is a fftw wisdom string.

        char *str = mxArrayToString(prhs[3]);
        wisdom_str += str;
        mxFree(str);

        //
        // If 4:th arg is a string then only a wisdom string is valid.
        //

        if (!fft.is_wisdom(wisdom_str)) {
          dream_err_msg("The string in arg 4 do not seem to be in a fftw wisdom format!");
        }
        else {
          load_wisdom = true;
        }

      } else { // 4th arg not a string then assume in-place mode.
        fft.forget_wisdom(); // Clear wisdom history (a new wisdom will be created below).
        if (nlhs > 0) {
          dream_err_msg("4th arg is not a fftw wisdom string and in-place mode is assumed. But then there should be no output args!");
        }
      }
    }
    break;

  case 5:  // In-place mode if >= 5 args.
    if ( mxIsChar(prhs[4]) ) { // 5th arg is a string ('=','+=','-=', or fftw wisdom).

      char *str = mxArrayToString(prhs[4]);
      std::string ip_mode(str);
      mxFree(str);

      // Valid strings are:
      //  '='  : In-place replace mode.
      //  '+=' : In-place add mode.
      //  '-=' : In-place sub mode.
      //  wisdom string.

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
        if (!fft.is_wisdom(ip_mode)) {
          dream_err_msg("Non-valid string in arg 4!");
        } else {
          wisdom_str = ip_mode;
          load_wisdom = true;
        }
      }
    } else { // 5th arg not a string
      dream_err_msg("Argument 5 is not a valid string format!");
    }
    break;

  case 6: // In-place mode if 6 input args.
    if ( mxIsChar(prhs[5]) ) { // 6:th arg is a string (fftw wisdom).

      // Read the wisdom string.
      char *str = mxArrayToString(prhs[5]);
      wisdom_str += str;
      mxFree(str);

      if (!fft.is_wisdom(wisdom_str) ) {
        dream_err_msg("The string in 5th arg do not seem to be in a FFTW wisdom format!");
      } else {
        load_wisdom = true;
      }
    } else { // 6th arg not a string
      dream_err_msg("Argument 6 is not a valid string format!");
    }
    break;

  default:
    {
      dream_err_msg("fftconv_ola requires 3 to 6 input arguments!");
    }
    break;
  }

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11)
  dream_idx_type nthreads = std::thread::hardware_concurrency();

  // Read DREAM_NUM_THREADS env var
  if(const char* env_p = std::getenv("DREAM_NUM_THREADS")) {
    dream_idx_type dream_threads = std::stoul(env_p);
    if (dream_threads < nthreads) {
      nthreads = dream_threads;
    }
  }

  // Check nthreads argument.
  if (nthreads > A_N) {
    // Add a -v verbose flag to display stuff like this!
    //mexPrintf("Warning: nthreads is larger then number of columns in first arg.\n");
    //mexPrintf("         Setting nthreads = # cols in 1st arg!\n");
    nthreads = A_N;
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

  // ********************************************
  //
  // Init the FFTW plans.
  //
  // ********************************************

  if(load_wisdom) {
    if (!fft.import_wisdom(wisdom_str)) {
      dream_err_msg("Failed to load FFTW wisdom!");
    }
  }

  //
  // Normal (non in-place) mode.
  //

  if (nrhs == 3 || (nrhs == 4 && load_wisdom)) { // Normal mode.

    // Allocate mamory for the output matrix.
    plhs[0] = mxCreateDoubleMatrix(A_M+B_M-1, A_N, mxREAL);
    Y = mxGetPr(plhs[0]);

    SIRData ymat(Y, A_M+B_M-1, A_N);
    ymat.clear();
  }

  //
  // In-place mode.
  //

  if ( (nrhs == 4 && !load_wisdom) || nrhs == 5 || nrhs == 6) {

    if ( mxGetM(prhs[3]) != A_M+B_M-1) {
      dream_err_msg("Wrong number of rows in argument 3!");
    }

    if ( mxGetN(prhs[3]) != A_N) {
      dream_err_msg("Wrong number of columns in argument 3!");
    }

    Y = mxGetPr(prhs[3]);
  }

  // Init thread data.
  threads = new std::thread[nthreads];
  if (!threads) {
    dream_err_msg("Failed to allocate memory for threads!");
  }

  // Allocate thread data.
  D = (DATA*) malloc(nthreads*sizeof(DATA));
  if (!D) {
    dream_err_msg("Failed to allocate memory for thread data!");
  }

  running = true;

  for (dream_idx_type thread_n=0; thread_n<nthreads; thread_n++) {

    dream_idx_type col_start = thread_n * A_N/nthreads;
    dream_idx_type col_stop = (thread_n+1) * A_N/nthreads;

    // Init local data.
    D[thread_n].col_start = col_start; // Local start index;
    D[thread_n].col_stop = col_stop; // Local stop index;
    D[thread_n].A = A;
    D[thread_n].A_M = A_M;
    D[thread_n].A_N = A_N;
    D[thread_n].B = B;
    D[thread_n].B_M = B_M;
    D[thread_n].B_N = B_N;
    D[thread_n].block_len = block_len;
    D[thread_n].Y = Y;
    D[thread_n].fft = &fft;
    D[thread_n].conv_mode = conv_mode;

    if (nthreads > 1) {
      // Start the threads.
      threads[thread_n] = std::thread(smp_dream_fftconv_ola, &D[thread_n]);
    } else {
      smp_dream_fftconv_ola(&D[0]);
    }

  } // for (thread_n = 0; thread_n < nthreads; thread_n++)

    //
    // Wait for all threads to finish.
    //

  if (nthreads > 1) {
    for (dream_idx_type thread_n=0; thread_n < nthreads; thread_n++) {
      threads[thread_n].join();
    }
  }

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

  // Return the FFTW Wisdom so that the plans can be re-used.
  if (return_wisdom) {
    std::string cmout = fft.get_wisdom();
    plhs[1] = mxCreateString(cmout.c_str());
  }

  return;
}
