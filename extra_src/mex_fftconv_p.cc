/***
*
* Copyright (C) 2006,2008,2009,2010,2012,2015,2016,2021,2022,2025 Fredrik Lingvall
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

#ifdef DEBUG
#include <mutex>
std::mutex print_mutex;         // Debug print mutex (C++11)
#endif

/***
 *
 *  Parallel (threaded) FFT based convolution.
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
  double *Y;
  FFT *fft;
  ConvMode conv_mode;
} DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//

void* smp_dream_fftconv(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

/***
 *
 * Thread function.
 *
 ***/

void* smp_dream_fftconv(void *arg)
{
  DATA D = *(DATA*) arg;
  dream_idx_type col_start=D.col_start, col_stop=D.col_stop, n;
  double *A = D.A, *B = D.B, *Y = D.Y;
  dream_idx_type A_M = D.A_M, B_M = D.B_M, B_N = D.B_N;
  FFT fft = *D.fft;
  ConvMode conv_mode = D.conv_mode;
  dream_idx_type fft_len = A_M+B_M-1;

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

  if (B_N > 1) {// B is a matrix.

    n=col_start;
    double *Ap = &A[0+n*A_M];
    double *Bp = &B[0+n*B_M];
    double *Yp = &Y[0+n*(A_M+B_M-1)];
    for (n=col_start; n<col_stop; n++) {

      fftconv(fft, Ap, A_M, Bp, B_M, Yp,
              a, b, c, af, bf, cf, conv_mode);

      Ap += A_M;
      Bp += B_M;
      Yp += fft_len;

      if (running==false) {
        std::cout << "fftconv_p: thread for column " << col_start+1 << " -> " << col_stop << " bailing out!\n";
        break;
      }

    } // end-for
  } else { // B is a vector.

    n=col_start;
    double *Ap = &A[0+n*A_M];
    double *Yp = &Y[0+n*(A_M+B_M-1)];
    for (n=col_start; n<col_stop; n++) {

      fftconv(fft, Ap, A_M, B, B_M, Yp, a,b,c,af,bf,cf,conv_mode);

      Ap += A_M;
      Yp += fft_len;

      if (running==false) {
        std::cout << "fftconv_p: thread for column " << col_start+1 << " -> " << col_stop << " bailing out!\n";
        break;
      }

    } // end-for
  } // end-if

#ifdef DEBUG
  {
    std::lock_guard<std::mutex> lk(print_mutex);
    std::cout << "col_start: " << col_start
                  << " col_stop: " << col_stop
                  << " k: " << k
                  << std::endl;
  }
#endif

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
 * Matlab (MEX) gateway function for FFTCONV_P.
 *
 ***/

extern void _main();

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *A=nullptr,*B=nullptr, *Y=nullptr;
  sighandler_t   old_handler, old_handler_abrt, old_handler_keyint;
  std::thread     *threads;
  dream_idx_type thread_n, nthreads;
  dream_idx_type col_start, col_stop, A_M, A_N, B_M, B_N;
  DATA   *D;
  int plan_method = 4; // Default to FFTW_ESTIMATE
  dream_idx_type fft_len;
  bool  return_wisdom = false, load_wisdom = false;
  bool is_set = false;
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
  if ( (nrhs < 2) ||  (nrhs > 5) ) {
    dream_err_msg("fftconv_p requires 2 to 5 input arguments!");
  }

  // Num outputs
  if (nlhs > 2) {
    dream_err_msg("Too many output arguments for fftconv_p!");
  }

  // 2nd output is wisdom string
  if (nlhs == 2) {
    return_wisdom = true;
  }

  A_M = mxGetM(prhs[0]);
  A_N = mxGetN(prhs[0]);
  A = mxGetPr(prhs[0]);

  B_M = mxGetM(prhs[1]);
  B_N = mxGetN(prhs[1]);
  B = mxGetPr(prhs[1]);

  // Check dims of arg 2.
  if ( B_M != 1 && B_N !=1 && B_N != A_N)
    dream_err_msg("Argument 2 must be a vector or a matrix with the same number of rows as arg 1!");

  if (  B_M == 1 || B_N == 1 ) { // B is a vector.
    B_M = B_M*B_N;
    B_N = 1;
  }

  fft_len = A_M+B_M-1;

  // NB We call the FFTW planners only from the main thread once
  // so no need to lock it with a mutex
  FFT fft(fft_len, nullptr, plan_method);

  std::string wisdom_str ="";

  switch (nrhs) {

  case 2:
    break;

  case 3:
    if (mxIsChar(prhs[2])) { // 3rd arg is a FFTW wisdom string.

      char *str = mxArrayToString(prhs[2]);
      wisdom_str += str;
      mxFree(str);

      //
      // If 3rd arg is a string then only a wisdom string is valid.
      //

      if (!fft.is_wisdom(wisdom_str)) {
        dream_err_msg("The string in arg 3 do not seem to be in a FFTW wisdom format!");
      }
      else {
        load_wisdom = true;
      }
    } else { // 3rd arg not a string then assume in-place mode.
      fft.forget_wisdom(); // Clear wisdom history (a new wisdom will be created below).
      if (nlhs > 0) {
        dream_err_msg("3rd arg is not a FFTW wisdom string and in-place mode is assumed. But then there should be no output args!");
      }
    }
    break;

  case 4: // In-place mode if >= 4 args.
    if (mxIsChar(prhs[3])) { // 5th arg is a string (=,+=,-=,or wisdom).

      char *str = mxArrayToString(prhs[3]);
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
    } else { // 4th arg not a string
      dream_err_msg("Argument 4 is not a valid string format!");
    }
    break;

  case 5: // In-place mode if input 5 args.
    if (mxIsChar(prhs[4])) { // 5:th arg is a string (FFTW wisdom).

      char *str = mxArrayToString(prhs[4]);
      wisdom_str += str;
      mxFree(str);

      if (!fft.is_wisdom(wisdom_str)) {
        dream_err_msg("The string in 5th arg do not seem to be in a FFTW wisdom format!");
      } else {
        load_wisdom = true;
      }

    } else { // 5th arg not a string
      dream_err_msg("Argument 5 is not a valid string format!");
    }
    break;

  default:
    dream_err_msg("fftconv_p requires 2 to 5 input arguments!");
    break;
  }

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11)
  nthreads = std::thread::hardware_concurrency();

  // Read DREAM_NUM_THREADS env var
  if(const char* env_p = std::getenv("DREAM_NUM_THREADS")) {
    dream_idx_type dream_threads = std::stoul(env_p);
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
    if (!fft.import_wisdom(wisdom_str)) {
      dream_err_msg("Failed to load FFTW wisdom!");
    }
  }

  //
  // Normal (non in-place) mode.
  //

  if (nrhs == 2 ||  (nrhs == 3 && load_wisdom)) {

    plhs[0] = mxCreateDoubleMatrix(A_M+B_M-1, A_N, mxREAL);
    Y = mxGetPr(plhs[0]);

    SIRData ymat(Y, A_M+B_M-1, A_N);
    ymat.clear();
  }

  //
  // In-place mode.
  //

  if ( (nrhs == 3 && !load_wisdom) || nrhs == 4 || nrhs == 5) {

    if ( mxGetM(prhs[2]) != A_M+B_M-1) {
      dream_err_msg("Wrong number of rows in argument 3!");
    }

    if ( mxGetN(prhs[2]) != A_N) {
      dream_err_msg("Wrong number of columns in argument 3!");
    }

    Y =  mxGetPr(prhs[2]);
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
  if (!threads) {
    dream_err_msg("Failed to allocate memory for threads!");
  }

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
    D[thread_n].fft = &fft;
    D[thread_n].conv_mode = conv_mode;

#ifdef DEBUG
    {
      std::cout << "Init thread_n: "  << thread_n
                << " col_start: " << col_start
                << " col_stop: " << col_stop
                << " nthreads: " << nthreads
                << std::endl;
    }
#endif
    if (nthreads > 1) {
      // Start the threads.
      threads[thread_n] = std::thread(smp_dream_fftconv, &D[thread_n]);
      set_dream_thread_affinity(thread_n, nthreads, threads);
    } else {
      smp_dream_fftconv(&D[0]);
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

  // Return the FFTW Wisdom so that the plans can be re-used.
  if (return_wisdom) {
    std::string cmout = fft.get_wisdom();
    plhs[1] = mxCreateString(cmout.c_str());
  }

  return;
}
