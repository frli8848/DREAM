/***
*
* Copyright (C) 2024 Fredrik Lingvall
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
#include <atomic>

#include "affinity.h"
#include "fftconv.h"
#include "fftconv_p.h"

std::mutex err_mutex_rect;
std::atomic<bool> running_conv;

void FFTConvP::abort(int signum)
{
  running_conv = false;
}

bool FFTConvP::is_running()
{
  return running_conv;
}

/***
 *
 *  Parallel (threaded) FFT based convolution.
 *
 ***/

//
// typedef:s
//

typedef struct
{
  dream_idx_type col_start;
  dream_idx_type col_stop;
  const double *A;
  dream_idx_type A_M;
  dream_idx_type A_N;
  const double *B;
  dream_idx_type B_M;
  dream_idx_type B_N;
  double *Y;
  FFT *fft;
  ConvMode conv_mode;
} DATA;

/***
 *
 * Thread function.
 *
 ***/

void* FFTConvP::smp_dream_fftconv(void *arg)
{
  DATA D = *(DATA*) arg;
  dream_idx_type col_start=D.col_start, col_stop=D.col_stop, n;
  const double *A = D.A, *B = D.B;
  double *Y = D.Y;
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
    const double *Ap = &A[0+n*A_M];
    const double *Bp = &B[0+n*B_M];
    double *Yp = &Y[0+n*(A_M+B_M-1)];
    for (n=col_start; n<col_stop; n++) {

      fftconv(fft, Ap, A_M, Bp, B_M, Yp,
              a, b, c, af, bf, cf, conv_mode);

      Ap += A_M;
      Bp += B_M;
      Yp += fft_len;

      if (running_conv == false) {
        std::cout << "fftconv_p: thread for column " << col_start+1 << " -> " << col_stop << " bailing out!\n";
        break;
      }

    } // end-for
  } else { // B is a vector.

    n=col_start;
    const double *Ap = &A[0+n*A_M];
    double *Yp = &Y[0+n*(A_M+B_M-1)];
    for (n=col_start; n<col_stop; n++) {

      fftconv(fft, Ap, A_M, B, B_M, Yp, a,b,c,af,bf,cf,conv_mode);

      Ap += A_M;
      Yp += fft_len;

      if (running_conv == false) {
        std::cout << "fftconv_p: thread for column " << col_start+1 << " -> " << col_stop << " bailing out!\n";
        break;
      }

    } // end-for
  } // end-if

  return(NULL);
}

/***
 *
 * fftconv_p
 *
 ***/

ErrorLevel FFTConvP::run(double *Y,
                         const double *A, dream_idx_type A_M, dream_idx_type A_N,
                         const double *B, dream_idx_type B_M, dream_idx_type B_N,
                         ConvMode mode)
{
  int plan_method = 4; // Default to FFTW_ESTIMATE
  bool return_wisdom = false, load_wisdom = false;
  ConvMode conv_mode = mode;
  ErrorLevel err = ErrorLevel::none;


  // If we want to save a plan (in the second output arg)
  // then use the more time-consuming MEAUSURE method.
  //if (nlhs == 2) {
  plan_method = 3; // 3 = MEASURE. This takes too long on long fft:s
  //plan_method = 4; // 4 = ESTIMATE.
  //}

  //
  // Check for proper inputs arguments.
  //

  //// 2nd output is wisdom string
  //if (nlhs == 2) {
  //  return_wisdom = true;
  //}


  // Check that arg 2 has the correct dim (matrix or vector allowed).
  //if ( B_M != 1 && B_N !=1 && B_N != A_N) {
  //  error("Argument 2 must be a vector or a matrix with the same number of columns as arg 1!");
  //}

  if (  B_M == 1 || B_N == 1 ) { // B is a vector.
    B_M = B_M*B_N;
    B_N = 1;
  }

  dream_idx_type fft_len = A_M+B_M-1;

  // NB We call the FFTW planners only from the main thread once
  // so no need to lock it with a mutex (FFTW planners are not thread
  // safe).
  FFT fft(fft_len, nullptr, plan_method);

  std::string wisdom_str ="";

  /* FIXME: Add support for reading/returning FFTW plans and set conv mode.

  switch (nrhs) {

  case 2:
    break;

  case 3:
    if ( args(2).is_string() ) { // 3rd arg is a FFTW wisdom string.

      wisdom_str = args(2).string_value();

      //
      // If 3rd arg is a string then only a wisdom string is valid.
      //

      if (!fft.is_wisdom(wisdom_str)) {
        error("The string in arg 3 do not seem to be in a FFTW wisdom format!");
        return oct_retval;
      }
      else {
        load_wisdom = true;
      }
    } else { // 3rd arg not a string then assume in-place mode.
      fft.forget_wisdom(); // Clear wisdom history (a new wisdom will be created below).
      if (nlhs > 0) {
        error("3rd arg is not a FFTW wisdom string and in-place mode is assumed. But then there should be no output args!");
        return oct_retval;
      }
    }
    break;

  case 4:  // In-place mode if >= 4 args.
    if ( args(3).is_string() ) { // 4th arg is a string ('=','+=','-=', or wisdom).

      std::string ip_mode = args(3).string_value();

      // Valid strings are:
      //  '='  : In-place replace mode.
      //  '+=' : In-place add mode.
      //  '-=' : In-place sub mode.
      //  wisdom string.

      bool is_set = false;

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
          error("Non-valid string in arg 4!");
          return oct_retval;
        } else {
          wisdom_str = ip_mode;
          load_wisdom = true;
        }
      }
    } else { // 4th arg not a string
      error("Argument 4 is not a valid string format!");
      return oct_retval;
    }
    break;

  case 5: // In-place mode if input 5 args.
    if ( args(4).is_string() ) { // 5:th arg is a string (FFTW wisdom).

      // Read the wisdom string.
      wisdom_str = args(4).string_value();

      if (!fft.is_wisdom(wisdom_str)) {
        error("The string in 5th arg do not seem to be in a FFTW wisdom format!");
        return oct_retval;
      } else {
        load_wisdom = true;
      }

    } else { // 5th arg not a string
      error("Argument 5 is not a valid string format!");
      return oct_retval;
    }
    break;

  default:
    error("fftconv_p requires 2 to 5 input arguments!");
    return oct_retval;
    break;
  }
  */

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11)
  dream_idx_type nthreads = std::thread::hardware_concurrency();

  // Read DREAM_NUM_THREADS env var
  if (const char* env_p = std::getenv("DREAM_NUM_THREADS")) {
    dream_idx_type dream_threads = std::stoul(env_p);
    if (dream_threads < nthreads) {
      nthreads = dream_threads;
    }
  }

  // nthreads can't be larger then the number of columns in the A matrix.
  if (nthreads > A_N) {
    nthreads = A_N;
  }

  //
  // Init the FFTW plans.
  //

  /*
  if(load_wisdom) {
    if (!fft.import_wisdom(wisdom_str)) {
      error("Failed to load FFTW wisdom!");
    }
  }
  */

  //
  // Call the CONV subroutine.
  //

  running_conv = true;

  // Allocate local data.
  DATA   *D = (DATA*) malloc(nthreads*sizeof(DATA));
  if (!D) {
    throw std::runtime_error("Failed to allocate memory for thread data!");
  }

  // Allocate mem for the threads.
  std::thread *threads = new std::thread[nthreads]; // Init thread data.
  if (!threads) {
    throw std::runtime_error("Failed to allocate memory for threads!");
  }

  for (dream_idx_type thread_n = 0; thread_n < nthreads; thread_n++) {

    dream_idx_type col_start = thread_n * A_N/nthreads;
    dream_idx_type col_stop =  (thread_n+1) * A_N/nthreads;

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

    if (nthreads > 1) {
      // Start the threads.
      threads[thread_n] =  conv_thread(&D[thread_n]);
      set_dream_thread_affinity(thread_n, nthreads, threads);
    } else {
      smp_dream_fftconv(&D[0]);
    }
  }

  if (nthreads > 1) {
    // Wait for all threads to finish.
    for (dream_idx_type thread_n = 0; thread_n < nthreads; thread_n++) {
      threads[thread_n].join();
    }
  }

  // Free memory.
  if (D) {
    free((void*) D);
  }

  // FIXME: Return the FFTW Wisdom so that the plans can be re-used.
  /*
  if (return_wisdom) {
    std::string cmout = fft.get_wisdom();
    oct_retval.append(cmout); // Add to output args.
  }
  */

  return err;
}
