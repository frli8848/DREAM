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

#include "dream.h"
#include "affinity.h"
#include "fftconv.h"

#include "arg_parser_julia.h"

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
 * Julia
 *
 ***/

namespace jl = jlcxx;

jl::ArrayRef<double, 2> jl_fftconv_p(jl::ArrayRef<double, 2> jl_A,
                                  jl::ArrayRef<double, 2> jl_B);

jl::ArrayRef<double, 2> jl_fftconv_p(jl::ArrayRef<double, 2> jl_A,
                                  jl::ArrayRef<double, 2> jl_B)
{
  ArgParser<double> ap;
  sighandler_t old_handler, old_handler_abrt, old_handler_keyint;
  DATA *D = nullptr;
  ConvMode conv_mode = ConvMode::equ;

  //
  // Input args
  //

  double *A = static_cast<double*>(ap.get_data(jl_A));
  dream_idx_type A_M = (dream_idx_type) ap.get_m(jl_A);
  dream_idx_type A_N = (dream_idx_type) ap.get_n(jl_A);

  double *B = static_cast<double*>(ap.get_data(jl_B));
  dream_idx_type B_M = (dream_idx_type) ap.get_m(jl_B);
  dream_idx_type B_N = (dream_idx_type) ap.get_n(jl_B);

  // Check that arg 2 has the correct dim (matrix or vector allowed).
  if ( B_M != 1 && B_N !=1 && B_N != A_N) {
    throw std::runtime_error("Argument 2 must be a vector or a matrix with the same number of columns as arg 1!");
  }

  if (  B_M == 1 || B_N == 1 ) { // B is a vector.
    B_M = B_M*B_N;
    B_N = 1;
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

  // nthreads can't be larger then the number of columns in the A matrix.
  if (nthreads > A_N) {
    nthreads = A_N;
  }

  //
  // FFTW Setup
  //

  int plan_method = 4; // Default to FFTW_ESTIMATE
  dream_idx_type fft_len = A_M+B_M-1;
  //bool return_wisdom = false, load_wisdom = false;

  FFT fft(fft_len, nullptr, plan_method);

  std::string wisdom_str ="";

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
  // Create an output matrix
  //

  jl_value_t *matrix_type = jl_apply_array_type((jl_value_t*) jl_float64_type,2);
  jl_array_t *jl_Y_array = jl_alloc_array_2d(matrix_type, A_M+B_M-1, A_N);
  double *Y = (double *) jl_array_data(jl_Y_array);
  auto jl_Y_mat = jl::ArrayRef<double, 2>(jl_Y_array);

  //
  // Call the FFTCONV subroutine.
  //

  running = true;

  // Allocate local data.
  D = (DATA*) malloc(nthreads*sizeof(DATA));
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
      threads[thread_n] = std::thread(smp_dream_fftconv, &D[thread_n]);
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

  return jl_Y_mat;
}

/***
 *
 *  Julia gateway function for (parallel) fftconv.
 *
 ***/

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.method("fftconv_p", jl_fftconv_p);
}
