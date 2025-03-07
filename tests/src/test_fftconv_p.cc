/***
*
* Copyright (C) 2022,2025 Fredrik Lingvall
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
    double *Ap = (double*) &A[0+n*A_M];
    double *Bp = (double*) &B[0+n*B_M];
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
    double *Ap = (double*) &A[0+n*A_M];
    double *Yp = &Y[0+n*(A_M+B_M-1)];
    for (n=col_start; n<col_stop; n++) {

      fftconv(fft, Ap, A_M, (double*) B, B_M, Yp, a,b,c,af,bf,cf,conv_mode);

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
 *
 *
 ***/

int fftconv_p(double *Y, std::string &wisdom_str_out, std::string ip_mode,
              const double *A, dream_idx_type A_M, dream_idx_type A_N,
              const double *B, dream_idx_type B_M, dream_idx_type B_N,
              std::string wisdom_str_in)
{
  sighandler_t    old_handler, old_handler_abrt, old_handler_keyint;
  std::thread     *threads;
  dream_idx_type  thread_n, nthreads;
  dream_idx_type col_start, col_stop;
  DATA   *D;
  int plan_method = 4; // Default to FFTW_ESTIMATE
  bool return_wisdom = false, load_wisdom = false;
  bool is_set = false;
  ConvMode conv_mode=ConvMode::equ;

  dream_idx_type fft_len = A_M+B_M-1;

  //std::mutex fft_mutex;
  //FFT fft(fft_len, &fft_mutex, plan_method);
  //
  // NB We call the FFTW planners only from the main thread once
  // so no need to lock it with a mutex
  FFT fft(fft_len, nullptr, plan_method);

  std::string wisdom_str ="";

  if ( wisdom_str_in != "" ) { // 3rd arg is a FFTW wisdom string.

    wisdom_str = wisdom_str_in;

    //
    // If 3rd arg is a string then only a wisdom string is valid.
    //

    if (!fft.is_wisdom(wisdom_str)) {
      std::cerr << "The string in arg 10 do not seem to be in a FFTW wisdom format!";
      return -1;
    } else {
      load_wisdom = true;
    }
  } else { // 3rd arg not a string then assume in-place mode.
    fft.forget_wisdom(); // Clear wisdom history (a new wisdom will be created below).
  }

  //
  // In-place mode
  //

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

  if (!is_set) {
    std::cerr << "The string in argument 4 is not a valid in-place format!";
    return -1;
  }

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11)
  nthreads = std::thread::hardware_concurrency();

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
      std::cerr << "Failed to load FFTW wisdom!";
      return -1;
    }
  }

  //
  // Call the CONV subroutine.
  //

  running = true;

  // Allocate local data.
  D = (DATA*) malloc(nthreads*sizeof(DATA));
  if (!D) {
    std::cerr << "Failed to allocate memory for thread data!";
    return -1;
  }

  // Allocate mem for the threads.
  threads = new std::thread[nthreads]; // Init thread data.
  if (!threads) {
    std::cerr << "Failed to allocate memory for threads!";
    return -1;
  }

  for (thread_n = 0; thread_n < nthreads; thread_n++) {

    col_start = thread_n * A_N/nthreads;
    col_stop =  (thread_n+1) * A_N/nthreads;

    // Init local data.
    D[thread_n].col_start = col_start; // Local start index;
    D[thread_n].col_stop = col_stop; // Local stop index;
    D[thread_n].A =  A;
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
    std::cerr << "CTRL-C pressed!\n"; // Bail out.
    return 0;
  }

  // Return the FFTW Wisdom so that the plans can be re-used.
  if (return_wisdom) {
    std::string cmout = fft.get_wisdom();
  }

  return 0;
}


int main(int argc, char *argv[])
{

  dream_idx_type A_M = 1000, A_N = 1000;
  dream_idx_type B_M = 1000, B_N = 1;
  SIRData Amat(A_M, A_N);
  SIRData Bmat(B_M, B_N);
  SIRData Ymat(A_M+B_M-1, A_N);

  std::string wisdom_str_out = "";
  std::string wisdom_str_in = "";
  std::string ip_mode = "=";

  int res = fftconv_p(Ymat.get(), wisdom_str_out, ip_mode,
                      Amat.get(), A_M, A_N,
                      Bmat.get(), B_M, B_N,
                      wisdom_str_in);


  return res;
}
