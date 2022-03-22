/***
*
* Copyright (C) 2006,2007,2008,2009,2010,2012,2014,2016,2021,2022 Fredrik Lingvall
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
#include <complex> // C++

//
// Octave headers.
//

#include <octave/oct.h>

//
// Macros
//

#ifdef DEBUG
#include <mutex>
std::mutex print_mutex;         // Debug print mutex (C++11)
#endif

#include "dream.h"
#include "affinity.h"
#include "fftconv.h"

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
  octave_idx_type col_start;
  octave_idx_type col_stop;
  double *A;
  octave_idx_type A_M;
  octave_idx_type A_N;
  double *B;
  octave_idx_type B_M;
  octave_idx_type B_N;
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
  octave_idx_type col_start=D.col_start, col_stop=D.col_stop, n;
  double *A = D.A, *B = D.B, *Y = D.Y;
  octave_idx_type A_M = D.A_M, B_M = D.B_M, B_N = D.B_N;
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
        octave_stdout << "fftconv_p: thread for column " << col_start+1 << " -> " << col_stop << " bailing out!\n";
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
        octave_stdout << "fftconv_p: thread for column " << col_start+1 << " -> " << col_stop << " bailing out!\n";
        break;
      }

    } // end-for
  } // end-if

#ifdef DEBUG
  {
    std::lock_guard<std::mutex> lk(print_mutex);
    octave_stdout << "col_start: " << col_start
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
 * Octave (oct) gateway function for FFTCONV_P.
 *
 ***/

DEFUN_DLD (fftconv_p, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} [Y,wisdom_str_out] = fftconv_p(A,B,wisdom_str_in);\n\
\n\
FFTCONV_P - Computes the one dimensional convolution of the columns in matrix A and the matrix (or vector) B.\n\
\n\
Input parameters:\n\
\n\
@table @code\n\
@item A\n\
An M x N matrix.\n\
@item B\n\
A K x N matrix or a K-length vector. If B is a vector each column in A is convolved with the vector B.\n\
@item wisdom_str_in\n\
Optional parameter. If the wisdom_str_in parameter is not supplied then @code{fftconv_p} calls FFTW wisdom plan\n\
functions before performing any frequency domain operations. This overhead can be avoided by supplying\n\
a pre-computed FFTW wisdom string wisdom_str_in. For more information see the FFTW user manunal\n\
available at @url{http://www.fftw.org}.\n\
@end table\n\
\n\
Output parameters:\n\
\n\
@table @code\n\
@item Y\n\
The (M+K-1) x N output matrix.\n\
@item wisdom_str_out\n\
Optional parameter. If the wisdom_str_out output parameter is supplied then @code{fftconv_p}\n\
will call FFTW wisdom plan functions and return the wisdom string which then can be used to\n\
speed up subsequent calls to @code{fftconv_p} by suppying the string as the input argument\n\
wisdom_str_in.\n\
@end tabl\n\
\n\
In place modes:\n\
\n\
fftconv_p(A,B,Y);\n\
\n\
fftconv_p(A,B,Y,mode);\n\
\n\
fftconv_p(A,B,Y,wisdom_str_in);\n\
\n\
fftconv_p(A,B,Y,mode,wisdom_str_in);\n\
\n\
where the 'mode' is a string which can be '=', '+=', or '-='.\n\
\n\
fftconv_p is a part of the DREAM Toolbox available at\n\
@url{http://www.signal.uu.se/Toolbox/dream/}.\n\
\n\
Copyright @copyright{} 2006-2022 Fredrik Lingvall.\n\
@seealso {conv, conv_p, fftconv, fftw_wisdom}\n\
@end deftypefn")
{
  double *A,*B, *Y;
  octave_idx_type col_start, col_stop, A_M, A_N, B_M, B_N, n;
  sighandler_t    old_handler, old_handler_abrt, old_handler_keyint;
  std::thread     *threads;
  octave_idx_type  thread_n, nthreads;
  DATA   *D;
  int plan_method = 4; // Default to FFTW_ESTIMATE
  octave_idx_type fft_len, return_wisdom = false, load_wisdom = false;
  int buflen, is_set = false;
  ConvMode conv_mode=ConvMode::equ;
  octave_value_list oct_retval;

  int nrhs = args.length ();

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
  if ( (nrhs < 2) ||  (nrhs < 2) ) {
    error("fftconv_p requires 2 to 5 input arguments!");
    return oct_retval;
  }

  // Num outputs
  if (nlhs > 2) {
    error("Too many output arguments for fftconv_p!");
    return oct_retval;
  }

  // 2nd output is wisdom string
  if (nlhs == 2) {
    return_wisdom = true;
  }

  const Matrix tmp0 = args(0).matrix_value();
  A_M = tmp0.rows();
  A_N = tmp0.cols();
  A = (double*) tmp0.fortran_vec();

  const Matrix tmp1 = args(1).matrix_value();
  B_M = tmp1.rows();
  B_N = tmp1.cols();
  B = (double*) tmp1.fortran_vec();

  // Check that arg 2.
  if ( B_M != 1 && B_N !=1 && B_N != A_N) {
    error("Argument 2 must be a vector or a matrix with the same number of columns as arg 1!");
    return oct_retval;
  }

  if (  B_M == 1 || B_N == 1 ) { // B is a vector.
    B_M = B_M*B_N;
    B_N = 1;
  }

  fft_len = A_M+B_M-1;

  //std::mutex fft_mutex;
  //FFT fft(fft_len, &fft_mutex, plan_method);
  //
  // NB We call the FFTW planners only from the main thread once
  // so no need to lock it with a mutex
  FFT fft(fft_len, nullptr, plan_method);

  std::string wisdom_str ="";

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
        if (fft.is_wisdom(ip_mode) < 0 ) {
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

      if (fft.is_wisdom(wisdom_str) < 0 ) {
        error("The string in 5th arg do not seem to be in a FFTW wisdom format!");
        return oct_retval;
      } else
        load_wisdom = true;

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

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11)
  nthreads = std::thread::hardware_concurrency();

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
      error("Failed to load FFTW wisdom!");
      return oct_retval;
    }
  }

  if (nrhs == 2 || (nrhs == 3 && load_wisdom)) { // Normal mode.

    //
    // Normal (non in-place) mode.
    //

    Matrix Ymat(A_M+B_M-1, A_N);
    Y = Ymat.fortran_vec();

    SIRData ymat(Y, A_M+B_M-1, A_N);
    ymat.clear();

    //
    // Call the CONV subroutine.
    //

    running = true;

    // Allocate local data.
    D = (DATA*) malloc(nthreads*sizeof(DATA));
    if (!D) {
      error("Failed to allocate memory for thread data!");
      return oct_retval;
    }

    // Allocate mem for the threads.
    threads = new std::thread[nthreads]; // Init thread data.
    if (!threads) {
      error("Failed to allocate memory for threads!");
      return oct_retval;
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

#ifdef DEBUG
      {
        octave_stdout << "Init thread_n: "  << thread_n
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
      error("CTRL-C pressed!\n"); // Bail out.
      return oct_retval;
    }

    // 1st output arg.
    oct_retval.append(Ymat);

    // Return the FFTW Wisdom so that the plans can be re-used.
    if (return_wisdom) {
      std::string cmout = fft.get_wisdom();
      oct_retval.append(cmout); // Add to output args.
    }

    return oct_retval;
  }

  if ( (nrhs == 3 && !load_wisdom) || nrhs == 4 || nrhs == 5) { // In-place mode.

    //
    // In-place mode.
    //

    if (  args(2).matrix_value().rows() != A_M+B_M-1) {
      error("Wrong number of rows in argument 3!");
      return oct_retval;
    }

    if (  args(2).matrix_value().cols() != A_N) {
      error("Wrong number of columns in argument 3!");
      return oct_retval;
    }

    const Matrix Ymat = args(2).matrix_value();
    Y = (double*) Ymat.fortran_vec(); // NB. Do  not clear data here!

    //
    // Call the CONV subroutine.
    //

    running = true;

    // Allocate local data.
    D = (DATA*) malloc(nthreads*sizeof(DATA));
    if (!D) {
      error("Failed to allocate memory for thread data!");
      return oct_retval;
    }

    // Allocate mem for the threads.
    threads = new std::thread[nthreads]; // Init thread data.
    if (!threads) {
      error("Failed to allocate memory for threads!");
      return oct_retval;
    }

    for (thread_n = 0; thread_n < nthreads; thread_n++) {

      col_start = thread_n * A_N/nthreads;
      col_stop =  (thread_n+1) * A_N/nthreads;

      // Init local data.
      D[thread_n].col_start = col_start; // Local start index;
      D[thread_n].col_stop = col_stop;   // Local stop index;
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
      } else {
        smp_dream_fftconv(&D[0]);
      }
    } // for (thread_n = 0; thread_n < nthreads; thread_n++)

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
      error("CTRL-C pressed!\n"); // Bail out.
      return oct_retval;
    }

    // Return the FFTW Wisdom so that the plans can be re-used.
    if (return_wisdom) {
      std::string cmout = fft.get_wisdom();
      oct_retval.append(cmout); // Add to output args.
    }

    return oct_retval;
  }

  return oct_retval; // Just to fix compiler warnings.
}
