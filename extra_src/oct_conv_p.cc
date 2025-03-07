/***
*
* copyright (C) 2006,2007,2008,2009,2014,2015,2016,2022,2023,2025 Fredrik Lingvall
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

#include <octave/oct.h>

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
  octave_idx_type col_start;
  octave_idx_type col_stop;
  double *A;
  octave_idx_type A_M;
  octave_idx_type A_N;
  double *B;
  octave_idx_type B_M;
  octave_idx_type B_N;
  double *Y;
  ConvMode conv_mode;
} DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//

void* smp_dream_conv_p(void *arg);
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
  octave_idx_type    col_start=D.col_start, col_stop=D.col_stop;
  double *A = D.A, *B = D.B, *Y = D.Y;
  octave_idx_type A_M = D.A_M, B_M = D.B_M, B_N = D.B_N;
  ConvMode conv_mode = D.conv_mode;

  // Do the convolution.

  for (dream_idx_type n=col_start; n<col_stop; n++) {

    if (B_N > 1) // B is a matrix.
      conv( &A[0+n*A_M], A_M, &B[0+n*B_M], B_M, &Y[0+n*(A_M+B_M-1)], conv_mode);
    else // B is a vector.
      conv( &A[0+n*A_M], A_M, B, B_M, &Y[0+n*(A_M+B_M-1)], conv_mode);

    if (running == false) {
      octave_stdout << "conv_p: thread for column " << col_start+1 << " -> " << col_stop << "bailing out!\n";
      break;
      //return(NULL);
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
 * Octave (oct) gateway function for CONV_P.
 *
 ***/

DEFUN_DLD (conv_p, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {}  [Y] = conv_p(A, B).\n\
\n\
CONV_P Computes one dimensional convolutions of the columns in the matrix A and the matrix (or vector) B.\n\
\n\
Input parameters:\n\
\n\
@table @samp\n\
@item A\n\
An MxN matrix.\n\
@item B\n\
A KxN matrix or a K-length vector. If B is a vector each column in A is convolved with the vector B.\n\
\n\
Output parameter:\n\
\n\
@table @samp\n\
@item Y\n\
The (M+K-1)xN output matrix.\n\
@end table\n\
\n\
conv_p is an oct-function that is a part of the DREAM Toolbox available at\n\
 @url{https://github.com/frli8848/DREAM}.\n\
\n\
Copyright @copyright{} 2006-2023 Fredrik Lingvall.\n\
@seealso {fftconv_p, fftconv, conv}\n\
@end deftypefn")
{
  double *A = nullptr, *B = nullptr, *Y = nullptr;
  sighandler_t old_handler, old_handler_abrt, old_handler_keyint;
  octave_idx_type col_start, col_stop, A_M, A_N, B_M, B_N;
  bool is_set = false;
  DATA *D = nullptr;
  ConvMode conv_mode = ConvMode::equ;

  octave_value_list oct_retval;

  int nrhs = args.length ();

  // Check for proper inputs arguments.

  switch (nrhs) {

  case 0:
  case 1:
    dream_err_msg("conv_p requires 2 to 4 input arguments!");
    return oct_retval;
    break;

  case 2:
    if (nlhs > 1) {
      dream_err_msg("Too many output arguments for conv_p!");
      return oct_retval;
    }
    break;

  case 3:
    if (nlhs > 0) {
      dream_err_msg("No output arguments required for conv_p in in-place operating mode!");
      return oct_retval;
    }
    break;

  case 4:
    if ( args(3).is_string() ) { // 4:th arg is a string '=', '+=', or '-='.

      std::string ip_mode = args(3).string_value();

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
        return oct_retval;
      }
    }
    break;

  default:
    dream_err_msg("conv_p requires 2 to 4 input arguments!");
    return oct_retval;
    break;
  }

  const Matrix tmp0 = args(0).matrix_value();
  A_M = tmp0.rows();
  A_N = tmp0.cols();
  A = (double*) tmp0.data();

  const Matrix tmp1 = args(1).matrix_value();
  B_M = tmp1.rows();
  B_N = tmp1.cols();
  B = (double*) tmp1.data();

  // Check that arg 2 has the correct dim (matrix or vector allowed).
  if ( B_M != 1 && B_N !=1 && B_N != A_N) {
    dream_err_msg("Argument 2 must be a vector or a matrix with the same number of columns as arg 1!");
    return oct_retval;
  }

  if (  B_M == 1 || B_N == 1 ) { // B is a vector.
    B_M = B_M*B_N;
    B_N = 1;
  }

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11)
  octave_idx_type nthreads = std::thread::hardware_concurrency();

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

  if (nrhs == 2) { // Test for in-place/normal mode.

    //
    // Normal (non in-place) mode.
    //

    Matrix Ymat(A_M+B_M-1, A_N);
    Y = (double*) Ymat.data();

    oct_retval.append(Ymat);

  } else {

    //
    // In-place mode.
    //

    if (args(2).matrix_value().rows() != A_M+B_M-1) {
      dream_err_msg("Wrong number of rows in argument 3!");
      return oct_retval;
    }

    if (args(2).matrix_value().cols() != A_N) {
      dream_err_msg("Wrong number of columns in argument 3!");
      return oct_retval;
    }

    const Matrix Ytmp = args(2).matrix_value();

    Y = (double*) Ytmp.data();
  }

  //
  // Call the CONV subroutine.
  //

  running = true;

  // Allocate local data.
  D = (DATA*) malloc(nthreads*sizeof(DATA));
  if (!D) {
    dream_err_msg("Failed to allocate memory for thread data!");
    return oct_retval;
  }

  // Allocate mem for the threads.
  std::thread *threads = new std::thread[nthreads]; // Init thread data.
  if (!threads) {
    dream_err_msg("Failed to allocate memory for threads!");
    return oct_retval;
  }

  for (octave_idx_type thread_n = 0; thread_n < nthreads; thread_n++) {

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
    for (octave_idx_type thread_n = 0; thread_n < nthreads; thread_n++) {
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
    return oct_retval;
  }

  return oct_retval;
}
