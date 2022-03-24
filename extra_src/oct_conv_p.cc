/***
*
* Copyright (C) 2006,2007,2008,2009,2014,2015,2016 Fredrik Lingvall
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

#include <iostream>
#include <csignal>
#include <thread>

#include "dream.h"
//#include "dream_error.h"

//
// Octave headers.
//

#include <octave/oct.h>

#define EQU 0
#define SUM 1
#define NEG 2

/***
 *
 *  Parallel (threaded) convolution.
 *
 ***/

//
// Globals
//

volatile int running;
volatile int in_place;
int mode = EQU;

//
// typedef:s
//

typedef struct
{
  octave_idx_type line_start;
  octave_idx_type line_stop;
  double *A;
  octave_idx_type A_M;
  octave_idx_type A_N;
  double *B;
  octave_idx_type B_M;
  octave_idx_type B_N;
  double *Y;
} DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//

void* smp_dream_conv_p(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

void conv(double *xr, octave_idx_type nx, double *yr, octave_idx_type ny, double *zr,
           int in_place, int mode);

/***
 *
 * Thread function.
 *
 ***/

void* smp_dream_conv_p(void *arg)
{
  DATA D = *(DATA *)arg;
  octave_idx_type    line_start=D.line_start, line_stop=D.line_stop, n;
  double *A = D.A, *B = D.B, *Y = D.Y;
  octave_idx_type A_M = D.A_M, B_M = D.B_M, B_N = D.B_N;

  // Do the convolution.

  for (n=line_start; n<line_stop; n++) {

    if (B_N > 1) // B is a matrix.
      conv( &A[0+n*A_M], A_M, &B[0+n*B_M], B_M, &Y[0+n*(A_M+B_M-1)],in_place,mode);
    else // B is a vector.
      conv( &A[0+n*A_M], A_M, B, B_M, &Y[0+n*(A_M+B_M-1)],in_place,mode);

    if (running == false) {
      octave_stdout << "conv_p: thread for column " << line_start+1 << " -> " << line_stop << "bailing out!\n";
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
Copyright @copyright{} 2006-2021 Fredrik Lingvall.\n\
@seealso {fftconv_p, fftconv, conv}\n\
@end deftypefn")
{
  double *A,*B, *Y;
  sighandler_t   old_handler, old_handler_abrt, old_handler_keyint;
  octave_idx_type line_start, line_stop, A_M, A_N, B_M, B_N, n;
  char   *the_str = NULL;
  int    buflen, is_set = false;
  DATA   *D;
  std::thread     *threads;
  octave_idx_type  thread_n, nthreads;
  octave_value_list oct_retval;

  in_place = false;

  int nrhs = args.length ();

  // Check for proper inputs arguments.

  switch (nrhs) {

  case 0:
  case 1:
    error("conv_p requires 2 to 4 input arguments!");
    return oct_retval;
    break;

  case 2:
    if (nlhs > 1) {
      error("Too many output arguments for conv_p!");
      return oct_retval;
    }
    break;

  case 3:
    if (nlhs > 0) {
      error("No output arguments required for conv_p in in-place operating mode!");
      return oct_retval;
    }
    break;

  case 4:
    if ( args(3).is_string() ) { // 4:th arg is a string '=', '+=', or '-='.
      std::string strin = args(3).string_value();
      buflen = strin.length();
      the_str = (char*) malloc(buflen * sizeof(char));
      for ( n=0; n<buflen; n++ ) {
        the_str[n] = strin[n];
      }

      // Valid strings are:
      //  '='  : In-place replace mode.
      //  '+=' : In-place add mode.
      //  '-=' : In-place sub mode.

      is_set = false;

      if (strcmp(the_str,"=") == 0) {
        mode = EQU;
        is_set = true;
      }

      if (strcmp(the_str,"+=") == 0) {
        mode = SUM;
        is_set = true;
      }

      if (strcmp(the_str,"-=") == 0) {
        mode = NEG;
        is_set = true;
      }

      if (is_set == false) {
        error("Non-valid string in arg 4!");
        return oct_retval;
      }
    }
    free(the_str);
    break;

  default:
    error("conv_p requires 2 to 4 input arguments!");
    return oct_retval;
    break;
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
    error("Argument 2 must be a vector or a matrix with the same number of rows as arg 1!");
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

    in_place = false;

    //
    // Create/get output matrix.
    //

    Matrix Ymat(A_M+B_M-1, A_N);
    Y = Ymat.fortran_vec();

    //
    // Call the CONV subroutine.
    //

    running = true;

    if (nthreads>1) { // Use threads

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
        D[thread_n].Y = Y;

        // Start the threads.
        threads[thread_n] = std::thread(smp_dream_conv_p, &D[thread_n]);

      } // for (thread_n = 0; thread_n < nthreads; thread_n++)

      // Wait for all threads to finish.
      for (thread_n = 0; thread_n < nthreads; thread_n++)
        threads[thread_n].join();

      // Free memory.
      if (D)
        free((void*) D);

    } else{			// Do not use threads

      if (B_N > 1) {		// B is a matrix.
        for (n=0; n<A_N; n++) {

          conv( &A[0+n*A_M], A_M, B, B_M, &Y[0+n*(A_M+B_M-1)],in_place,mode);

          if (running==false) {
            printf("conv_p: bailing out!\n");
            break;
          }
        }
      } else {			// B is a vector.
        for (n=0; n<A_N; n++) {

          conv( &A[0+n*A_M], A_M, B, B_M, &Y[0+n*(A_M+B_M-1)],in_place,mode);

          if (running==false) {
            printf("conv_p: bailing out!\n");
            break;
          }
        }
      }
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

    oct_retval.append(Ymat);
    return oct_retval;

  } else {

    //
    // In-place mode.
    //

    in_place = true;

    if (args(2).matrix_value().rows() != A_M+B_M-1) {
      error("Wrong number of rows in argument 3!");
      return oct_retval;
    }

    if (args(2).matrix_value().cols() != A_N) {
      error("Wrong number of columns in argument 3!");
      return oct_retval;
    }

    const Matrix Ytmp = args(2).matrix_value();

    Y = (double*) Ytmp.fortran_vec();

    //
    // Call the CONV subroutine.
    //

    running = true;

    if (nthreads>1) { // Use threads

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
        D[thread_n].Y = Y;

        // Start the threads.
        threads[thread_n] = std::thread(smp_dream_conv_p, &D[thread_n]);
      }

      // Wait for all threads to finish.
      for (thread_n = 0; thread_n < nthreads; thread_n++)
        threads[thread_n].join();

      // Free memory.
      if (D)
        free((void*) D);

    } else { // Do not use threads

      if (B_N > 1) { // B is a matrix.
        for (n=0; n<A_N; n++) {

          conv( &A[0+n*A_M], A_M, &B[0+n*B_M], B_M, &Y[0+n*(A_M+B_M-1)],in_place,mode);

          if (running==false) {
            printf("conv_p: bailing out!\n");
            break;
          }
        }
      } else {// B is a vector.

        for (n=0; n<A_N; n++) {

          conv( &A[0+n*A_M], A_M, B, B_M, &Y[0+n*(A_M+B_M-1)],in_place,mode);

          if (running==false) {
            printf("conv_p: bailing out!\n");
            break;
          }
        }
      }
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

    return oct_retval;
  }

}
