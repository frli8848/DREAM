/***
*
* Copyright (C) 2021,2023,2025 Fredrik Lingvall
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

//
// Octave headers.
//

#include <octave/oct.h>

#include "dream.h"
#include "affinity.h"

/***
 *
 *  Parallel (threaded) matrix and sub-matrix addition.
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
  octave_idx_type r;
  octave_idx_type k;
  octave_idx_type len;
} DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//

void* smp_dream_add_p(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

/***
 *
 * Thread function.
 *
 ***/

void* smp_dream_add_p(void *arg)
{
  DATA D = *(DATA *)arg;
  octave_idx_type    col_start=D.col_start, col_stop=D.col_stop, n;
  double *A = (double*) D.A, *B = D.B;
  octave_idx_type A_M = D.A_M, B_M = D.B_M, r = D.r, k = D.k, len = D.len;

  for (n=col_start; n<col_stop; n++) {

    for (octave_idx_type l=0; l<len; l++) {
      A[r+(k+n)*A_M + l] += B[l+n*B_M];
    }

    if (running==false) {
      octave_stdout << "add_p: thread for column " << col_start+1 << " -> " << col_stop << " bailing out!\n";
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
 * Octave (OCT) gateway function for add_p.
 *
 ***/

DEFUN_DLD (add_p, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} add_p(A,col_idx,row_idx,B);\n\
\n\
ADD_P Performs parallel (inplace) add of the matrix B and the matrix A using threads.\n\
\n\
Input parameters:\n\
\n\
@table @samp\n\
@item A\n\
An MxN matrix.\n\
@item row_idx\n\
row_idx is a two element vector where row_idx(1) is the start index \n\
@item col_idx\n\
col_idx is a two element vector where col_idx(1) is the start index \n\
@item B\n\
An (row_idx(2)-row_idx(1)) x (col_idx(2)-col_idx(1)) matrix.\n\
@end table\n\
\n\
Output parameter:\n\
\n\
@table @samp\n\
@item -\n\
There are no output arguments for @code{add_p}.\n\
@end table\n\
\n\
@code{add_p} is an oct-function that is a part of the DREAM Toolbox available at @url{https://github.com/frli8848/DREAM}.\n\
\n\
Addright @addright{} 2006-2023 Fredrik Lingvall.\n\
@end deftypefn")
{
  double         *A,*B;
  sighandler_t   old_handler, old_handler_abrt, old_handler_keyint;
  octave_idx_type col_start, col_stop, A_M, A_N, B_M, B_N;
  octave_idx_type r_M, r_N, k_M, k_N;
  double         *r, *k;
  DATA           *D=nullptr;
  std::thread     *threads;
  octave_idx_type  thread_n, nthreads;
  octave_value_list oct_retval;

  int nrhs = args.length ();

  // Check for proper number of arguments

  if (nrhs != 4) {
    error("add_p requires 4 input arguments!");
    return oct_retval;
  }

  if (nlhs > 0) {
    error("Too many output arguments for add_p, none required!");
    return oct_retval;
  }

  const Matrix tmp = args(0).matrix_value();
  A_M = tmp.rows();
  A_N = tmp.cols();
  A = (double*) tmp.data();

  const Matrix tmp2 = args(1).matrix_value();
  r_M = tmp2.rows();
  r_N = tmp2.cols();
  r = (double*) tmp2.data();

  const Matrix tmp3 = args(2).matrix_value();
  k_M = tmp3.rows();
  k_N = tmp3.cols();
  k = (double*) tmp3.data();

  const Matrix tmp4 = args(3).matrix_value();
  B_M = tmp4.rows();
  B_N = tmp4.cols();
  B = (double*) tmp4.data();

  // Check that arg 2.
  if ( r_M * r_N !=2 ) {
    error("Argument 2 must be a 2 element vector!");
    return oct_retval;
  }

  if ( (octave_idx_type) r[1] < 1) {
    error("1st element of argument 2 must be > 1!");
    return oct_retval;
  }

  if (( (octave_idx_type) r[0] > A_M) || ( (octave_idx_type) r[1] > A_M)) {
    error("One element of argument 2 exceeds row dimension of arg 1!");
    return oct_retval;
  }

  // Check that arg 3.
  if ( k_M * k_N !=2 ) {
    error("Argument 3 must be a 2 element vector!");
    return oct_retval;
  }

  if ( (octave_idx_type) k[1] < 1) {
    error("1st element of argument 3 must be > 1!");
    return oct_retval;
  }

  if (( (octave_idx_type) k[0] > A_N) || ( (octave_idx_type) k[1] > A_N)) {
    error("One element of argument 3 exceeds colomn dimension of arg 1!");
    return oct_retval;
  }

  // Check that arg 4.
  if (B_M < (octave_idx_type) (r[1]-r[0])+1 ) {
    error("Argument 4 has not the number of rows as given by arg 2!");
    return oct_retval;
  }

  if (B_N != (octave_idx_type) (k[1]-k[0])+1 ) {
    error("Argument 4 has not the number of columns as given by arg 3!");
    return oct_retval;
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
  // Call the add subroutine.
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

    col_start = ((octave_idx_type) thread_n) * A_N/nthreads;
    col_stop =  ((octave_idx_type) thread_n+1) * A_N/nthreads;

    // Init local data.
    D[thread_n].col_start = col_start; // Local start index;
    D[thread_n].col_stop = col_stop; // Local stop index;
    D[thread_n].A = A;
    D[thread_n].A_M = A_M;
    D[thread_n].A_N = A_N;
    D[thread_n].B = B;
    D[thread_n].B_M = B_M;
    D[thread_n].B_N = B_N;
    D[thread_n].r = (octave_idx_type) (r[0]-1);
    D[thread_n].k = (octave_idx_type) (k[0]-1);
    D[thread_n].len = (octave_idx_type) (r[1]-r[0]+1);

    // Start the threads.
    threads[thread_n] = std::thread(smp_dream_add_p, &D[thread_n]);
    set_dream_thread_affinity(thread_n, nthreads, threads);
  }

  // Wait for all threads to finish.
  for (thread_n = 0; thread_n < nthreads; thread_n++)
    threads[thread_n].join();

  // Free memory.
  free((void*) D);

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
