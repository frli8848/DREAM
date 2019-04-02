/***
*
* Copyright (C) 2006,2007,2008,2009,2010,2012,2014,2016,2019 Fredrik Lingvall
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


#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <thread>
#include <signal.h>

//#include <complex.h> // C (iso C99).
#include <complex> // C++
#include <stdio.h>
#include <fftw3.h>

//
// Octave headers.
//

#include <octave/oct.h>

//
// Macros
//

#define mxGetM(N)   args(N).matrix_value().rows()
#define mxGetN(N)   args(N).matrix_value().cols()
#define mxIsChar(N) args(N).is_string()

#define EQU 0
#define SUM 1
#define NEG 2

// ToDo: Add check if input pars have zero length
// ToDo: Add help text for in-place mode.
// We pobably need fftw >= 3.2.x to handle 64-bit array indexing.
// ToDO If we have fftw 3.2.x we can use fftw_plan_dft_r2c_1d_64 etc.

//#include <mutex>
//std::mutex print_mutex;         // Debug print mutex (C++11)


#include "dream.h"
#include "affinity.h"

/***
 *
 *  Parallel (threaded) FFTW based convolution.
 *
 ***/

//
// Globals
//

volatile int running;
volatile int in_place;
int mode = EQU;
int plan_method = 4; // Default to ESTIMATE method.

// FFTW plans.
//static  fftw_plan    p_forward   = NULL;
//static  fftw_plan    p_backward  = NULL;
fftw_plan    p_forward;
fftw_plan    p_backward;

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
} DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//
void* smp_fftconv_process(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

void fftconv(double *xr, octave_idx_type nx, double *yr, octave_idx_type ny, double *zr,
             double      *a,  double *b, double *c,
             std::complex<double> *af, std::complex<double> *bf, std::complex<double> *cf);

/***
 *
 * Thread function.
 *
 ***/

void* smp_fftconv_process(void *arg)
{
  DATA D = *(DATA *)arg;
  octave_idx_type    col_start=D.col_start, col_stop=D.col_stop, n;
  double *A = D.A, *B = D.B, *Y = D.Y;
  octave_idx_type A_M = D.A_M, B_M = D.B_M, B_N = D.B_N;

  // Input vectors.
  double *a  = NULL, *b  = NULL;
  double *c  = NULL;

  // Fourier Coefficients.
  //fftw_complex *af  = NULL, *bf  = NULL, *cf  = NULL;
  std::complex<double> *af  = NULL, *bf  = NULL, *cf  = NULL;
  octave_idx_type fft_len;

  // Allocate space for input vectors.

  fft_len = A_M+B_M-1;


  //
  // fftw_malloc may not be thread safe! Check this!!!
  //

  a = (double*) fftw_malloc(sizeof(double)*2*(fft_len/2+1));
  //af = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));
  af = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));

  b = (double*) fftw_malloc(sizeof(double)*2*(fft_len/2+1));
  bf = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));

  // Allocate space for output vector.

  c = (double*) fftw_malloc(sizeof(double)*2*(fft_len/2+1));
  //cf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));
  cf = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));

  //
  // Do the convolution.
  //
  size_t k = 0;

  if (B_N > 1) {// B is a matrix.

    n=col_start;
    double *Ap = &A[0+n*A_M];
    double *Bp = &B[0+n*B_M];
    double *Yp = &Y[0+n*(A_M+B_M-1)];
    for (n=col_start; n<col_stop; n++) {

      fftconv( Ap, A_M, Bp, B_M, Yp, a,b,c,af,bf,cf);

      Ap += A_M;
      Bp += B_M;
      Yp += (A_M+B_M-1);

      if (running==false) {
        octave_stdout << "fftconv_p: thread for column " << col_start+1 << " -> " << col_stop << " bailing out!\n";
        break;
      }

      k++;

    } // end-for
  } else { // B is a vector.

    n=col_start;
    double *Ap = &A[0+n*A_M];
    double *Yp = &Y[0+n*(A_M+B_M-1)];
    for (n=col_start; n<col_stop; n++) {

      fftconv( Ap, A_M, B, B_M, Yp, a,b,c,af,bf,cf);

      Ap += A_M;
      Yp += (A_M+B_M-1);

      if (running==false) {
        octave_stdout << "fftconv_p: thread for column " << col_start+1 << " -> " << col_stop << " bailing out!\n";
        break;
      }

      k++;

    } // end-for
  } // end-if

  /*
  {
    std::lock_guard<std::mutex> lk(print_mutex);
    octave_stdout << "col_start: " << col_start
                  << " col_stop: " << col_stop
                  << " k: " << k
                  << std::endl;

  }
  */

  //
  //  Cleanup
  //

  // Free buffer memory.
  if (a)
    fftw_free(a);
  else {
    error("Error in freeing memory in fftconv_p thread!!");
  }

  if(af)
    fftw_free(af);
  else {
    error("Error in freeing memory in  fftconv_p thread!!");
  }

  if(b)
    fftw_free(b);
  else {
    error("Error in freeing memory in  fftconv_p thread!!");
  }

  if (bf)
    fftw_free(bf);
  else {
    error("Error in freeing memory in  fftconv_p thread!!");
  }

  if (c)
    fftw_free(c);
  else {
    error("Error in freeing memory in  fftconv_p thread!!");
  }

  if (cf)
    fftw_free(cf);
  else {
    error("Error in freeing memory in  fftconv_p thread!!");
  }

  return(NULL);
}



/***
 *
 * Convolution of two vectors.
 *
 ***/

void fftconv(double *xr, octave_idx_type nx, double *yr, octave_idx_type ny, double *zr,
             double      *a,  double *b, double *c,
             std::complex<double> *af, std::complex<double> *bf, std::complex<double> *cf)
{
  octave_idx_type n, fft_len;

  fft_len = nx+ny-1;

  //
  // Copy and zero-pad.
  //
  for (n=0; n < nx; n++) {
    a[n] = xr[n];
  }

  for (n=nx; n < fft_len; n++)
    a[n] = 0.0; // Zero-pad.

  for (n=0; n < ny; n++)
    b[n] = yr[n];
  for (n=ny; n < fft_len; n++)
    b[n] = 0.0; // Zero-pad.

  // Fourier transform xr.
  fftw_execute_dft_r2c(p_forward,a,reinterpret_cast<fftw_complex*>(af));

  // Fourier transform yr.
  fftw_execute_dft_r2c(p_forward,b,reinterpret_cast<fftw_complex*>(bf));

  // Do the filtering.
  for (n = 0; n < fft_len; n++) {
    cf[n] = (af[n] * bf[n])  / ((double) (fft_len));
  }

  //
  // Compute the inverse DFT of the filtered data.
  //

  fftw_execute_dft_c2r(p_backward,reinterpret_cast<fftw_complex*>(cf),c);

  // Copy data to output matrix.
  if (in_place == false) {

    //for (n = 0; n < fft_len; n++)
    //  zr[n] = c[n];
    memcpy(zr,c,fft_len*sizeof(double));

  } else { // in-place

    switch (mode) {

    case EQU:
      // in-place '=' operation.
      //for (n = 0; n < fft_len; n++) {
      //zr[n] = c[n];
      //}
      memcpy(zr,c,fft_len*sizeof(double));
      break;

    case SUM:
      // in-place '+=' operation.
      for (n = 0; n < fft_len; n++) {
        zr[n] += c[n];
      }
      break;

    case NEG:
      // in-place '-=' operation.
      for (n = 0; n < fft_len; n++) {
        zr[n] -= c[n];
      }
      break;

    default:
      // in-place '=' operation.
      //for (n = 0; n < fft_len; n++) {
      //zr[n] = c[n];
      //}
      memcpy(zr,c,fft_len*sizeof(double));
      break;
    }
  }
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
 * Octave (oct) gateway function for FFTCONV.
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
Optional parameter. If the wisdom_str_in parameter is not supplied then @code{fftconv_p} calls fftw wisdom plan\n\
functions before performing any frequency domain operations. This overhead can be avoided by supplying\n\
a pre-computed fftw wisdom string wisdom_str_in. For more information see the fftw user manunal\n\
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
will call fftw wisdom plan functions and return the wisdom string which then can be used to\n\
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
NOTE: fftconv_p requires the FFTW library version 3 @url{http://www.fftw.org}.\n\
\n\
fftconv_p is a part of the DREAM Toolbox available at\n\
@url{http://www.signal.uu.se/Toolbox/dream/}.\n\
\n\
Copyright @copyright{} 2006-2016 Fredrik Lingvall.\n\
@seealso {conv_p, fftconv, conv, fftw_wisdom}\n\
@end deftypefn")
{
  double *A,*B, *Y;
  octave_idx_type col_start, col_stop, A_M, A_N, B_M, B_N, n;
  sighandler_t    old_handler, old_handler_abrt, old_handler_keyint;
  std::thread     *threads;
  octave_idx_type  thread_n, nthreads;
  DATA   *D;
  // Input vectors (only used for creating fftw plans).
  double *a  = NULL, *b  = NULL;
  double *c  = NULL;
  // Fourier Coefficients (only used for creating fftw plans).
  //fftw_complex *af  = NULL, *bf  = NULL, *cf  = NULL;
  std::complex<double> *af  = NULL, *bf  = NULL, *cf  = NULL;
  octave_idx_type fft_len, return_wisdom = false, load_wisdom = false;
  char *the_str = NULL;
  int buflen, is_set = false;
  octave_value_list oct_retval;

  in_place = false;

  int nrhs = args.length ();

  //
  // Set the method which fftw computes plans
  //

  // If we want to save a plan (in the second output arg)
  // then use the more time-consuming MEAUSURE method.
  if (nlhs == 2) {
    plan_method = 3; // 3 = MEASURE. This takes too long on long fft:s
    //plan_method = 4; // 4 = ESTIMATE.
  }

  // Check for proper inputs arguments.

  switch (nrhs) {

  case 0:
  case 1:
    error("fftconv_p requires 2 to 5 input arguments!");
    return oct_retval;
    break;

  case 2:
    if (nlhs > 2) {
      error("Too many output arguments for fftconv_p!");
      return oct_retval;
    }
    if (nlhs == 2)
      return_wisdom = true;

    break;

  case 3:
    if ( args(2).is_string() ) { // 3rd arg is a fftw wisdom string.
      std::string strin = args(2).string_value();

      buflen = strin.length();

      //buflen = mxGetM(prhs[3])*mxGetN(prhs[3]);
      the_str = (char*) fftw_malloc(buflen * sizeof(char));
      //mxGetString(prhs[3], the_str, buflen); // Obsolete in Matlab 7.x
      for ( n=0; n<buflen; n++ ) {
        the_str[n] = strin[n];
      }

      //
      // If 3rd arg is a string then only a wisdom string is valid..
      //

      if (!strcmp("fftw_wisdom",the_str)) {
        error("The string in arg 3 do not seem to be in a fftw wisdom format!");
        return oct_retval;
      }
      else {
        load_wisdom = true;
      }
    } else { // 3rd arg not a string then assume in-place mode.
      fftw_forget_wisdom(); // Clear wisdom history (a new wisdom will be created below).
      if (nlhs > 0) {
        error("3rd arg is not a fftw wisdom string and in-place mode is assumed. But then there should be no output args!");
        return oct_retval;
      }
    }
    break;

  case 4:  // In-place mode if >= 4 args.
    if ( args(3).is_string() ) { // 4th arg is a string ('=','+=','-=', or wisdom).

      std::string strin = args(3).string_value();
      buflen = strin.length();
      the_str = (char*) fftw_malloc(buflen * sizeof(char));
      for ( n=0; n<buflen; n++ ) {
        the_str[n] = strin[n];
      }

      // Valid strings are:
      //  '='  : In-place replace mode.
      //  '+=' : In-place add mode.
      //  '-=' : In-place sub mode.
      //  wisdom string.

      is_set = false;

      if (strncmp(the_str,"=",1) == 0) {
        mode = EQU;
        is_set = true;
      }

      if (strncmp(the_str,"+=",2) == 0) {
        mode = SUM;
        is_set = true;
      }

      if (strncmp(the_str,"-=",2) == 0) {
        mode = NEG;
        is_set = true;
      }

      if (is_set == false) {
        if (strcmp("fftw_wisdom",the_str) < 0 ) {
          error("Non-valid string in arg 4!");
          return oct_retval;
        } else {
          load_wisdom = true;
        }
      }
    } else { // 4th arg not a string
      error("Argument 4 is not a valid string format!");
      return oct_retval;
    }
    break;


  case 5: // In-place mode if input 5 args.
    if ( args(4).is_string() ) { // 5:th arg is a string (fftw wisdom).

      // Read the wisdom string.
      std::string strin = args(4).string_value();
      buflen = strin.length();
      the_str = (char*) fftw_malloc(buflen * sizeof(char));
      for ( n=0; n<buflen; n++ ) {
        the_str[n] = strin[n];
      }

      if (strcmp("fftw_wisdom",the_str) < 0 ) {
        error("The string in 5th arg do not seem to be in a fftw wisdom format!");
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

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11)
  nthreads = std::thread::hardware_concurrency();

  if(const char* env_p = std::getenv("OMP_NUM_THREADS")) {
    int omp_threads = std::stoul(env_p);
    if (omp_threads < nthreads) {
      nthreads = omp_threads;
    }
  }

  // nthreads can't be larger then the number of columns in the A matrix.
  if (nthreads > A_N) {
    nthreads = A_N;
  }

  /*
  {
    std::lock_guard<std::mutex> lk(print_mutex);
    octave_stdout << "A_M: " <<  A_M
                  << " A_N: " <<  A_N
                  << " B_M: " <<  B_M
                  << " A_N: " <<  A_N
                  << " nthreads: "  << nthreads <<std::endl;
  }
  */

  //
  // Register signal handlers.
  //

  if ((old_handler = signal(SIGTERM, &sighandler)) == SIG_ERR) {
    printf("Couldn't register signal handler.\n");
  }

  if ((old_handler_abrt=signal(SIGABRT, &sighandler)) == SIG_ERR) {
    printf("Couldn't register signal handler.\n");
  }

  if ((old_handler_keyint=signal(SIGINT, &sighandler)) == SIG_ERR) {
    printf("Couldn't register signal handler.\n");
  }

  // Allocate space for (temp) input/output vectors (only used for creating plans).

  fft_len = A_M+B_M-1;

  a = (double*) fftw_malloc(sizeof(double)*2*(fft_len/2+1));
  //af = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));
  af = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));

  b = (double*) fftw_malloc(sizeof(double)*2*(fft_len/2+1));
  //bf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));
  bf = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));

  c = (double*) fftw_malloc(sizeof(double)*2*(fft_len/2+1));
  //cf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));
  cf = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));

  //
  // Init the FFTW plans.
  //

  if(load_wisdom) {

    if (!fftw_import_wisdom_from_string(the_str)) {
      error("Failed to load fftw wisdom!");
      return oct_retval;
    } else
      fftw_free(the_str); // Clean up.
  }

  // 1) Very slow
  if (plan_method == 1) {
    p_forward = fftw_plan_dft_r2c_1d(fft_len, a, reinterpret_cast<fftw_complex*>(af), FFTW_EXHAUSTIVE);
    p_backward = fftw_plan_dft_c2r_1d(fft_len, reinterpret_cast<fftw_complex*>(cf), c, FFTW_EXHAUSTIVE);
  }

  // 2) Slow
  if (plan_method == 2) {
    p_forward = fftw_plan_dft_r2c_1d(fft_len, a, reinterpret_cast<fftw_complex*>(af), FFTW_PATIENT);
    p_backward = fftw_plan_dft_c2r_1d(fft_len, reinterpret_cast<fftw_complex*>(cf), c, FFTW_PATIENT);
  }

  // 3) Too slow on long FFTs.
  if (plan_method == 3) {
    p_forward = fftw_plan_dft_r2c_1d(fft_len, a, reinterpret_cast<fftw_complex*>(af), FFTW_MEASURE);
    p_backward = fftw_plan_dft_c2r_1d(fft_len, reinterpret_cast<fftw_complex*>(cf), c, FFTW_MEASURE);
  }

  // 4)
  if (plan_method == 4) {
    p_forward = fftw_plan_dft_r2c_1d(fft_len, a, reinterpret_cast<fftw_complex*>(af),  FFTW_ESTIMATE);
    p_backward = fftw_plan_dft_c2r_1d(fft_len, reinterpret_cast<fftw_complex*>(cf), c,  FFTW_ESTIMATE);
  }

  if (nrhs == 2 || (nrhs == 3 && load_wisdom)) { // Normal mode.

    in_place = false;

    //
    // Normal (non in-place) mode.
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

        /*
        {
          std::lock_guard<std::mutex> lk(print_mutex);
          octave_stdout << "col_start: " << col_start
                        << " col_stop: " << col_stop
                        << " thread_n: "  << thread_n <<std::endl;
        }
        */

        // Start the threads.
        threads[thread_n] = std::thread(smp_fftconv_process, &D[thread_n]);
        set_dream_thread_affinity(thread_n, nthreads, threads);
      }

      //
      // Wait for all threads to finish.
      //

      for (thread_n = 0; thread_n < nthreads; thread_n++)
        threads[thread_n].join();

      // Free memory.
      if (D) {
        free((void*) D);
      }

    } else { // Do not use threads.

      if (B_N > 1) {// B is a matrix.

        for (n=0; n<A_N; n++) {

          fftconv( &A[0+n*A_M], A_M, &B[0+n*B_M], B_M, &Y[0+n*(A_M+B_M-1)],
                   a,b,c,af,bf,cf);

          if (running==false) {
            printf("fftconv_p: bailing out!\n");
            break;
          }

        } // end-for
      } else { // B is a vector.

        for (n=0; n<A_N; n++) {

          fftconv( &A[0+n*A_M], A_M, B, B_M, &Y[0+n*(A_M+B_M-1)],
                   a,b,c,af,bf,cf);

          if (running==false) {
            printf("fftconv_p: bailing out!\n");
            break;
          }

        } // end-for
      } // end-if

    }

    //
    // Restore old signal handlers.
    //

    if (signal(SIGTERM, old_handler) == SIG_ERR) {
      printf("Couldn't register old signal handler.\n");
    }

    if (signal(SIGABRT,  old_handler_abrt) == SIG_ERR) {
      printf("Couldn't register signal handler.\n");
    }

    if (signal(SIGINT, old_handler_keyint) == SIG_ERR) {
      printf("Couldn't register signal handler.\n");
    }

    if (!running) {
      error("CTRL-C pressed!\n"); // Bail out.
      return oct_retval;
    }

    // 1st output arg.
    oct_retval.append(Ymat);

    // Return the FFTW Wisdom so that the plans can be re-used.
    if (return_wisdom) {
      the_str = fftw_export_wisdom_to_string();
      buflen = strlen(the_str);

      std::string cmout( buflen, ' ' );
      cmout.insert( buflen, (const char*) the_str);

      // Add to output args.
      oct_retval.append( cmout);

      fftw_free(the_str);
    }

    // Clear temp vectors used for the FFTW plans.

    fftw_free(a);
    fftw_free(af);

    fftw_free(b);
    fftw_free(bf);

    fftw_free(c);
    fftw_free(cf);

    // Cleanup the plans.
    fftw_destroy_plan(p_forward);
    fftw_destroy_plan(p_backward);

    // Clean up FFTW
    //fftw_cleanup(); This seems to put Octave in an unstable state. Calling fftconv will crash Octave.

    return oct_retval;
  }


  if ( (nrhs == 3 && !load_wisdom) || nrhs == 4 || nrhs == 5) { // In-place mode.

    in_place = true;

    //
    //
    // In-place mode.
    //
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
    Y = (double*) Ymat.fortran_vec();

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

        // Start the threads.
        threads[thread_n] = std::thread(smp_fftconv_process, &D[thread_n]);

      } // for (thread_n = 0; thread_n < nthreads; thread_n++)

      // Wait for all threads to finish.
      for (thread_n = 0; thread_n < nthreads; thread_n++)
        threads[thread_n].join();

      // Free memory.
      if (D) {
        free((void*) D);
      }

    } else { // Do not use threads.

      if (B_N > 1) {// B is a matrix.

        for (n=0; n<A_N; n++) {

          fftconv( &A[0+n*A_M], A_M, &B[0+n*B_M], B_M, &Y[0+n*(A_M+B_M-1)],
                   a,b,c,af,bf,cf);

          if (running==false) {
            printf("fftconv_p: bailing out!\n");
            break;
          }

        } // end-for
      } else { // B is a vector.

        for (n=0; n<A_N; n++) {

          fftconv( &A[0+n*A_M], A_M, B, B_M, &Y[0+n*(A_M+B_M-1)],
                   a,b,c,af,bf,cf);

          if (running==false) {
            printf("fftconv_p: bailing out!\n");
            break;
          }

        } // end-for
      } // end-if

    }

    //
    // Restore old signal handlers.
    //

    if (signal(SIGTERM, old_handler) == SIG_ERR) {
      printf("Couldn't register old signal handler.\n");
    }

    if (signal(SIGABRT,  old_handler_abrt) == SIG_ERR) {
      printf("Couldn't register signal handler.\n");
    }

    if (signal(SIGINT, old_handler_keyint) == SIG_ERR) {
      printf("Couldn't register signal handler.\n");
    }

    if (!running) {
      error("CTRL-C pressed!\n"); // Bail out.
      return oct_retval;
    }

    // Return the FFTW Wisdom so that the plans can be re-used.
    if (return_wisdom) {
      the_str = fftw_export_wisdom_to_string();
      buflen = strlen(the_str);

      std::string cmout( buflen, ' ' );
      cmout.insert( buflen, (const char*) the_str);

      // Add to output args.
      oct_retval.append( cmout);

      fftw_free(the_str);
    }

    // Clear temp vectors used for the FFTW plans.

    fftw_free(a);
    fftw_free(af);

    fftw_free(b);
    fftw_free(bf);

    fftw_free(c);
    fftw_free(cf);

    // Cleanup the plans.
    fftw_destroy_plan(p_forward);
    fftw_destroy_plan(p_backward);

    // Clean up FFTW
    //fftw_cleanup(); This seems to put Octave in an unstable state. Calling fftconv will crash Octave.


    return oct_retval;
  }

  return oct_retval; // Just to fix compiler warnings.
}
