/***
*
* Copyright (C) 2006,2007,2008,2009,2014,2015,2016,2021 Fredrik Lingvall
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

// FIXME: Move to the FFT class?

#include <iostream>
#include <csignal>
#include <thread>
#include <complex> // C++

#include <fftw3.h>

#include "dream.h"
#include "affinity.h"

//
// Octave headers.
//

#include <octave/oct.h>

/***
 *
 *  Parallel (threaded) FFTW based convolution and summation.
 *
 ***/

//
// Globals
//
volatile int running;
volatile int in_place;
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
  octave_idx_type line_start;
  octave_idx_type line_stop;
  int L;
  double **A;
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

void* smp_dream_sum_fftconv(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

void add_fftconv(double **xr, octave_idx_type nidx, octave_idx_type nx, double *yr, octave_idx_type ny, double *zr, int L,
                 double      *a,  double *b, double *c,
                 std::complex<double> *af, std::complex<double> *bf, std::complex<double> *cf);


/***
 *
 * Thread function.
 *
 ***/

void* smp_dream_sum_fftconv(void *arg)
{
  DATA D = *(DATA *)arg;
  octave_idx_type    line_start=D.line_start, line_stop=D.line_stop, n;
  double **A = D.A, *B = D.B, *Y = D.Y;
  octave_idx_type A_M = D.A_M, B_M = D.B_M; //, B_N = D.B_N;
  int L = D.L;

  // Input vectors.
  double *a  = NULL, *b  = NULL;
  double *c  = NULL;

  // Fourier Coefficients.
  //fftw_complex *af  = NULL, *bf  = NULL, *cf  = NULL;
  std::complex<double> *af  = NULL, *bf  = NULL, *cf  = NULL;
  octave_idx_type fft_len;

  // Allocate space for input vectors.

  fft_len = A_M+B_M-1;

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

  for (n=line_start; n<line_stop; n++) {
    add_fftconv( A, n, A_M, B, B_M, &Y[0+n*(A_M+B_M-1)], L,
                 a,b,c,af,bf,cf);

    if (running==false) {
      octave_stdout << "sum_fftconv: thread for column " << line_start+1 << " -> " << line_stop << " bailing out!\n";
      break;
    }
  }

  /***
      if (B_N > 1) {// B is a matrix.

      for (n=line_start; n<line_stop; n++) {

      //add_fftconv( &A[0+n*A_M]), A_M, &B[0+n*B_M], B_M, &Y[0+n*(A_M+B_M-1)],
      //	   a,b,c,af,bf,cf);

      add_fftconv( A, n, A_M, &B[0+n*B_M], B_M, &Y[0+n*(A_M+B_M-1)], L,
      a,b,c,af,bf,cf);

      if (running==false) {
      printf("sum_fftconv: thread for column %d -> %d bailing out!\n",line_start+1,line_stop);
      break;
      }

      } // end-for
      } else { // B is a vector.

      for (n=line_start; n<line_stop; n++) {

      add_fftconv( A, n, A_M, B, B_M, &Y[0+n*(A_M+B_M-1)], L,
      a,b,c,af,bf,cf);

      if (running==false) {
      printf("sum_fftconv: thread for column %d -> %d bailing out!\n",line_start+1,line_stop);
      break;
      }

      } // end-for

      } // end.if
  ***/

  //
  //  Cleanup
  //

  // Free buffer memory.
  if (a)
    fftw_free(a);
  else {
    error("Error in freeing memory in sum_fftconv thread!!");
  }

  if(af)
    fftw_free(af);
  else {
    error("Error in freeing memory in  sum_fftconv thread!!");
  }

  if(b)
    fftw_free(b);
  else {
    error("Error in freeing memory in  sum_fftconv thread!!");
  }

  if (bf)
    fftw_free(bf);
  else {
    error("Error in freeing memory in  sum_fftconv thread!!");
  }

  if (c)
    fftw_free(c);
  else {
    error("Error in freeing memory in  sum_fftconv thread!!");
  }

  if (cf)
    fftw_free(cf);
  else {
    error("Error in freeing memory in  sum_fftconv thread!!");
  }

  return(NULL);
}

/***
 *
 * Convolution of two vectors and summation of the result.
 *
 ***/

void add_fftconv(double **xr, octave_idx_type nidx, octave_idx_type nx, double *yr, octave_idx_type ny, double *zr, int L,
             double      *a,  double *b, double *c,
             std::complex<double> *af, std::complex<double> *bf, std::complex<double> *cf)
{
  octave_idx_type n, k, fft_len;

  fft_len = nx+ny-1;

  // Clear output Fourier coefficeints.
  for (n=0; n < fft_len; n++)
      cf[n] = 0.0;

  // Loop aver all L inputs.
  for (k=0; k<L; k++) {

    // Copy and zero-pad arg 1.
    for (n=0; n < nx; n++) {
      a[n] = (xr[k])[n+nidx*nx];
    }
    for (n=nx; n < fft_len; n++)
      a[n] = 0.0; // Zero-pad.

    // Copy k:th column in arg 2 and zero-pad.
    for (n=0; n < ny; n++)
      b[n] = yr[n+k*ny];
    for (n=ny; n < fft_len; n++)
      b[n] = 0.0; // Zero-pad.

    // Fourier transform xr.
    fftw_execute_dft_r2c(p_forward,a,reinterpret_cast<fftw_complex*>(af));

    // Fourier transform yr.
    fftw_execute_dft_r2c(p_forward,b,reinterpret_cast<fftw_complex*>(bf));

    // Do the filtering and add the results
    for (n = 0; n < fft_len; n++) {
      cf[n] += (af[n] * bf[n])  / ((double) (fft_len));
    }
  } // for - K

  //
  // Compute the inverse DFT of the summed and filtered data.
  //

  fftw_execute_dft_c2r(p_backward,reinterpret_cast<fftw_complex*>(cf),c);

  // Copy data to output matrix.
  if (in_place == false) {

    //for (n = 0; n < fft_len; n++)
    //  zr[n] = c[n];
    memcpy(zr,c,fft_len*sizeof(double));

  } else { // in-place add operation.

    for (n = 0; n < fft_len; n++)
      zr[n] += c[n];

    // in-place '=' operation.
    //for (n = 0; n < fft_len; n++)
    //  zr[n] = c[n];
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

// void fftw_forget_wisdom(void);

/***
 *
 * Octave (oct) gateway function for SUM_FFTCONV.
 *
 ***/

DEFUN_DLD (sum_fftconv, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {}  Y = sum_fftconv(A,B,wisdom_str);\n\
or sum_fftconv(A,B,Y,wisdom_str);\n\
\n\
SUM_FFTCONV - Computes (using parallel threaded processing) the sum of one dimensional\n\
convolutions of the columns in each matrix in the 3D matrix A\n\
with the corresponding columns in the matrix B.\n\
\n\
@unnumberedsec Normal mode\n\
\n\
In normal mode sum_fftconv performs an operation similar to:\n\
@verbatim\n\
\n\
YF = zeros(M+K-1,N);\n\
for l=1:L\n\
  for n=1:N\n\
    YF(:,n) = YF(:,n) + fft(A(:,n,l),M+K-1).* fft(B(:,l),M+K-1);\n\
  end\n\
end\n\
Y = real(ifft(Y));\n\
\n\
@end verbatim\n\
using threaded processing. The computations are performed using FFT:s from the FFTW library.\n\
\n\
Input parameters:\n\
\n\
@table @code\n\
@item A\n\
  An MxNxL 3D matrix.\n\
@item B\n\
  A KxL matrix.\n\
@item wisdom_str\n\
Optional parameter. If the wisdom_str parameter is not supplied then fftconv calls fftw wisdom plan\n\
functions before performing any frequency domain operations. This overhead can be avoided by supplying\n \
a pre-computed fftw wisdom string. For more information see the fftw user manunal\n\
available at @url{http://www.fftw.org}.\n\
@end table\n\
\n\
The wisdom_str can be obtained using the fftconv_p function. A typical example is,\n\
@verbatim\n\
\n\
 % Compute a new fftw wisdom string.\n\
[tmp,wisdom_str]  = fftconv_p(A(:,1,1),B(:,1));\n\
\n\
for i=1:N\n\
\n\
  % Do some stuff here.\n\
\n\
  Y = sum_fftconv(A,B,wisdom_str);\n\
end\n\
\n\
@end verbatim\n\
where the overhead of calling fftw plan functions is now avoided inside the for loop.\n\
\n\
Output parameter:\n\
\n\
@table @code\n\
@item Y\n\
  The (M+K-1)xN output matrix.\n\
@end table\n\
\n\
@unnumberedsec In-place mode\n\
\n\
In in-place mode sum_fftconv performs the operations in-place on a pre-allocated\n\
matrix. Here sum_fftconv do not have any output arguments and the\n\
results are instead stored directly in the pre-allocated (M+K-1)xN input matrix Y. A typical usage is:\n\
@verbatim\n\
\n\
Y = zeros(M+K-1,N);% Allocate space for Y.\n\
\n\
for i=1:N\n\
\n\
  % Do some stuff here.\n\
\n\
   sum_fftconv(A,B,Y,wisdom_str);\n\
end\n\
\n\
@end verbatim\n\
where memory allocation of the (possible large) matrix Y now is avoided inside the for loop.\n\
\n\
NOTE: A side-effect of using in-place mode is that if a copy Y2 of Y is made\n\
then both Y2 and Y will be altered by sum_fftconv. That is, by performing,\n\
@verbatim\n\
\n\
Y  = zeros(M+K-1,N);\n\
Y2 = Y; % No actual copy of data here.\n\
sum_fftconv(A,B,Y,wisdom_str);\n\
\n\
@end verbatim\n\
@noindent then both Y and Y2 will be changed (since Octave do not make a new copy of the data\n\
in Y2 unless Y2 is changed before the sum_fftconv call).\n\
\n\
sum_fftconv is a part of the DREAM Toolbox available at\n\
@url{http://www.signal.uu.se/Toolbox/dream/}.\n\
\n\
Copyright @copyright{} 2006-2021 Fredrik Lingvall.\n\
@seealso {fftconv, conv, fftconv, conv, fftw_wisdom}\n\
@end deftypefn")
{
  double **A,*B, *Y;
  sighandler_t   old_handler, old_handler_abrt, old_handler_keyint;
  octave_idx_type line_start, line_stop, A_M, A_N, B_M, B_N, n ,k;
  int             A_L;
  std::thread    *threads;
  octave_idx_type thread_n, nthreads;
  DATA   *D;
  // Input vectors (only used for creating fftw plans).
  double *a  = NULL, *b  = NULL;
  double *c  = NULL;
  // Fourier Coefficients (only used for creating fftw plans).
  //fftw_complex *af  = NULL, *bf  = NULL, *cf  = NULL;
  std::complex<double> *af  = NULL, *bf  = NULL, *cf  = NULL;
  octave_idx_type fft_len, return_wisdom = false, load_wisdom = false;
  char *wisdom_str = NULL;
  int buflen;
  octave_value_list oct_retval;
  dim_vector dv;
  int dims;

  in_place = false;

  int nrhs = args.length ();

  //
  // Set the method which fftw computes plans
  //

  // If we want to save a plan (in the second output arg)
  // then use the more time-consuming MEAUSURE method.
  if (nlhs == 2) {
    plan_method = 3; // 3 = MEASURE.
  }

  // Check for proper inputs arguments.

  switch (nrhs) {

  case 0:
  case 1:
    error("sum_fftconv requires 2 to 4 input arguments!");
    return oct_retval;
    break;

  case 2:
    if (nlhs > 1) {
      error("Too many output arguments for sum_fftconv!");
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
      wisdom_str = (char*) fftw_malloc(buflen * sizeof(char));
      //mxGetString(prhs[3], wisdom_str, buflen); // Obsolete in Matlab 7.x
      for ( n=0; n<buflen; n++ ) {
        wisdom_str[n] = strin[n];
      }

      if (!strcmp("wisdom",wisdom_str)) {
        error("The string in arg 4 do not seem to be in fftw wisdom format!");
        return oct_retval;
      }
      else {
        load_wisdom = true;
      }
    } else { // 4th arg not a string
      fftw_forget_wisdom(); // Clear wisdom history (a new wisdom will be created below).
      if (nlhs > 0) {
        error("No output arguments required for sum_fftconv in in-place operating mode!");
        return oct_retval;
      }
    }
    break;

  case 4:
    if ( args(3).is_string() ) { // 4th arg is a fftw wisdom string.

      std::string strin = args(3).string_value();
      buflen = strin.length();
      wisdom_str = (char*) fftw_malloc(buflen * sizeof(char));
      for ( n=0; n<buflen; n++ ) {
        wisdom_str[n] = strin[n];
      }

      if (!strcmp("wisdom",wisdom_str)) {
        error("The string in arg 4 do not seem to be in fftw wisdom format!");
        return oct_retval;
      }
      else {
        load_wisdom = true;
      }
    }
    else { // 4th arg not a string
      error("Argument 4 do not seem to be in fftw wisdom string format!");
      return oct_retval;
    }
    break;

  default:
    error("sum_fftconv requires 2 to 4 input arguments!");
    return oct_retval;
    break;
  }

  const NDArray tmp0 = args(0).array_value();
  dims = args(0).ndims();
  if (dims != 3) {
    error("Argument 1 should be a 3D Matrix\n");
    return oct_retval;
  }

  dv  = args(0).dims();
  A_M = dv(0);
  A_N = dv(1);
  A_L = dv(2);

  // Store pointers to the L A-matrices in a vector.
  A = (double**) malloc(A_L*sizeof(double*));
  for (k=0; k<A_L; k++)
    A[k] = (double*) &(tmp0.fortran_vec()[A_M*A_N*k]);

  const Matrix tmp1 = args(1).matrix_value();
  B_M = tmp1.rows();
  B_N = tmp1.cols();

  B = (double*) tmp1.fortran_vec();

  if (A_L != B_N) {
    error("3rd dimension of arg 1 must match the number of columns in arg 2\n");
    return oct_retval;
  }

  // Check that arg 2.
  //if ( B_M != 1 && B_N !=1 && B_N != A_N) {
  //  error("Argument 2 must be a vector or a matrix with the same number of rows as arg 1!");
  //  return oct_retval;
  //}

  if (  B_M == 1 || B_N == 1 ) { // B is a vector.
    B_M = B_M*B_N;
    B_N = 1;
  }

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11)
  nthreads = std::thread::hardware_concurrency();

  // We cannot have more threads than the number of observation points.
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
    if (!fftw_import_wisdom_from_string(wisdom_str)) {
      error("Failed to load fftw wisdom!");
      return oct_retval;
    } else
      fftw_free(wisdom_str); // Clean up.
  }

  // 1)
  if (plan_method == 1) {
    p_forward = fftw_plan_dft_r2c_1d(fft_len, a, reinterpret_cast<fftw_complex*>(af), FFTW_EXHAUSTIVE);
    p_backward = fftw_plan_dft_c2r_1d(fft_len, reinterpret_cast<fftw_complex*>(cf), c,FFTW_EXHAUSTIVE);
  }

  // 2)
  if (plan_method == 2) {
    p_forward = fftw_plan_dft_r2c_1d(fft_len, a, reinterpret_cast<fftw_complex*>(af), FFTW_PATIENT);
    p_backward = fftw_plan_dft_c2r_1d(fft_len, reinterpret_cast<fftw_complex*>(cf), c, FFTW_PATIENT);
  }

  // 3)
  if (plan_method == 3) {
    p_forward = fftw_plan_dft_r2c_1d(fft_len, a, reinterpret_cast<fftw_complex*>(af), FFTW_MEASURE);
    p_backward = fftw_plan_dft_c2r_1d(fft_len, reinterpret_cast<fftw_complex*>(cf), c, FFTW_MEASURE);
  }

  // 4)
  if (plan_method == 4) {
    p_forward = fftw_plan_dft_r2c_1d(fft_len, a, reinterpret_cast<fftw_complex*>(af),  FFTW_ESTIMATE);
    p_backward = fftw_plan_dft_c2r_1d(fft_len, reinterpret_cast<fftw_complex*>(cf), c,  FFTW_ESTIMATE);
  }

  //
  // Normal (non in-place) mode.
  //

  if (nrhs == 2 ||  (nrhs == 3 && load_wisdom)) {

    in_place = false;

    //
    //
    // Normal (non in-place) mode.
    //
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
        D[thread_n].L =  A_L;

        // Start the threads.
        threads[thread_n] = std::thread(smp_dream_sum_fftconv, &D[thread_n]);

      } // for (thread_n = 0; thread_n < nthreads; thread_n++)

      // Wait for all threads to finish.
      for (thread_n = 0; thread_n < nthreads; thread_n++)
        threads[thread_n].join();

      // Free memory.
      if (D)
        free((void*) D);

    } else { // Do not use threads.

      for (n=0; n<A_N; n++) {

        add_fftconv( A,n, A_M, B, B_M, &Y[0+n*(A_M+B_M-1)], A_L,
                     a,b,c,af,bf,cf);

        if (running==false) {
          printf("sum_fftconv: bailing out!\n");
          break;
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

    free(A);

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

    oct_retval.append(Ymat);
    return oct_retval;
  }

  //
  // In-place mode.
  //

  if ( (nrhs == 3 && !load_wisdom) || nrhs == 4 ) {

    in_place = true;


    if (args(2).matrix_value().rows() != A_M+B_M-1) {
      error("Wrong number of rows in argument 3!");
      return oct_retval;
    }

    if (args(2).matrix_value().cols() != A_N) {
      error("Wrong number of columns in argument 3!");
      return oct_retval;
    }

    const Matrix Ytmp = args(3).matrix_value();
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
        D[thread_n].L =  A_L;

        // Start the threads.
        threads[thread_n] = std::thread(smp_dream_sum_fftconv, &D[thread_n]);
        set_dream_thread_affinity(thread_n, nthreads, threads);
      }

      // Wait for all threads to finish.
      for (thread_n = 0; thread_n < nthreads; thread_n++)
        threads[thread_n].join();

      // Free memory.
      if (D)
        free((void*) D);

    } else { // Do not use threads.

      for (n=0; n<A_N; n++) {

        add_fftconv( A,n, A_M, B, B_M, &Y[0+n*(A_M+B_M-1)], A_L,
                     a,b,c,af,bf,cf);

        if (running==false) {
          printf("sum_fftconv: bailing out!\n");
          break;
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

    free(A);

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

    // Return the FFTW Wisdom so that the plans can be re-used.
    if (return_wisdom) {
      wisdom_str = fftw_export_wisdom_to_string();
      buflen = strlen(wisdom_str);

      std::string cmout( buflen, ' ' );
      cmout.insert( buflen, (const char*) wisdom_str);

      // Add to output args.
      oct_retval.append( cmout);

      fftw_free(wisdom_str);
    }

    return oct_retval;
  }

  return oct_retval; // Just to fix compiler warnings.
}
