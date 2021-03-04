/***
*
* Copyright (C) 2006,2007,2008,2009,2016 Fredrik Lingvall
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

#include <signal.h>

#include <string.h>
#include <stdlib.h>
#include <pthread.h>

#include <uchar.h>

#include <complex> // C++
#include <stdio.h>
#include <fftw3.h>
#include <math.h>
//#include <unistd.h>

#include "dream.h"

#include "mex.h"

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
  dream_idx_type line_start;
  dream_idx_type line_stop;
  int L;
  double **A;
  dream_idx_type A_M;
  dream_idx_type A_N;
  double *B;
  dream_idx_type B_M;
  dream_idx_type B_N;
  double *Y;
} DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//

void* smp_dream_sum_fftconv_p_msvc(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

void add_fftconv(double **xr, dream_idx_type nidx, dream_idx_type nx, double *yr, dream_idx_type ny, double *zr, int L,
                 double      *a,  double *b, double *c,
                 std::complex<double> *af, std::complex<double> *bf, std::complex<double> *cf);

/***
 *
 * Thread function.
 *
 ***/

void* smp_dream_sum_fftconv_p_msvc(void *arg)
{
  DATA D = *(DATA *)arg;
  dream_idx_type    line_start=D.line_start, line_stop=D.line_stop, n;
  double **A = D.A, *B = D.B, *Y = D.Y;
  dream_idx_type A_M = D.A_M, B_M = D.B_M, B_N = D.B_N;
  int L = D.L;

  // Input vectors.
  double *a  = NULL, *b  = NULL;
  double *c  = NULL;

  // Fourier Coefficients.
  fftw_complex *af  = NULL, *bf  = NULL, *cf  = NULL;
  //std::complex<double> *af  = NULL, *bf  = NULL, *cf  = NULL;
  int fft_len;

  // Allocate space for input vectors.

  fft_len = A_M+B_M-1;

  a = (double*) fftw_malloc(sizeof(double)*2*(fft_len/2+1));
  af = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));
  //af = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));

  b = (double*) fftw_malloc(sizeof(double)*2*(fft_len/2+1));
  bf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));
  //bf = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));

  // Allocate space for output vector.
  c = (double*) fftw_malloc(sizeof(double)*2*(fft_len/2+1));
  cf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));
  //cf = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));

  //
  // Do the convolution.
  //

  for (n=line_start; n<line_stop; n++) {
    add_fftconv( A, n, A_M, B, B_M, &Y[0+n*(A_M+B_M-1)], L,
                 a,b,c,
                 reinterpret_cast<std::complex<double>*>(af),
                 reinterpret_cast<std::complex<double>*>(bf),
                 reinterpret_cast<std::complex<double>*>(cf));

    if (running==false) {
      printf("sum_fftconv_p: thread for column %d -> %d bailing out!\n",line_start+1,line_stop);
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
      printf("sum_fftconv_p: thread for column %d -> %d bailing out!\n",line_start+1,line_stop);
      break;
      }

      } // end-for
      } else { // B is a vector.

      for (n=line_start; n<line_stop; n++) {

      add_fftconv( A, n, A_M, B, B_M, &Y[0+n*(A_M+B_M-1)], L,
      a,b,c,af,bf,cf);

      if (running==false) {
      printf("sum_fftconv_p: thread for column %d -> %d bailing out!\n",line_start+1,line_stop);
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
    mexErrMsgTxt("Error in freeing memory in sum_fftconv_p thread!!");
  }

  if(af)
    fftw_free(af);
  else {
    mexErrMsgTxt("Error in freeing memory in  sum_fftconv_p thread!!");
  }

  if(b)
    fftw_free(b);
  else {
    mexErrMsgTxt("Error in freeing memory in  sum_fftconv_p thread!!");
  }

  if (bf)
    fftw_free(bf);
  else {
    mexErrMsgTxt("Error in freeing memory in  sum_fftconv_p thread!!");
  }

  if (c)
    fftw_free(c);
  else {
    mexErrMsgTxt("Error in freeing memory in  sum_fftconv_p thread!!");
  }

  if (cf)
    fftw_free(cf);
  else {
    mexErrMsgTxt("Error in freeing memory in  sum_fftconv_p thread!!");
  }

  return(NULL);
}



/***
 *
 * Convolution of two vectors and summation of the result.
 *
 ***/

void add_fftconv(double **xr, dream_idx_type nidx, dream_idx_type nx, double *yr, dream_idx_type ny, double *zr, int L,
             double      *a,  double *b, double *c,
             std::complex<double> *af, std::complex<double> *bf, std::complex<double> *cf)
{
  dream_idx_type n, k, fft_len;

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
 * Matlab (mex) gateway function for SUM_FFTCONV_P.
 *
 ***/
extern "C"
void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double **A,*B, *Y;
  sighandler_t   old_handler, old_handler_abrt, old_handler_keyint;
  pthread_t *threads;
  dream_idx_type line_start, line_stop, A_M, A_N, B_M, B_N, N, err, n, k;
  int    thread_n,  A_L;
  void   *retval;
  DATA   *D;
  // Input vectors (only used for creating fftw plans).
  double *a  = NULL, *b  = NULL;
  double *c  = NULL;
  // Fourier Coefficients (only used for creating fftw plans).
  fftw_complex *af  = NULL, *bf  = NULL, *cf  = NULL;
  dream_idx_type fft_len, return_wisdom = false, load_wisdom = false;
  char *wisdom_str;
  int buflen;
  const mwSize *dv; // This don't compile in XP with MinGW (not defined in matrix.h?!).
                    // Needs Matlab R2006b (eg. R2006b's matlab/extern/include/matrix.h)
  //const long *dv; // This gives a runtime error.
  int dims;


  in_place = false;

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
  case 2:
    mexErrMsgTxt("sum_fftconv_p requires 3 to 5 input arguments!");
    break;

  case 3:
    if (nlhs > 1) {
      mexErrMsgTxt("Too many output arguments for sum_fftconv_p!");
    }
    if (nlhs == 2)
      return_wisdom = true;

    break;

  case 4:
    if (mxIsChar(prhs[3])) { // 4th arg is a string.
      //wisdom_str = (char*) mxGetChars(prhs[3]);
      buflen = mxGetM(prhs[3])*mxGetN(prhs[3])+1;
      wisdom_str = (char*) fftw_malloc(buflen * sizeof(char));
      mxGetString(prhs[3], wisdom_str, buflen); // Obsolete in Matlab 7.x

      //
      // If 4:th arg is a string then only a wisdom string is valid.
      //

      if (strcmp("wisdom",wisdom_str) < 0) {
        mexErrMsgTxt("The string in arg 4 do not seem to be in fftw wisdom format!");
      }
      else {
        load_wisdom = true;
      }
    } else { // 4th arg not a string then assume in-place mode.
      fftw_forget_wisdom(); // Clear wisdom history (a new wisdom will be created below).
      if (nlhs > 0) {
        mexErrMsgTxt("No output arguments required for fftconv_p in in-place operating mode!");
      }
    }
    break;

  case 5:
    if (mxIsChar(prhs[4])) { // 5th arg is a string.
      //wisdom_str = (char*) mxGetChars(prhs[4]);
      buflen = mxGetM(prhs[4])*mxGetN(prhs[4])+1;
      wisdom_str = (char*) fftw_malloc(buflen * sizeof(char));
      mxGetString(prhs[4], wisdom_str, buflen); // Obsolete in Matlab 7.x ?

      if (!strcmp("wisdom",wisdom_str)) {
        mexErrMsgTxt("The string in arg 5 do not seem to be in fftw wisdom format!");

      }
      else {
        load_wisdom = true;
      }
    }
    else { // 5th arg not a string
      mexErrMsgTxt("Argument 5 do not seem to be in fftw wisdom string format!");

    }
    break;

  default:
    mexErrMsgTxt("sum_fftconv_p requires 3 to 5 input arguments!");

    break;
  }

  dims = mxGetNumberOfDimensions(prhs[0]);

  if (dims != 3) {
    mexErrMsgTxt("Argument 1 should be a 3D Matrix\n");
  }

  dv = mxGetDimensions(prhs[0]);
  A_M = dv[0];
  A_N = dv[1];
  A_L = dv[2];

  //A_M = mxGetM(prhs[0]);
  //A_N = mxGetN(prhs[0]);

  // Store pointers to the L A-matrices in a vector.
  A = (double**) malloc(A_L*sizeof(double*));
  double* tmp_p = mxGetPr(prhs[0]);
  for (k=0; k<A_L; k++)
    A[k] = &(tmp_p[A_M*A_N*k]);

  B_M = mxGetM(prhs[1]);
  B_N = mxGetN(prhs[1]);
  B = mxGetPr(prhs[1]);

  //printf(" A_M = %d, A_N = %d, A_L = %d, B_N = %d\n",  A_M,A_N,A_L,B_M);
  //if (nrhs >= 4)
  //  printf(" P_M = %d, P_N = %d\n", mxGetM(prhs[3]),mxGetN(prhs[3]));

  //  printf(" A_M = %d, A_N = %d, A_L = %d\n", (mxGetDimensions(prhs[0]))[0],(mxGetDimensions(prhs[0]))[1],(mxGetDimensions(prhs[0]))[2]);

  if (A_L != B_N) {
    mexErrMsgTxt("3rd dimension of arg 1 must match the number of columns in arg 2\n");

  }

  // Check that arg 2.
  //if ( B_M != 1 && B_N !=1 && B_N != A_N) {
  //  mexErrMsgTxt("Argument 2 must be a vector or a matrix with the same number of rows as arg 1!");
  //
  //}

  if (  B_M == 1 || B_N == 1 ) { // B is a vector.
    B_M = B_M*B_N;
    B_N = 1;
  }

  //
  // Number of threads.
  //

  // Check that arg 3 is a scalar.
  if ( mxGetM(prhs[2]) * mxGetN(prhs[2]) !=1)
    mexErrMsgTxt("Argument 3 must be a scalar!");

  //N = (int) args(2).matrix_value().fortran_vec()[0];
  N = (int) mxGetPr(prhs[2])[0];
  if (N < 1) {
    mexErrMsgTxt("Argument 3 must be larger or equal to 1 (min number of CPUs)!");

  }

  // Check n_cpus argument.
  if (N > A_N) {
    // Add a -v verbose flag to display stuff like this!
    //mexPrintf("Warning: n_cpus is larger then number of columns in first arg.\n");
    //mexPrintf("         Setting n_cpus = # cols in 1st arg!\n");
    N = A_N; //
  }

  //
  // Create/get output matrix.
  //

  if (nrhs == 3 ||  (nrhs == 4 && load_wisdom)) { // Normal mode.

    in_place = false;
    plhs[0] = mxCreateDoubleMatrix(A_M+B_M-1, A_N, mxREAL);
    Y = mxGetPr(plhs[0]);
  }
  //else { // in-place mode.

  if ( (nrhs == 4 && !load_wisdom) || nrhs == 5 ) { // In-place mode.

    in_place = true;

    if ( mxGetM(prhs[3]) != A_M+B_M-1)
      mexErrMsgTxt("Wrong number of rows in argument 4!");

    if ( mxGetN(prhs[3]) != A_N)
      mexErrMsgTxt("Wrong number of columns in argument 4!");

    Y = mxGetPr(prhs[3]);
  }

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
  af = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));
  //af = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));

  b = (double*) fftw_malloc(sizeof(double)*2*(fft_len/2+1));
  bf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));
  //bf = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));

  c = (double*) fftw_malloc(sizeof(double)*2*(fft_len/2+1));
  cf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));
  //cf = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));

  //
  // Init the FFTW plans.
  //
  if(load_wisdom) {
    if (!fftw_import_wisdom_from_string(wisdom_str)) {
      mexErrMsgTxt("Failed to load fftw wisdom!");

    } else
      fftw_free(wisdom_str); // Clean up.
  }

  // 1)
  //p_forward = fftw_plan_dft_r2c_1d(fft_len, a, af, FFTW_EXHAUSTIVE);
  //p_backward = fftw_plan_dft_c2r_1d(fft_len, cf, c,FFTW_EXHAUSTIVE);

  // 2)
  //p_forward = fftw_plan_dft_r2c_1d(fft_len, a, af, FFTW_PATIENT);
  //p_backward = fftw_plan_dft_c2r_1d(fft_len, cf, c, FFTW_PATIENT);
  //p_forward = fftw_plan_dft_r2c_1d(fft_len, a, reinterpret_cast<fftw_complex*>(af), FFTW_PATIENT);
  //p_backward = fftw_plan_dft_c2r_1d(fft_len, reinterpret_cast<fftw_complex*>(cf), c, FFTW_PATIENT);

  // 3)
  if (plan_method == 3) {
    p_forward = fftw_plan_dft_r2c_1d(fft_len, a, af, FFTW_MEASURE);
    p_backward = fftw_plan_dft_c2r_1d(fft_len, cf, c, FFTW_MEASURE);
  }

  // 4)
  if (plan_method == 4) {
    p_forward = fftw_plan_dft_r2c_1d(fft_len, a, af, FFTW_ESTIMATE);
    p_backward = fftw_plan_dft_c2r_1d(fft_len, cf, c, FFTW_ESTIMATE);
  }

  //
  // Return the FFTW Wisdom so that the plans can be re-used.
  //
  if (return_wisdom) {
    wisdom_str = fftw_export_wisdom_to_string();
    plhs[1] = mxCreateString(wisdom_str);
    fftw_free(wisdom_str);
  }

   //
  // Call the CONV subroutine.
  //

  running = true;

  if (N>1) { // Use threads

    // Allocate local data.
    D = (DATA*) malloc(N*sizeof(DATA));
    if (!D) {
      mexErrMsgTxt("Failed to allocate memory for thread data!");

    }

    // Allocate mem for the threads.
    threads = (pthread_t*) malloc(N*sizeof(pthread_t));

    if (!threads) {
      mexErrMsgTxt("Failed to allocate memory for threads!");

    }

    for (thread_n = 0; thread_n < N; thread_n++) {

      line_start = thread_n * A_N/N;
      line_stop =  (thread_n+1) * A_N/N;

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
      // Starts the threads.
      err = pthread_create(&threads[thread_n], NULL, smp_dream_sum_fftconv_p_msvc, &D[thread_n]);
      if (err != 0) {
        mexErrMsgTxt("Error when creating a new thread!\n");

      }
    }

    // Wait for all threads to finish.
    for (thread_n = 0; thread_n < N; thread_n++) {
      err = pthread_join(threads[thread_n], &retval);
      if (err != 0) {
        mexErrMsgTxt("Error when joining a thread!\n");

      }
    }

    // Free memory.
    if (D) {
      free((void*) D);
    }

    if (threads) {
      free((void*) threads);
    }

  } else { // Do not use threads.

      for (n=0; n<A_N; n++) {

        add_fftconv( A,n, A_M, B, B_M, &Y[0+n*(A_M+B_M-1)], A_L,
                     a,b,c,
                     reinterpret_cast<std::complex<double>*>(af),
                     reinterpret_cast<std::complex<double>*>(bf),
                     reinterpret_cast<std::complex<double>*>(cf));

        if (running==false) {
          printf("sum_fftconv_p: bailing out!\n");
          break;
        }
      }
      /***
          if (B_N > 1) {// B is a matrix.

          for (n=0; n<A_N; n++) {

          add_fftconv( A,n, A_M, &B[0+n*B_M], B_M, &Y[0+n*(A_M+B_M-1)], A_L,
          a,b,c,af,bf,cf);

          if (running==false) {
          printf("sum_fftconv_p: bailing out!\n");
          break;
          }

          } // end-for
          } else { // B is a vector.

          for (n=0; n<A_N; n++) {

          add_fftconv( A, n, A_M, B, B_M, &Y[0+n*(A_M+B_M-1)], A_L,
          a,b,c,af,bf,cf);

          if (running==false) {
          printf("sum_fftconv_p: bailing out!\n");
          break;
          }

          } // end-for

          } // end.if
      ***/
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
    mexErrMsgTxt("CTRL-C pressed!\n"); // Bail out.
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
  fftw_cleanup();

  return;
}
