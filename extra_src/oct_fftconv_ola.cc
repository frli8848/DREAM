/***
*
* Copyright (C) 2010,2011,2012,2014,2015,2016 Fredrik Lingvall
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

// $Revision: 909 $ $Date: 2016-11-25 13:13:43 +0100 (Fri, 25 Nov 2016) $ $LastChangedBy: frli8848 $

#include <string.h>
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

#define DOUBLE_PRECISION 0
#define SINGLE_PRECISION 1

#define EQU 0
#define SUM 1
#define NEG 2
#define OLA 3

// TODO: Add check if input pars have zero length
// TODO: Add help text for in-place mode.
// We pobably need fftw >= 3.2.x to handle 64-bit array indexing.
// TODO: If we have fftw 3.2.x we can use fftw_plan_dft_r2c_1d_64 etc.

#include "dream.h"

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
int mode = OLA;
int plan_method = 4; // Default to ESTIMATE method.

// FFTW plans.
fftw_plan    p_forward;
fftw_plan    p_backward;
fftwf_plan  sp_forward;
fftwf_plan  sp_backward;

//
// typedef:s
//

typedef struct
{
  dream_idx_type line_start;
  dream_idx_type line_stop;
  const double *A;
  dream_idx_type A_M;
  dream_idx_type A_N;
  const double *B;
  dream_idx_type B_M;
  dream_idx_type B_N;
  dream_idx_type block_len;
  double *Y;
} DATA;

typedef struct
{
  dream_idx_type line_start;
  dream_idx_type line_stop;
  const float *A;
  dream_idx_type A_M;
  dream_idx_type A_N;
  const float *B;
  dream_idx_type B_M;
  dream_idx_type B_N;
  dream_idx_type block_len;
  float *Y;
} S_DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//

void* smp_process(void *arg);
void* smp_s_process(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

void fftconv(const double *xr, dream_idx_type nx, const double *yr, dream_idx_type ny, double *zr,
             double      *a,  double *b, double *c,
             std::complex<double> *af, std::complex<double> *bf, std::complex<double> *cf);

void sfftconv(const float *xr, dream_idx_type nx, const float *yr, dream_idx_type ny, float *zr,
             float      *a,  float *b, float *c,
             std::complex<float> *af, std::complex<float> *bf, std::complex<float> *cf);

/***
 *
 * Thread function (double precision).
 *
 ***/

void* smp_process(void *arg)
{
  DATA D = *(DATA *)arg;
  dream_idx_type    line_start=D.line_start, line_stop=D.line_stop, n;
  const double *A = D.A, *B = D.B;
  double *Y = D.Y;
  dream_idx_type A_M = D.A_M, B_M = D.B_M, B_N = D.B_N;
  dream_idx_type l,k, block_len = D.block_len;

  // Input vectors.
  double *a  = NULL, *b  = NULL;
  double *c  = NULL;

  // Fourier Coefficients.
  std::complex<double> *af  = NULL, *bf  = NULL, *cf  = NULL;
  dream_idx_type fft_len;

  // Allocate space for input vectors.

  //fft_len = A_M+B_M-1;
  fft_len = block_len+B_M-1;

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

  double *y = c; // We use the same data space here.
  double *a_block = (double*) fftw_malloc(block_len*sizeof(double));

  // First clear the output vector since the overlap-and-add
  // method don't overwrite data - it just adds results to the existing array.
  memset(Y,0,(line_stop-line_start+1)*(A_M+B_M-1)*sizeof(double));

  //
  // Do the convolution.
  //

  if (B_N > 1) {// B is a matrix.

    for (n=line_start; n<line_stop; n++) {

      //
      // The overlap-and-add algorithm.
      //

      k = 0;
      while (k < A_M ) {

        // Convolve one block
        if (k+block_len < A_M)

          fftconv( &A[k+n*A_M], block_len, &B[0+n*B_M], B_M, y, a,b,c,af,bf,cf);

        else {

          // We must do the convolution (eg. the FFT) of the same lenght here (=fft_len) but the first
          // arg is now shorter than above so copy and zero-pad the first arg to fftconv.

          // Copy.
          //for (l=0; l<A_M-k; l++)
          //  a_block[l] = A[k+l + n*A_M];
          memcpy(a_block,&A[k + n*A_M],(A_M-k)*sizeof(double));

          // Zero-pad.
          //for (l=A_M-k; l<block_len; l++)
          //  a_block[l] = 0.0; // Zero-pad.
          memset(&a_block[A_M-k],0,(block_len-(A_M-k))*sizeof(double));

          fftconv( a_block, block_len, &B[0+n*B_M], B_M, y, a,b,c,af,bf,cf);
        }

        // Add the overlap.
        if (k+fft_len <  A_M+B_M-1) {
          for (l=0; l<fft_len; l++)
            Y[k+l + n*(A_M+B_M-1)] += y[l];
        } else {
          for (l=0; l<(A_M+B_M-1)-k; l++)
            Y[k+l + n*(A_M+B_M-1)] += y[l];
        }

        k += block_len;
      } // while

      if (running==false) {
        octave_stdout << "fftconv_ola: thread for column " << line_start+1 << " -> " <<
          line_stop << " bailing out at column = " << n <<"!\n";
        break;
      }

    } // end-for
  } else { // B is a vector.

    for (n=line_start; n<line_stop; n++) {

      //
      // The overlap-and-add algorithm.
      //

      k = 0;
      while (k < A_M ) {

        // Convolve one block
        if (k+block_len < A_M)

          fftconv( &A[k+n*A_M], block_len, B, B_M, y, a,b,c,af,bf,cf);

        else {

          // We must do the convolution (eg. the FFT) of the same lenght here (=fft_len) but the first
          // arg is now shorter than above so copy and zero-pad the first arg to fftconv.

          // Copy.
          //for (l=0; l<(A_M-k); l++)
          //  a_block[l] = A[k+l + n*A_M];
          memcpy(a_block,&A[k + n*A_M],(A_M-k)*sizeof(double));

          // Zero-pad.
          //for (l=A_M-k; l<block_len; l++)
          //  a_block[l] = 0.0; // Zero-pad.
          memset(&a_block[A_M-k],0,(block_len-(A_M-k))*sizeof(double));

          fftconv( a_block, block_len, B, B_M, y, a,b,c,af,bf,cf);

        }

        // Add the overlap.
        if (k+fft_len < A_M+B_M-1) {

          for (l=0; l<fft_len; l++)
            Y[k+l + n*(A_M+B_M-1)] += y[l];

        } else {

          for (l=0; l<(A_M+B_M-1)-k; l++)
            Y[k+l + n*(A_M+B_M-1)] += y[l];

        }

        k += block_len;
      } // while

      if (running==false) {
        octave_stdout << "fftconv_ola: thread for column " << line_start+1 << " -> " <<
          line_stop << " bailing out at column = " << n <<"!\n";
        break;
      }

    } // for
  } // if

  //
  //  Cleanup
  //

  // Free buffer memory.
  if (a_block)
    fftw_free(a_block);
  else
    error("smp_process a_block memory free failed in fftconv_ola thread!!");

  if (a)
    fftw_free(a);
  else {
    error("smp_process a memory free failed in fftconv_ola thread!!");
  }

  if(af)
    fftw_free(af);
  else {
    error("smp_process af memory free failed in fftconv_ola thread!!");
  }

  if(b)
    fftw_free(b);
  else {
    error("smp_process b memory free failed in fftconv_ola thread!!");
  }

  if (bf)
    fftw_free(bf);
  else {
    error("smp_process bf memory free failed in fftconv_ola thread!!");
  }

  if (c)
    fftw_free(c);
  else {
    error("smp_process c memory free failed in fftconv_ola thread!!");
  }

  if (cf)
    fftw_free(cf);
  else {
    error("smp_process cf memory free failed in fftconv_ola thread!!");
  }

  return(NULL);
}

/***
 *
 * Thread function (single precision).
 *
 ***/

void* smp_s_process(void *arg)
{
  S_DATA D = *(S_DATA *)arg;
  dream_idx_type    line_start=D.line_start, line_stop=D.line_stop, n;
  const float *A = D.A, *B = D.B;
  float *Y = D.Y;
  dream_idx_type A_M = D.A_M, B_M = D.B_M, B_N = D.B_N;
  dream_idx_type l,k, block_len = D.block_len;

  // Input vectors.
  float *a  = NULL, *b  = NULL;
  float *c  = NULL;

  // Fourier Coefficients.
  std::complex<float> *af  = NULL, *bf  = NULL, *cf  = NULL;
  dream_idx_type fft_len;

  // Allocate space for input vectors.

  //fft_len = A_M+B_M-1;
  fft_len = block_len+B_M-1;

  //
  // fftw_malloc may not be thread safe! Check this!!!
  //

  a = (float*) fftw_malloc(sizeof(float)*2*(fft_len/2+1));
  //af = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));
  af = (std::complex<float>*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));

  b = (float*) fftw_malloc(sizeof(float)*2*(fft_len/2+1));
  bf = (std::complex<float>*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));

  // Allocate space for output vector.

  c = (float*) fftw_malloc(sizeof(float)*2*(fft_len/2+1));
  //cf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));
  cf = (std::complex<float>*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));

  float *y = c; // We use the same data space here.
  float *a_block = (float*) fftw_malloc(block_len*sizeof(float));

  // First clear the output vector since the overlap-and-add
  // method don't overwrite data - it just adds results to the existing array.
  memset(Y,0,(line_stop-line_start+1)*(A_M+B_M-1)*sizeof(float));

  //
  // Do the convolution.
  //

  if (B_N > 1) {// B is a matrix.

    for (n=line_start; n<line_stop; n++) {

      //
      // The overlap-and-add algorithm.
      //

      k = 0;
      while (k < A_M ) {

        // Convolve one block
        if (k+block_len < A_M)

          sfftconv( &A[k+n*A_M], block_len, &B[0+n*B_M], B_M, y, a,b,c,af,bf,cf);

        else {

          // We must do the convolution (eg. the FFT) of the same lenght here (=fft_len) but the first
          // arg is now shorter than above so copy and zero-pad the first arg to fftconv.

          // Copy.
          //for (l=0; l<A_M-k; l++)
          //  a_block[l] = A[k+l + n*A_M];
          memcpy(a_block,&A[k + n*A_M],(A_M-k)*sizeof(float));

          // Zero-pad.
          //for (l=A_M-k; l<block_len; l++)
          //  a_block[l] = 0.0; // Zero-pad.
          memset(&a_block[A_M-k],0,(block_len-(A_M-k))*sizeof(float));

          sfftconv( a_block, block_len, &B[0+n*B_M], B_M, y, a,b,c,af,bf,cf);
        }

        // Add the overlap.
        if (k+fft_len <  A_M+B_M-1) {
          for (l=0; l<fft_len; l++)
            Y[k+l + n*(A_M+B_M-1)] += y[l];
        } else {
          for (l=0; l<(A_M+B_M-1)-k; l++)
            Y[k+l + n*(A_M+B_M-1)] += y[l];
        }

        k += block_len;
      } // while

      if (running==false) {
        octave_stdout << "fftconv_ola: thread for column " << line_start+1 << " -> " <<
          line_stop << " bailing out at column = " << n <<"!\n";
        break;
      }

    } // end-for
  } else { // B is a vector.

    for (n=line_start; n<line_stop; n++) {

      //
      // The overlap-and-add algorithm.
      //

      k = 0;
      while (k < A_M ) {

        // Convolve one block
        if (k+block_len < A_M)

          sfftconv( &A[k+n*A_M], block_len, B, B_M, y, a,b,c,af,bf,cf);

        else {

          // We must do the convolution (eg. the FFT) of the same lenght here (=fft_len) but the first
          // arg is now shorter than above so copy and zero-pad the first arg to fftconv.

          // Copy.
          //for (l=0; l<(A_M-k); l++)
          //  a_block[l] = A[k+l + n*A_M];
          memcpy(a_block,&A[k + n*A_M],(A_M-k)*sizeof(float));

          // Zero-pad.
          //for (l=A_M-k; l<block_len; l++)
          //  a_block[l] = 0.0; // Zero-pad.
          memset(&a_block[A_M-k],0,(block_len-(A_M-k))*sizeof(float));

          sfftconv( a_block, block_len, B, B_M, y, a,b,c,af,bf,cf);

        }

        // Add the overlap.
        if (k+fft_len < A_M+B_M-1) {

          for (l=0; l<fft_len; l++)
            Y[k+l + n*(A_M+B_M-1)] += y[l];

        } else {

          for (l=0; l<(A_M+B_M-1)-k; l++)
            Y[k+l + n*(A_M+B_M-1)] += y[l];

        }

        k += block_len;
      } // while

      if (running==false) {
        octave_stdout << "fftconv_ola: thread for column " << line_start+1 << " -> " <<
          line_stop << " bailing out at column = " << n <<"!\n";
        break;
      }

    } // for
  } // if

  //
  //  Cleanup
  //

  // Free buffer memory.
  if (a_block)
    fftw_free(a_block);
  else
    error("smp_s_process a_block memory free failed in fftconv_ola thread!!");

  if (a)
    fftw_free(a);
  else {
    error("smp_s_process a memory free failed in fftconv_ola thread!!");
  }

  if(af)
    fftw_free(af);
  else {
    error("smp_s_process af memory free failed in fftconv_ola thread!!");
  }

  if(b)
    fftw_free(b);
  else {
    error("smp_s_process b memory free failed in fftconv_ola thread!!");
  }

  if (bf)
    fftw_free(bf);
  else {
    error("smp_s_process bf memory free failed in fftconv_ola thread!!");
  }

  if (c)
    fftw_free(c);
  else {
    error("smp_s_process c memory free failed in fftconv_ola thread!!");
  }

  if (cf)
    fftw_free(cf);
  else {
    error("smp_s_process a_block memory free failed in fftconv_ola thread!!");
  }

  return(NULL);
}


/***
 *
 * Convolution of two (double precision) vectors.
 *
 ***/

void fftconv(const double *xr, dream_idx_type nx, const double *yr, dream_idx_type ny, double *zr,
             double      *a,  double *b, double *c,
             std::complex<double> *af, std::complex<double> *bf, std::complex<double> *cf)
{
  dream_idx_type n, fft_len;

  fft_len = nx+ny-1;

  //
  // Copy and zero-pad.
  //

  //for (n=0; n < nx; n++)
  //  a[n] = xr[n];
  memcpy(a,xr,nx*sizeof(double));
  //for (n=nx; n < fft_len; n++)
  //  a[n] = 0.0; // Zero-pad.
  memset(&a[nx],0,(fft_len-nx+1)*sizeof(double));

  //for (n=0; n < ny; n++)
  //  b[n] = yr[n];
  memcpy(b,yr,ny*sizeof(double));
  //for (n=ny; n < fft_len; n++)
  //  b[n] = 0.0; // Zero-pad.
  memset(&b[ny],0,(fft_len-ny+1)*sizeof(double));

  //
  // Do the forward FFT transforms.
  //

  // Fourier transform xr.
  fftw_execute_dft_r2c(p_forward,a,reinterpret_cast<fftw_complex*>(af));

  // Fourier transform yr.
  fftw_execute_dft_r2c(p_forward,b,reinterpret_cast<fftw_complex*>(bf));

  //
  // Do the filtering.
  //

  for (n = 0; n < fft_len; n++) {
    cf[n] = (af[n] * bf[n])  / ((double) (fft_len));
  }

  //
  // Compute the inverse DFT of the filtered data.
  //

  fftw_execute_dft_c2r(p_backward,reinterpret_cast<fftw_complex*>(cf),c);

  if (mode != OLA) {

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

        // case OLA:
        // Do nothing since the copy is made in the overlap update.
        //break;

      default:
        // in-place '=' operation.
        //for (n = 0; n < fft_len; n++) {
        //zr[n] = c[n];
        //}
        memcpy(zr,c,fft_len*sizeof(double));
        break;
      }
    }
  } // if (mode != OLA)
}


/***
 *
 * Convolution of two single precision vectors.
 *
 ***/

void sfftconv(const float *xr, dream_idx_type nx, const float *yr, dream_idx_type ny, float *zr,
             float      *a,  float *b, float *c,
             std::complex<float> *af, std::complex<float> *bf, std::complex<float> *cf)
{
  dream_idx_type n, fft_len;

  fft_len = nx+ny-1;

  //
  // Copy and zero-pad.
  //

  //for (n=0; n < nx; n++)
  //  a[n] = xr[n];
  memcpy(a,xr,nx*sizeof(float));
  //for (n=nx; n < fft_len; n++)
  //  a[n] = 0.0; // Zero-pad.
  memset(&a[nx],0,(fft_len-nx+1)*sizeof(float));

  //for (n=0; n < ny; n++)
  //  b[n] = yr[n];
  memcpy(b,yr,ny*sizeof(float));
  //for (n=ny; n < fft_len; n++)
  //  b[n] = 0.0; // Zero-pad.
  memset(&b[ny],0,(fft_len-ny+1)*sizeof(float));

  //
  // Do the forward FFT transforms.
  //

  // Fourier transform xr (single precision).
  fftwf_execute_dft_r2c(sp_forward,a,reinterpret_cast<fftwf_complex*>(af));

  // Fourier transform yr (single precision).
  fftwf_execute_dft_r2c(sp_forward,b,reinterpret_cast<fftwf_complex*>(bf));

  //
  // Do the filtering.
  //

  for (n = 0; n < fft_len; n++) {
    cf[n] = (af[n] * bf[n])  / ((float) (fft_len));
  }

  //
  // Compute the inverse (single precision) DFT of the filtered data.
  //

  fftwf_execute_dft_c2r(sp_backward,reinterpret_cast<fftwf_complex*>(cf),c);

  if (mode != OLA) {

    // Copy data to output matrix.
    if (in_place == false) {

      //for (n = 0; n < fft_len; n++)
      //  zr[n] = c[n];
      memcpy(zr,c,fft_len*sizeof(float));

    } else { // in-place

      switch (mode) {

      case EQU:
        // in-place '=' operation.
        //for (n = 0; n < fft_len; n++) {
        //zr[n] = c[n];
        //}
        memcpy(zr,c,fft_len*sizeof(float));
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

        // case OLA:
        // Do nothing since the copy is made in the overlap update.
        //break;

      default:
        // in-place '=' operation.
        //for (n = 0; n < fft_len; n++) {
        //zr[n] = c[n];
        //}
        memcpy(zr,c,fft_len*sizeof(float));
        break;
      }
    }
  } // if (mode != OLA)
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
 * Octave (oct) gateway function for FFTCONV_OLA.
 *
 ***/

DEFUN_DLD (fftconv_ola, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} [Y,wisdom_str_out] = fftconv_ola(A,B,block_len,wisdom_str_in);\n\
\n\
FFTCONV_OLA - Computes the one dimensional convolution of the columns in matrix A and the matrix (or vector) B \n\
using the overlap-and-add method.\n\
\n\
Input parameters:\n\
\n\
@table @code\n\
@item A\n\
An M x N matrix.\n\
@item B\n\
A K x N matrix or a K-length vector. If B is a vector each column in A is convolved with the vector B.\n\
@item wisdom_str_in\n\
Optional parameter. If the wisdom_str_in parameter is not supplied then @code{fftconv_ola} calls fftw wisdom plan\n\
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
Optional parameter. If the wisdom_str_out output parameter is supplied then @code{fftconv_ola}\n\
will call fftw wisdom plan functions and return the wisdom string which then can be used to\n\
speed up subsequent calls to @code{fftconv_ola} by suppying the string as the input argument\n\
wisdom_str_in.\n\
@end table\n\
\n\
In place modes:\n\
\n\
fftconv_ola(A,B,block_len,Y);\n\
\n\
fftconv_ola(A,B,block_len,Y,mode);\n\
\n\
fftconv_ola(A,B,block_len,Y,wisdom_str_in);\n\
\n\
fftconv_ola(A,B,block_len,Y,mode,wisdom_str_in);\n\
\n\
where the 'mode' is a string which can be '=', '+=', or '-='.\n\
\n\
NOTE: fftconv_ola requires the FFTW library version 3 @url{http://www.fftw.org}.\n\
\n\
fftconv_ola is a part of the DREAM Toolbox available at\n\
@url{http://www.signal.uu.se/Toolbox/dream/}.\n\
\n\
Copyright @copyright{} 2010-2016 Fredrik Lingvall.\n\
@seealso {conv_p, fftconv_p, fftconv, conv, fftw_wisdom}\n\
@end deftypefn")
{
  int data_format;

  const double *A=NULL,*B=NULL;	// Double precision.
  double *Y=NULL, *a_block=NULL;
  const float *sA=NULL,*sB=NULL; // Single precision.
  float *sY=NULL, *s_a_block=NULL;
  dream_idx_type line_start, line_stop, A_M=0, A_N=0, B_M=0, B_N=0, n;
  sighandler_t   old_handler, old_handler_abrt, old_handler_keyint;
  std::thread     *threads;
  octave_idx_type  thread_n, nthreads;
  DATA    *D=NULL;		// Double precision.
  S_DATA *sD=NULL;		// Single precision.

  // Input vectors (only used for creating fftw plans).
  double *a=NULL, *b=NULL;	// Double precision.
  double *c=NULL;
  float *sa=NULL, *sb=NULL;	// Single precision.
  float *sc=NULL;

  // Fourier Coefficients (only used for creating fftw plans).
  //fftw_complex *af=NULL, *bf=NULL, *cf=NULL;
  //fftwf_complex *af=NULL, *bf=NULL, *cf=NULL;
  std::complex<double> *af=NULL, *bf=NULL, *cf=NULL;
  std::complex<float> *saf=NULL, *sbf=NULL, *scf=NULL;

  dream_idx_type fft_len, return_wisdom = false, load_wisdom = false;
  char *the_str = NULL;
  int buflen, is_set = false;
  dream_idx_type l, k, block_len;
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

  mode = OLA;

  // Check for proper inputs arguments.
  switch (nrhs) {

  case 0:
  case 1:			// A matrix/vector
  case 2:			// B matrix/vector
    error("fftconv_ola requires 3 to 6 input arguments!");
    return oct_retval;
    break;

  case 3:			// OLA bleock length

    if (nlhs > 2) {
      error("Too many output arguments for fftconv_ola!");
      return oct_retval;
    }

    if (nlhs == 2)
      return_wisdom = true;

    break;

  case 4: // Output arg in in-place mode or a wisdom string in normal mode.
    if ( args(3).is_string() ) { // 4th arg is a fftw wisdom string.

      std::string strin = args(3).string_value();
      buflen = strin.length();
      the_str = (char*) fftw_malloc(buflen * sizeof(char));
      for ( n=0; n<buflen; n++ ) {
        the_str[n] = strin[n];
      }

      //
      // If 4:th arg is a string then only a wisdom string is valid.
      //

      if (!strcmp("fftw_wisdom",the_str) && (!strcmp("fftwf_wisdom",the_str))) {
        error("The string in arg 4 do not seem to be in a fftw wisdom format!");
        return oct_retval;
      }
      else
        load_wisdom = true;


    } else { // 4th arg not a string then assume in-place mode.
      fftw_forget_wisdom(); // Clear wisdom history (a new wisdom will be created below).
      if (nlhs > 0) {
        error("4th arg is not a fftw wisdom string and in-place mode is assumed. But then there should be no output args!");
        return oct_retval;
      }
    }
    break;

  case 5:  // In-place mode if >= 5 args.
    if ( args(4).is_string() ) { // 5th arg is a string ('=','+=','-=', or fftw wisdom).

      std::string strin = args(4).string_value();
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
        if ( (strcmp("fftw_wisdom",the_str) < 0) && (strcmp("fftwf_wisdom",the_str) < 0) ) {
          error("Non-valid string in arg 5!");
          return oct_retval;
        } else {
          load_wisdom = true;
        }
      }
    } else { // 5th arg not a string
      error("Argument 5 is not a valid string format!");
      return oct_retval;
    }
    break;


  case 6: // In-place mode if 6 input args.
    if ( args(5).is_string() ) { // 6:th arg is a string (fftw wisdom).

      // Read the wisdom string.
      std::string strin = args(5).string_value();
      buflen = strin.length();
      the_str = (char*) fftw_malloc(buflen * sizeof(char));
      for ( n=0; n<buflen; n++ ) {
        the_str[n] = strin[n];
      }

      if ( (strcmp("fftw_wisdom",the_str) < 0) && (strcmp("fftwf_wisdom",the_str) < 0) ) {
        error("The string in 6th arg do not seem to be in a fftw wisdom format!");
        return oct_retval;
      } else
        load_wisdom = true;

    } else { // 6th arg not a string
      error("Argument 6 is not a valid string format!");
      return oct_retval;
    }
    break;

  default:
    error("fftconv_ola requires 3 to 6 input arguments!");
    return oct_retval;
    break;
  }

  //
  // Check data format
  //

  // Double precision input data.
  if(args(0).is_double_type() && args(1).is_double_type() ) {

    data_format = DOUBLE_PRECISION;

    const Matrix tmp0 = args(0).matrix_value();
    const Matrix tmp1 = args(1).matrix_value();

    A_M = tmp0.rows();
    A_N = tmp0.cols();

    B_M = tmp1.rows();
    B_N = tmp1.cols();
  }

  // Single precision input data.
  if(args(0).is_single_type() && args(1).is_single_type() ) {

    data_format = SINGLE_PRECISION;

    const FloatMatrix tmp0 = args(0).float_matrix_value();
    const FloatMatrix tmp1 = args(1).float_matrix_value();

    A_M = tmp0.rows();
    A_N = tmp0.cols();

    B_M = tmp1.rows();
    B_N = tmp1.cols();

    //A = tmp0.data();
    //B = tmp1.data();
  }

  // Check if arg 1 and 2 have different precisions.
  if( (args(0).is_double_type() && args(1).is_single_type()) ||
      (args(0).is_single_type() && args(1).is_double_type()) ) {
    error("Input arg 1 and 2 don't have the same precision!");
    return oct_retval;
  }

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
  // The segment length (block length).
  //

  if (mxGetM(2)*mxGetN(2) ==1)
    block_len = (dream_idx_type) args(2).matrix_value().data()[0];
  else {
    error("Argument 3 must be a scalar!");
    return oct_retval;
  }

  if (block_len < 0 || block_len > A_M) {
    error("Argument 3 is out of bounds! Must be > 0 and less than number of rows of arg 1!");
    return oct_retval;
  }

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11)
  nthreads = std::thread::hardware_concurrency();

  // Check nthreads argument.
  if (nthreads > A_N) {
    // Add a -v verbose flag to display stuff like this!
    //mexPrintf("Warning: nthreads is larger then number of columns in first arg.\n");
    //mexPrintf("         Setting nthreads = # cols in 1st arg!\n");
    nthreads = A_N;
  }

  //
  // Register signal handlers.
  //

  if ((old_handler = signal(SIGTERM, &sighandler)) == SIG_ERR) {
    printf("Couldn't register signal handler.\n");
  }

  if ((old_handler_abrt = signal(SIGABRT, &sighandler)) == SIG_ERR) {
    printf("Couldn't register signal handler.\n");
  }

  if ((old_handler_keyint = signal(SIGINT, &sighandler)) == SIG_ERR) {
    printf("Couldn't register signal handler.\n");
  }


  // ********************************************
  //
  // Init the FFTW plans.
  //
  // ********************************************

  if(load_wisdom) {

    if (data_format == DOUBLE_PRECISION) {
      if (!fftw_import_wisdom_from_string(the_str)) {
        error("Failed to load (double precision) fftw wisdom!");
        return oct_retval;
      } else
        fftw_free(the_str); // Clean up.
    }

    if (data_format == SINGLE_PRECISION) {
      if (!fftwf_import_wisdom_from_string(the_str)) {
        error("Failed to load (double precision) fftw wisdom!");
        return oct_retval;
      } else
        fftw_free(the_str); // Clean up.
    }

  }

  //
  // Allocate space for (temp) input/output vectors (only used for creating plans).
  //

  //fft_len = A_M+B_M-1;
  fft_len = block_len+B_M-1;

  if (data_format == DOUBLE_PRECISION) {

    a = (double*) fftw_malloc(sizeof(double)*2*(fft_len/2+1));
    //af = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));
    af = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));

    b = (double*) fftw_malloc(sizeof(double)*2*(fft_len/2+1));
    //bf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));
    bf = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));

    c = (double*) fftw_malloc(sizeof(double)*2*(fft_len/2+1));
    //cf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));
    cf = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*2*(fft_len/2+1));


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

  } //  if (data_format == DOUBLE_PRECISION)

  if (data_format == SINGLE_PRECISION) {

    sa = (float*) fftw_malloc(sizeof(float)*2*(fft_len/2+1));
    saf = (std::complex<float>*) fftw_malloc(sizeof(fftwf_complex)*2*(fft_len/2+1));

    sb = (float*) fftw_malloc(sizeof(float)*2*(fft_len/2+1));
    sbf = (std::complex<float>*) fftw_malloc(sizeof(fftwf_complex)*2*(fft_len/2+1));

    sc = (float*) fftw_malloc(sizeof(float)*2*(fft_len/2+1));
    scf = (std::complex<float>*) fftw_malloc(sizeof(fftwf_complex)*2*(fft_len/2+1));

    // 1) Very slow
    if (plan_method == 1) {
      sp_forward = fftwf_plan_dft_r2c_1d(fft_len, sa, reinterpret_cast<fftwf_complex*>(saf), FFTW_EXHAUSTIVE);
      sp_backward = fftwf_plan_dft_c2r_1d(fft_len, reinterpret_cast<fftwf_complex*>(scf), sc, FFTW_EXHAUSTIVE);
    }

    // 2) Slow
    if (plan_method == 2) {
      sp_forward = fftwf_plan_dft_r2c_1d(fft_len, sa, reinterpret_cast<fftwf_complex*>(saf), FFTW_PATIENT);
      sp_backward = fftwf_plan_dft_c2r_1d(fft_len, reinterpret_cast<fftwf_complex*>(scf), sc, FFTW_PATIENT);
    }

    // 3) Too slow on long FFTs.
    if (plan_method == 3) {
      sp_forward = fftwf_plan_dft_r2c_1d(fft_len, sa, reinterpret_cast<fftwf_complex*>(saf), FFTW_MEASURE);
      sp_backward = fftwf_plan_dft_c2r_1d(fft_len, reinterpret_cast<fftwf_complex*>(scf), sc, FFTW_MEASURE);
    }

    // 4)
    if (plan_method == 4) {
      sp_forward = fftwf_plan_dft_r2c_1d(fft_len, sa, reinterpret_cast<fftwf_complex*>(saf),  FFTW_ESTIMATE);
      sp_backward = fftwf_plan_dft_c2r_1d(fft_len, reinterpret_cast<fftwf_complex*>(scf), sc,  FFTW_ESTIMATE);
    }

  } // if (data_format == SINGLE_PRECISION)


  // ********************************************
  //
  // Double Precision
  //
  // ********************************************

  if (data_format == DOUBLE_PRECISION) {

    // Get pointers to input data.
    const Matrix tmp0 = args(0).matrix_value();
    const Matrix tmp1 = args(1).matrix_value();
    A = tmp0.data();
    B = tmp1.data();

    if (nrhs == 3 || (nrhs == 4 && load_wisdom)) { // Normal mode.

      //
      //
      // Normal (non in-place) mode.
      //
      //

      in_place = false;

      running = true;

      // Allocate mamory for the output matrix.
      Matrix Ymat(A_M+B_M-1, A_N);
      Y = Ymat.fortran_vec();

      if (nthreads>1) { // Use threads

        // Allocate mem for the threads.
        threads = new std::thread[nthreads]; // Init thread data.
        if (!threads) {
          error("Failed to allocate memory for threads!");
          return oct_retval;
        }

        // Allocate local data.
        D = (DATA*) malloc(nthreads*sizeof(DATA));
        if (!D) {
          error("Failed to allocate memory for thread data!");
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
          D[thread_n].block_len = block_len;
          D[thread_n].Y = Y;

          // Start the threads.
          threads[thread_n] = std::thread(smp_process, &D[thread_n]);

        } // for (thread_n = 0; thread_n < nthreads; thread_n++)

        //
        // Wait for all threads to finish.
        //

        for (thread_n = 0; thread_n < nthreads; thread_n++)
          threads[thread_n].join();

        // Free memory.
        if (D)
          free((void*) D);

      } else { // Do not use threads.

        double *y = c; // We use the same data space here.
        a_block = (double*) fftw_malloc(block_len*sizeof(double));

        if (B_N > 1) {// B is a matrix.

          for (n=0; n<A_N; n++) {

            //
            // The overlap-and-add algorithm.
            //

            k = 0;
            while (k < A_M ) {

              // Convolve one block
              if (k+block_len < A_M)

                fftconv( &A[k+n*A_M], block_len, &B[0+n*B_M], B_M, y, a,b,c,af,bf,cf);

              else {

                // We must do the convolution (eg. the FFT) of the same lenght here (=fft_len) but the first
                // arg is now shorter than above so copy and zero-pad the first arg to fftconv.

                // Copy.
                //for (l=0; l<(A_M-k); l++)
                // a_block[l] = A[k+l + n*A_M];
                memcpy(a_block,&A[k + n*A_M],(A_M-k)*sizeof(double));

                // Zero-pad.
                //for (l=A_M-k; l<block_len; l++)
                // a_block[l] = 0.0; // Zero-pad.
                memset(&a_block[A_M-k],0,(block_len-(A_M-k))*sizeof(double));

                fftconv( a_block, block_len, &B[0+n*B_M], B_M, y, a,b,c,af,bf,cf);
              }

              // Add the overlap.
              if (k+fft_len <  A_M+B_M-1) {
                for (l=0; l<fft_len; l++)
                  Y[k+l + n*(A_M+B_M-1)] += y[l];
              } else {
                for (l=0; l<(A_M+B_M-1)-k; l++)
                  Y[k+l + n*(A_M+B_M-1)] += y[l];
              }

              k += block_len;
            } // while

            if (running==false) {
              printf("fftconv_ola: bailing out!\n");
              break;
            }

          } // end-for

        } else { // B is a vector.

          for (n=0; n<A_N; n++) {

            //
            // The overlap-and-add algorithm.
            //

            k = 0;
            while (k < A_M ) {

              // Convolve one block
              if (k+block_len < A_M)

                fftconv( &A[k+n*A_M], block_len, B, B_M, y, a,b,c,af,bf,cf);

              else {

                // We must do the convolution (eg. the FFT) of the same lenght here (=fft_len) but the first
                // arg is now shorter than above so copy and zero-pad the first arg to fftconv.

                // Copy.
                //for (l=0; l<(k+block_len)-A_M; l++)
                // a_block[l] = A[k+l + n*A_M];
                memcpy(a_block,&A[k + n*A_M],(A_M-k)*sizeof(double));

                // Zero-pad.
                //for (l=A_M-k; l<block_len; l++)
                // a_block[l] = 0.0; // Zero-pad.
                memset(&a_block[A_M-k],0,(block_len-(A_M-k))*sizeof(double));

                fftconv( a_block, block_len, B, B_M, y, a,b,c,af,bf,cf);
              }

              // Add the overlap.
              if (k+fft_len <  A_M+B_M-1) {
                for (l=0; l<fft_len; l++)
                  Y[k+l + n*(A_M+B_M-1)] += y[l];
              } else {
                for (l=0; l<(A_M+B_M-1)-k; l++)
                  Y[k+l + n*(A_M+B_M-1)] += y[l];
              }

              k += block_len;
            } // while

            if (running==false) {
              printf("fftconv_ola: bailing out!\n");
              break;
            }

          } // end-for
        } // end-if

        // Free buffer memory.
        if (a_block)
          fftw_free(a_block);
        else
          error("a_block memory free in fftconv_ola failed!");

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
      if (a)
        fftw_free(a);
      else
        error("a memory free failed in fftconv_ola!");

      if (af)
        fftw_free(af);
      else
        error("af memory free failed in fftconv_ola!");

      if (b)
        fftw_free(b);
      else
        error("b memory free failed in fftconv_ola!");

      if (bf)
        fftw_free(bf);
      else
        error("bf memory free failed in fftconv_ola!");

      if (c)
        fftw_free(c);
      else
        error("c memory free failed in fftconv_ola!");

      if (cf)
        fftw_free(cf);
      else
        error("cf memory free failed in fftconv_ola!");

      // Cleanup the plans.
      fftw_destroy_plan(p_forward);
      fftw_destroy_plan(p_backward);

      // Clean up FFTW
      //fftw_cleanup(); This seems to put Octave in an unstable state. Calling fftconv will crash Octave.

      return oct_retval;

    }


    if ( (nrhs == 4 && !load_wisdom) || nrhs == 5 || nrhs == 6) { // In-place mode.

      in_place = true;

      //
      //
      // In-place mode.
      //
      //

      if (args(3).matrix_value().rows() != A_M+B_M-1) {
        error("Wrong number of rows in argument 5!");
        return oct_retval;
      }

      if (args(3).matrix_value().cols() != A_N) {
        error("Wrong number of columns in argument 5!");
        return oct_retval;
      }

      // Check datatype of (in-place) output arg.
      if (!args(3).is_double_type()) {
        error("The output matrix (4th arg) must ba a double precision matrix when the input matricies are double precision!");
        return oct_retval; // FIXME : Do we need to clear memory here?
      }

      const Matrix Ytmp = args(3).matrix_value();
      Y = (double*) Ytmp.data();

      //
      // Call the CONV subroutine.
      //

      running = true;

      if (nthreads>1) { // Use threads

        // Allocate mem for the threads.
        threads = new std::thread[nthreads]; // Init thread data.
        if (!threads) {
          error("Failed to allocate memory for threads!");
          return oct_retval;
        }

        //
        // Double precision
        //


        // Allocate local data.
        D = (DATA*) malloc(nthreads*sizeof(DATA));
        if (!D) {
          error("Failed to allocate memory for thread data!");
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
          D[thread_n].block_len = block_len;
          D[thread_n].Y = Y;

          // Start the threads.
          threads[thread_n] = std::thread(smp_process, &D[thread_n]);

        } // for (thread_n = 0; thread_n < nthreads; thread_n++)

        // Wait for all threads to finish.
        for (thread_n = 0; thread_n < nthreads; thread_n++)
          threads[thread_n].join();

        // Free memory.
        if (D)
          free((void*) D);

      } else { // Do not use threads.

        double *y = c; // We use the same data space here.
        a_block = (double*) fftw_malloc(block_len*sizeof(double));

        if (B_N > 1) {// B is a matrix.

          for (n=0; n<A_N; n++) {

            //
            // The overlap-and-add algorithm.
            //

            k = 0;
            while (k < A_M ) {

              // Convolve one block
              if (k+block_len < A_M)

                fftconv( &A[k+n*A_M], block_len, &B[0+n*B_M], B_M, y, a,b,c,af,bf,cf);

              else {

                // We must do the convolution (eg. the FFT) of the same lenght here (=fft_len) but the first
                // arg is now shorter than above so copy and zero-pad the first arg to fftconv.

                // Copy.
                //for (l=0; l<(A_M-k); l++)
                //  a_block[l] = A[k+l + n*A_M];
                memcpy(a_block,&A[k + n*A_M],(A_M-k)*sizeof(double));

                // Zero-pad.
                //for (l=A_M-k; l<block_len; l++)
                // a_block[l] = 0.0; // Zero-pad.
                memset(&a_block[A_M-k],0,(block_len-(A_M-k))*sizeof(double));

                fftconv( a_block, block_len, &B[0+n*B_M], B_M, y, a,b,c,af,bf,cf);
              }

              // Add the overlap.
              if (k+fft_len <  A_M+B_M-1) {
                for (l=0; l<fft_len; l++)
                  Y[k+l + n*(A_M+B_M-1)] += y[l];
              } else {
                for (l=0; l<(A_M+B_M-1)-k; l++)
                  Y[k+l + n*(A_M+B_M-1)] += y[l];
              }

              k += block_len;
            } // while

            if (running==false) {
              printf("fftconv_ola: bailing out!\n");
              break;
            }

          } // end-for
        } else { // B is a vector.

          for (n=0; n<A_N; n++) {

            //
            // The overlap-and-add algorithm.
            //

            k = 0;
            while (k < A_M ) {

              // Convolve one block
              if (k+block_len < A_M)

                fftconv( &A[k+n*A_M], block_len, B, B_M, y, a,b,c,af,bf,cf);

              else {

                // We must do the convolution (eg. the FFT) of the same lenght here (=fft_len) but the first
                // arg is now shorter than above so copy and zero-pad the first arg to fftconv.

                // Copy.
                //for (l=0; l<(A_M-k); l++)
                // a_block[l] = A[k+l + n*A_M];
                memcpy(a_block,&A[k + n*A_M],(A_M-k)*sizeof(double));

                // Zero-pad.
                //for (l=0; l<(A_M-k); l++)
                // a_block[l] = 0.0; // Zero-pad.
                memset(&a_block[A_M-k],0,(block_len-(A_M-k))*sizeof(double));

                fftconv( a_block, block_len, B, B_M, y, a,b,c,af,bf,cf);
              }

              // Add the overlap.
              if (k+fft_len <  A_M+B_M-1) {
                for (l=0; l<fft_len; l++)
                  Y[k+l + n*(A_M+B_M-1)] += y[l];
              } else {
                for (l=0; l<(A_M+B_M-1)-k; l++)
                  Y[k+l + n*(A_M+B_M-1)] += y[l];
              }

              k += block_len;
            } // while

            if (running==false) {
              printf("fftconv_ola: bailing out!\n");
              break;
            }

          } // end-for
        } // end-if


        // Free buffer memory.
        if (a_block)
          fftw_free(a_block);
        else
          error("a_block memory free failed in fftconv_ola!");

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

      if(a)
        fftw_free(a);
      else
        error("a memory free failed in fftconv_ola!");

      if (af)
        fftw_free(af);
      else
        error("af memory free failed in fftconv_ola!");

      if (b)
        fftw_free(b);
      else
        error("b memory free failed in fftconv_ola!");

      if (bf)
        fftw_free(bf);
      else
        error("bf memory free failed in fftconv_ola!");

      if (c)
        fftw_free(c);
      else
        error("c memory free failed in fftconv_ola!");

      if (cf)
        fftw_free(cf);
      else
        error("cf memory free failed in fftconv_ola!");

      // Cleanup the plans.
      fftw_destroy_plan(p_forward);
      fftw_destroy_plan(p_backward);

      // Clean up FFTW
      //fftw_cleanup(); This seems to put Octave in an unstable state. Calling fftconv will crash Octave.

      return oct_retval;
    }

  } // if (data_format == DOUBLE_PRECISION)


  // ********************************************
  //
  // Single Precision
  //
  // ********************************************

  if (data_format == SINGLE_PRECISION) {

    // Get pointers to input data.
    const FloatMatrix tmp0 = args(0).float_matrix_value();
    const FloatMatrix tmp1 = args(1).float_matrix_value();
    sA = tmp0.data();
    sB = tmp1.data();

    if (nrhs == 4 || (nrhs == 5 && load_wisdom)) { // Normal mode.

      //
      //
      // Normal (non in-place) mode.
      //
      //

      in_place = false;

      running = true;

      // Allocate memory for the single precision
      // output matrix.
      FloatMatrix Ymat(A_M+B_M-1, A_N);
      sY = Ymat.fortran_vec();

      if (nthreads>1) { // Use threads

        // Allocate mem for the threads.
        threads = new std::thread[nthreads]; // Init thread data.
        if (!threads) {
          error("Failed to allocate memory for threads!");
          return oct_retval;
        }

        // Allocate local data.
        sD = (S_DATA*) malloc(nthreads*sizeof(S_DATA));
        if (!sD) {
          error("Failed to allocate memory for thread data!");
          return oct_retval;
        }

        for (thread_n = 0; thread_n < nthreads; thread_n++) {

          line_start = thread_n * A_N/nthreads;
          line_stop =  (thread_n+1) * A_N/nthreads;

          // Init local data.
          sD[thread_n].line_start = line_start; // Local start index;
          sD[thread_n].line_stop = line_stop; // Local stop index;
          sD[thread_n].A = sA;
          sD[thread_n].A_M = A_M;
          sD[thread_n].A_N = A_N;
          sD[thread_n].B = sB;
          sD[thread_n].B_M = B_M;
          sD[thread_n].B_N = B_N;
          sD[thread_n].block_len = block_len;
          sD[thread_n].Y = sY;

          // Start the threads.
          threads[thread_n] = std::thread(smp_process, &D[thread_n]);

        } // for (thread_n = 0; thread_n < nthreads; thread_n++)

        //
        // Wait for all threads to finish.
        //

        for (thread_n = 0; thread_n < nthreads; thread_n++)
          threads[thread_n].join();

        // Free memory.
        if (sD)
          free((void*) sD);

      } else { // Do not use threads.

        float *y = sc; // We use the same data space here.
        s_a_block = (float*) fftw_malloc(block_len*sizeof(float));

        if (B_N > 1) {// B is a matrix.

          for (n=0; n<A_N; n++) {

            //
            // The overlap-and-add algorithm.
            //

            k = 0;
            while (k < A_M ) {

              // Convolve one block
              if (k+block_len < A_M)

                sfftconv( &sA[k+n*A_M], block_len, &sB[0+n*B_M], B_M, y, sa,sb,sc,saf,sbf,scf);

              else {

                // We must do the convolution (eg. the FFT) of the same lenght here (=fft_len) but the first
                // arg is now shorter than above so copy and zero-pad the first arg to fftconv.

                // Copy.
                //for (l=0; l<(A_M-k); l++)
                // a_block[l] = A[k+l + n*A_M];
                memcpy(s_a_block,&sA[k + n*A_M],(A_M-k)*sizeof(float));

                // Zero-pad.
                //for (l=A_M-k; l<block_len; l++)
                // a_block[l] = 0.0; // Zero-pad.
                memset(&a_block[A_M-k],0,(block_len-(A_M-k))*sizeof(float));

                sfftconv( s_a_block, block_len, &sB[0+n*B_M], B_M, y, sa,sb,sc,saf,sbf,scf);
              }

              // Add the overlap.
              if (k+fft_len <  A_M+B_M-1) {
                for (l=0; l<fft_len; l++)
                  sY[k+l + n*(A_M+B_M-1)] += y[l];
              } else {
                for (l=0; l<(A_M+B_M-1)-k; l++)
                  sY[k+l + n*(A_M+B_M-1)] += y[l];
              }

              k += block_len;
            } // while

            if (running==false) {
              printf("fftconv_ola: bailing out!\n");
              break;
            }

          } // end-for

        } else { // B is a vector.

          for (n=0; n<A_N; n++) {

            //
            // The overlap-and-add algorithm.
            //

            k = 0;
            while (k < A_M ) {

              // Convolve one block
              if (k+block_len < A_M)

                sfftconv( &sA[k+n*A_M], block_len, sB, B_M, y, sa,sb,sc,saf,sbf,scf);

              else {

                // We must do the convolution (eg. the FFT) of the same lenght here (=fft_len) but the first
                // arg is now shorter than above so copy and zero-pad the first arg to fftconv.

                // Copy.
                //for (l=0; l<(k+block_len)-A_M; l++)
                // a_block[l] = A[k+l + n*A_M];
                memcpy(s_a_block,&sA[k + n*A_M],(A_M-k)*sizeof(float));

                // Zero-pad.
                //for (l=A_M-k; l<block_len; l++)
                // a_block[l] = 0.0; // Zero-pad.
                memset(&s_a_block[A_M-k],0,(block_len-(A_M-k))*sizeof(float));

                sfftconv( s_a_block, block_len, sB, B_M, y, sa,sb,sc,saf,sbf,scf);
              }

              // Add the overlap.
              if (k+fft_len <  A_M+B_M-1) {
                for (l=0; l<fft_len; l++)
                  sY[k+l + n*(A_M+B_M-1)] += y[l];
              } else {
                for (l=0; l<(A_M+B_M-1)-k; l++)
                  sY[k+l + n*(A_M+B_M-1)] += y[l];
              }

              k += block_len;
            } // while

            if (running==false) {
              printf("fftconv_ola: bailing out!\n");
              break;
            }

          } // end-for
        } // end-if

        // Free buffer memory.
        if (s_a_block)
          fftw_free(s_a_block);
        else
          error("s_a_block memory free failed in fftconv_ola!");

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
        the_str = fftwf_export_wisdom_to_string();
        buflen = strlen(the_str);

        std::string cmout( buflen, ' ' );
        cmout.insert( buflen, (const char*) the_str);

        // Add to output args.
        oct_retval.append( cmout);

        fftw_free(the_str);
      }

      // Clear temp vectors used for the FFTW plans.

      if (sa)
        fftw_free(sa);
      else
        error("sa memory free failed in fftconv_ola!");

      if (saf)
        fftw_free(saf);
      else
        error("saf memory free failed in fftconv_ola!");

      if (sb)
        fftw_free(sb);
      else
        error("sb memory free failed in fftconv_ola!");

      if (sbf)
        fftw_free(sbf);
      else
        error("sbf memory free failed in fftconv_ola!");

      if (sc)
        fftw_free(sc);
      else
        error("sc memory free failed in fftconv_ola!");

      if (scf)
        fftw_free(scf);
      else
        error("scf memory free failed in fftconv_ola!");

      // Cleanup the plans.
      fftwf_destroy_plan(sp_forward);
      fftwf_destroy_plan(sp_backward);

      // Clean up FFTW
      //fftw_cleanup(); This seems to put Octave in an unstable state. Calling fftconv will crash Octave.

      return oct_retval;

    }


    if ( (nrhs == 4 && !load_wisdom) || nrhs == 5 || nrhs == 6) { // In-place mode.

      in_place = true;

      //
      //
      // In-place mode.
      //
      //

      if (  args(3).matrix_value().rows() != A_M+B_M-1) {
        error("Wrong number of rows in argument 4!");
        return oct_retval;
      }

      if (  args(3).matrix_value().cols() != A_N) {
        error("Wrong number of columns in argument 4!");
        return oct_retval;
      }

      // Check datatype of (iun-place) output arg.
      if (!args(3).is_single_type()) {
        error("The output matrix (4th arg) must be single precision matrix when the input matrcies are single precision!");
        return oct_retval; // FIXME : Do we need to clear memory here?
      }

      const FloatMatrix Ytmp = args(3).matrix_value();
      sY = (float*) Ytmp.data();

      //
      // Call the CONV subroutine.
      //

      running = true;

      if (nthreads>1) { // Use threads

        // Allocate mem for the threads.
        threads = new std::thread[nthreads]; // Init thread data.
        if (!threads) {
          error("Failed to allocate memory for threads!");
          return oct_retval;
        }

        //
        // Double precision
        //


        // Allocate local data.
        sD = (S_DATA*) malloc(nthreads*sizeof(S_DATA));
        if (!sD) {
          error("Failed to allocate memory for thread data!");
          return oct_retval;
        }

        for (thread_n = 0; thread_n < nthreads; thread_n++) {

          line_start = thread_n * A_N/nthreads;
          line_stop =  (thread_n+1) * A_N/nthreads;

          // Init local data.
          sD[thread_n].line_start = line_start; // Local start index;
          sD[thread_n].line_stop = line_stop; // Local stop index;
          sD[thread_n].A = sA;
          sD[thread_n].A_M = A_M;
          sD[thread_n].A_N = A_N;
          sD[thread_n].B = sB;
          sD[thread_n].B_M = B_M;
          sD[thread_n].B_N = B_N;
          sD[thread_n].block_len = block_len;
          sD[thread_n].Y = sY;

          // Start the threads.
          threads[thread_n] = std::thread(smp_process, &D[thread_n]);

        }  // for (thread_n = 0; thread_n < nthreads; thread_n++)

        // Wait for all threads to finish.
        for (thread_n = 0; thread_n < nthreads; thread_n++)
          threads[thread_n].join();

        // Free memory.
        if (sD)
          free((void*) sD);

      } else { // Do not use threads.

        float *y = sc; // We use the same data space here.
        s_a_block = (float*) fftw_malloc(block_len*sizeof(float));

        if (B_N > 1) {// B is a matrix.

          for (n=0; n<A_N; n++) {

            //
            // The overlap-and-add algorithm.
            //

            k = 0;
            while (k < A_M ) {

              // Convolve one block
              if (k+block_len < A_M)

                sfftconv( &sA[k+n*A_M], block_len, &sB[0+n*B_M], B_M, y, sa,sb,sc,saf,sbf,scf);

              else {

                // We must do the convolution (eg. the FFT) of the same lenght here (=fft_len) but the first
                // arg is now shorter than above so copy and zero-pad the first arg to fftconv.

                // Copy.
                //for (l=0; l<(A_M-k); l++)
                //  a_block[l] = A[k+l + n*A_M];
                memcpy(s_a_block,&sA[k + n*A_M],(A_M-k)*sizeof(float));

                // Zero-pad.
                //for (l=A_M-k; l<block_len; l++)
                // a_block[l] = 0.0; // Zero-pad.
                memset(&s_a_block[A_M-k],0,(block_len-(A_M-k))*sizeof(float));

                sfftconv( s_a_block, block_len, &sB[0+n*B_M], B_M, y, sa,sb,sc,saf,sbf,scf);
              }

              // Add the overlap.
              if (k+fft_len <  A_M+B_M-1) {
                for (l=0; l<fft_len; l++)
                  sY[k+l + n*(A_M+B_M-1)] += y[l];
              } else {
                for (l=0; l<(A_M+B_M-1)-k; l++)
                  sY[k+l + n*(A_M+B_M-1)] += y[l];
              }

              k += block_len;
            } // while

            if (running==false) {
              printf("fftconv_ola: bailing out!\n");
              break;
            }

          } // end-for
        } else { // B is a vector.

          for (n=0; n<A_N; n++) {

            //
            // The overlap-and-add algorithm.
            //

            k = 0;
            while (k < A_M ) {

              // Convolve one block
              if (k+block_len < A_M)

                sfftconv( &sA[k+n*A_M], block_len, sB, B_M, y, sa,sb,sc,saf,sbf,scf);

              else {

                // We must do the convolution (eg. the FFT) of the same lenght here (=fft_len) but the first
                // arg is now shorter than above so copy and zero-pad the first arg to fftconv.

                // Copy.
                //for (l=0; l<(A_M-k); l++)
                // a_block[l] = A[k+l + n*A_M];
                memcpy(s_a_block,&sA[k + n*A_M],(A_M-k)*sizeof(float));

                // Zero-pad.
                //for (l=0; l<(A_M-k); l++)
                // a_block[l] = 0.0; // Zero-pad.
                memset(&s_a_block[A_M-k],0,(block_len-(A_M-k))*sizeof(float));

                sfftconv( s_a_block, block_len, sB, B_M, y, sa,sb,sc,saf,sbf,scf);
              }

              // Add the overlap.
              if (k+fft_len <  A_M+B_M-1) {
                for (l=0; l<fft_len; l++)
                  sY[k+l + n*(A_M+B_M-1)] += y[l];
              } else {
                for (l=0; l<(A_M+B_M-1)-k; l++)
                  sY[k+l + n*(A_M+B_M-1)] += y[l];
              }

              k += block_len;
            } // while

            if (running==false) {
              printf("fftconv_ola: bailing out!\n");
              break;
            }

          } // end-for
        } // end-if


        // Free buffer memory.
        if (s_a_block)
          fftw_free(s_a_block);
        else
          error("s_a_block memory free failed in fftconv_ola!");

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

      if (sa)
        fftw_free(sa);
      else
        error("sa memory free failed in fftconv_ola!");

      if (saf)
        fftw_free(saf);
      else
        error("saf memory free failed in fftconv_ola!");

      if (sb)
        fftw_free(sb);
      else
        error("sb memory free failed in fftconv_ola!");

      if (sbf)
        fftw_free(sbf);
      else
        error("sbf memory free failed in fftconv_ola!");

      if (sc)
        fftw_free(sc);
      else
        error("sc memory free failed in fftconv_ola!");

      if (scf)
        fftw_free(scf);
      else
        error("scf memory free failed in fftconv_ola!");

      // Cleanup the (single precision) plans.
      fftwf_destroy_plan(sp_forward);
      fftwf_destroy_plan(sp_backward);

      // Clean up FFTW
      fftw_cleanup();

      return oct_retval;
    }

  } // if (data_format == SINGLE_PRECISION)

  return oct_retval;
}
