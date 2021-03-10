/***
*
* Copyright (C) 2010,2011,2012,2014,2015,2016,2021 Fredrik Lingvall
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

#include "dream.h"
#include "fftconv.h"

#include <fftw3.h>

//
// Octave headers.
//

#include <octave/oct.h>

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
int plan_method = 4; // Default to ESTIMATE method.

// FFTW plans.
fftw_plan    p_forward;
fftw_plan    p_backward;

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
  dream_idx_type block_len;
  double *Y;
  FFT *fft;
  ConvMode conv_mode;
} DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//

void* smp_dream_fftconv_ola(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

/***
 *
 * Thread function
 *
 ***/

void* smp_dream_fftconv_ola(void *arg)
{
  DATA D = *(DATA *)arg;
  dream_idx_type col_start=D.col_start, col_stop=D.col_stop;
  double *A = D.A, *B = D.B;
  double *Y = D.Y;
  FFT fft = *D.fft;
  dream_idx_type A_M = D.A_M, B_M = D.B_M, B_N = D.B_N;
  dream_idx_type block_len = D.block_len;
  ConvMode conv_mode = D.conv_mode;
  dream_idx_type fft_len = block_len + B_M - 1;

  // Input vectors.
  FFTVec a_v(fft_len);
  FFTVec b_v(fft_len);
  FFTVec c_v(fft_len);
  double *a = a_v.get(), *b = b_v.get(), *c = c_v.get();

  FFTVec y_v(fft_len);
  double *y = y_v.get();

  // Fourier Coefficients.
  FFTCVec af_v(fft_len);
  FFTCVec bf_v(fft_len);
  FFTCVec cf_v(fft_len);
  std::complex<double> *af = af_v.get(), *bf = bf_v.get(), *cf = cf_v.get();

  FFTVec a_block_v(fft_len);
  double *a_block = a_block_v.get();

  //
  // Do the convolution.
  //

  if (B_N > 1) {// B is a matrix.

    dream_idx_type n=col_start;
    double *Ap = &A[0+n*A_M];
    double *Bp = &B[0+n*B_M];

    for (n=col_start; n<col_stop; n++) {

      //
      // The overlap-and-add algorithm.
      //

      dream_idx_type k = 0;
      dream_idx_type len = 0;
      while (k < A_M ) {

        //
        // Convolve one block
        //

        if (k+block_len < A_M) {

          fftconv(fft, Ap, block_len, Bp, B_M, y,
                  a, b, c, af, bf, cf, conv_mode);

          len = block_len;

        } else {

          // We must do the convolution (eg. the FFT) of the same lenght here (=fft_len) but the first
          // arg is now shorter than above so copy and zero-pad the first arg to fftconv.

          len = A_M-k;

          // Copy.
          memcpy(a_block, &A[k + n*A_M], len*sizeof(double));

          // Zero-pad.
          memset(&a_block[A_M-k], 0, (block_len - len)*sizeof(double));

          fftconv(fft, a_block, block_len, Bp, B_M, y,
                  a, b, c, af, bf, cf, conv_mode);
        }

        Ap += len;

        //
        // Add the overlap.
        //

        if (k+fft_len <  A_M+B_M-1) {
          for (dream_idx_type l=0; l<fft_len; l++) {
            Y[k+l + n*(A_M+B_M-1)] += y[l];
          }
        } else {
          for (dream_idx_type l=0; l<(A_M+B_M-1)-k; l++) {
            Y[k+l + n*(A_M+B_M-1)] += y[l];
          }
        }

        k += block_len;
      }

      if (running==false) {
        octave_stdout << "fftconv_ola: thread for column " << col_start+1 << " -> " <<
          col_stop << " bailing out at column = " << n <<"!\n";
        break;
      }

      Bp += B_M;

    } // end-for

  } else { // B is a vector.

    dream_idx_type n=col_start;
    double *Ap = &A[0+n*A_M];

    for (n=col_start; n<col_stop; n++) {

      //
      // The overlap-and-add algorithm.
      //

      dream_idx_type k = 0;
      dream_idx_type len = 0;
      while (k < A_M ) {

        // Convolve one block
        if (k+block_len < A_M) {

          fftconv(fft, Ap, block_len, B, B_M, y,
                  a, b, c, af, bf, cf, conv_mode);

          len = block_len;

        } else {

          // We must do the convolution (eg. the FFT) of the same lenght here (=fft_len) but the first
          // arg is now shorter than above so copy and zero-pad the first arg to fftconv.

          len = A_M-k;

          // Copy.
          memcpy(a_block, &A[k + n*A_M], len*sizeof(double));

          // Zero-pad.
          memset(&a_block[A_M-k], 0, (block_len - len)*sizeof(double));

          fftconv(fft, a_block, block_len, B, B_M, y,
                  a, b, c, af, bf, cf, conv_mode);

        }

        Ap += len;

        // Add the overlap.
        if (k+fft_len < A_M+B_M-1) {
          for (dream_idx_type l=0; l<fft_len; l++) {
            Y[k+l + n*(A_M+B_M-1)] += y[l];
          }
        } else {
          for (dream_idx_type l=0; l<(A_M+B_M-1)-k; l++) {
            Y[k+l + n*(A_M+B_M-1)] += y[l];
          }
        }

        k += block_len;
      } // while

      if (running==false) {
        octave_stdout << "fftconv_ola: thread for column " << col_start+1 << " -> " <<
          col_stop << " bailing out at column = " << n <<"!\n";
        break;
      }

      Ap += A_M;

    } // for
  } // if

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
Copyright @copyright{} 2010-2021 Fredrik Lingvall.\n\
@seealso {conv_p, fftconv_p, fftconv, conv, fftw_wisdom}\n\
@end deftypefn")
{
  sighandler_t   old_handler, old_handler_abrt, old_handler_keyint;
  std::thread     *threads;
  dream_idx_type fft_len, return_wisdom = false, load_wisdom = false;
  char *the_str = nullptr;
  int buflen, is_set = false;
  dream_idx_type block_len;
  ConvMode conv_mode=ConvMode::equ;

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
  case 1:			// A matrix/vector
  case 2:			// B matrix/vector
    {
      error("fftconv_ola requires 3 to 6 input arguments!");
      return oct_retval;
    }
    break;

  case 3:			// OLA bleock length
    {
      if (nlhs > 2) {
        error("Too many output arguments for fftconv_ola!");
        return oct_retval;
      }

      if (nlhs == 2) {
        return_wisdom = true;
      }
    }
    break;

  case 4: // Output arg in in-place mode or a wisdom string in normal mode.
    {
      if ( args(3).is_string() ) { // 4th arg is a fftw wisdom string.

        std::string strin = args(3).string_value();
        buflen = strin.length();
        the_str = (char*) fftw_malloc(buflen * sizeof(char));
        for (dream_idx_type n=0; n<buflen; n++ ) {
          the_str[n] = strin[n];
        }

        //
        // If 4:th arg is a string then only a wisdom string is valid.
        //

        if (!strcmp("fftw_wisdom",the_str) && (!strcmp("fftwf_wisdom",the_str))) {
          error("The string in arg 4 do not seem to be in a fftw wisdom format!");
          return oct_retval;
        }
        else {
          load_wisdom = true;
        }

      } else { // 4th arg not a string then assume in-place mode.
        fftw_forget_wisdom(); // Clear wisdom history (a new wisdom will be created below).
        if (nlhs > 0) {
          error("4th arg is not a fftw wisdom string and in-place mode is assumed. But then there should be no output args!");
          return oct_retval;
        }
      }
    }
    break;

  case 5:  // In-place mode if >= 5 args.
    {
      if ( args(4).is_string() ) { // 5th arg is a string ('=','+=','-=', or fftw wisdom).

        std::string strin = args(4).string_value();
        buflen = strin.length();
        the_str = (char*) fftw_malloc(buflen * sizeof(char));
        for (dream_idx_type n=0; n<buflen; n++ ) {
          the_str[n] = strin[n];
        }

        // Valid strings are:
        //  '='  : In-place replace mode.
        //  '+=' : In-place add mode.
        //  '-=' : In-place sub mode.
        //  wisdom string.

        is_set = false;

        if (strncmp(the_str,"=",1) == 0) {
          conv_mode = ConvMode::equ;
          is_set = true;
        }

        if (strncmp(the_str,"+=",2) == 0) {
          conv_mode = ConvMode::sum;
          is_set = true;
        }

        if (strncmp(the_str,"-=",2) == 0) {
          conv_mode = ConvMode::neg;
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
    }
    break;

  case 6: // In-place mode if 6 input args.
    {
      if ( args(5).is_string() ) { // 6:th arg is a string (fftw wisdom).

        // Read the wisdom string.
        std::string strin = args(5).string_value();
        buflen = strin.length();
        the_str = (char*) fftw_malloc(buflen * sizeof(char));
        for (dream_idx_type n=0; n<buflen; n++ ) {
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
    }
    break;

  default:
    {
      error("fftconv_ola requires 3 to 6 input arguments!");
      return oct_retval;
    }
    break;
  }

  const Matrix tmp0 = args(0).matrix_value();
  dream_idx_type A_M = tmp0.rows();
  dream_idx_type A_N = tmp0.cols();
  double *A = (double*) tmp0.fortran_vec();

  const Matrix tmp1 = args(1).matrix_value();
  dream_idx_type B_M = tmp1.rows();
  dream_idx_type B_N = tmp1.cols();
  double *B = (double*) tmp1.fortran_vec();

  // Check that arg 2.
  if ( B_M != 1 && B_N !=1 && B_N != A_N) {
    error("Argument 2 must be a vector or a matrix with the same number of rows as arg 1!");
    return oct_retval;
  }

  if ( B_M == 1 || B_N == 1 ) { // B is a vector.
    B_M = B_M*B_N;
    B_N = 1;
  }

  //
  // The segment length (block length).
  //

  if (mxGetM(2)*mxGetN(2) == 1) {
    block_len = (dream_idx_type) args(2).matrix_value().data()[0];
  } else {
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
  dream_idx_type nthreads = std::thread::hardware_concurrency();

  // Read DREAM_NUM_THREADS env var
  if(const char* env_p = std::getenv("DREAM_NUM_THREADS")) {
    dream_idx_type dream_threads = std::stoul(env_p);
    if (dream_threads < nthreads) {
      nthreads = dream_threads;
    }
  }

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

  if ((old_handler = std::signal(SIGTERM, &sighandler)) == SIG_ERR) {
    std::cerr << "Couldn't register SIGTERM signal handler!" << std::endl;
  }

  if ((old_handler_abrt = std::signal(SIGABRT, &sighandler)) == SIG_ERR) {
    std::cerr << "Couldn't register SIGABRT signal handler!" << std::endl;
  }

  if ((old_handler_keyint = std::signal(SIGINT, &sighandler)) == SIG_ERR) {
    std::cerr << "Couldn't register SIGINT signal handler!" << std::endl;
  }

  // ********************************************
  //
  // Init the FFTW plans.
  //
  // ********************************************

  if(load_wisdom) {
    if (!fftw_import_wisdom_from_string(the_str)) {
      error("Failed to load (double precision) fftw wisdom!");
      return oct_retval;
    } else {
      fftw_free(the_str); // Clean up.
    }
  }

  //std::mutex fft_mutex;
  //FFT fft(fft_len, &fft_mutex, plan_method);
  //
  // NB We call the FFTW planners only from the main thread once
  // so no need to lock it with a mutex

  fft_len = block_len + B_M - 1;
  FFT fft(fft_len, nullptr, plan_method);

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
    double *Y = Ymat.fortran_vec();

    SIRData ymat(Y, A_M+B_M-1, A_N);
    ymat.clear();

    // Init thread data.
    threads = new std::thread[nthreads];
    if (!threads) {
      error("Failed to allocate memory for threads!");
      return oct_retval;
    }

    // Allocate thread data.
    DATA *D = (DATA*) malloc(nthreads*sizeof(DATA));
    if (!D) {
      error("Failed to allocate memory for thread data!");
      return oct_retval;
    }

    for (dream_idx_type thread_n=0; thread_n<nthreads; thread_n++) {

      dream_idx_type col_start = thread_n * A_N/nthreads;
      dream_idx_type col_stop = (thread_n+1) * A_N/nthreads;

      // Init local data.
      D[thread_n].col_start = col_start; // Local start index;
      D[thread_n].col_stop = col_stop; // Local stop index;
      D[thread_n].A = A;
      D[thread_n].A_M = A_M;
      D[thread_n].A_N = A_N;
      D[thread_n].B = B;
      D[thread_n].B_M = B_M;
      D[thread_n].B_N = B_N;
      D[thread_n].block_len = block_len;
      D[thread_n].Y = Y;
      D[thread_n].fft = &fft;
      D[thread_n].conv_mode = conv_mode;

      if (nthreads > 1) {
        // Start the threads.
        threads[thread_n] = std::thread(smp_dream_fftconv_ola, &D[thread_n]);
      } else {
        smp_dream_fftconv_ola(&D[0]);
      }

    } // for (thread_n = 0; thread_n < nthreads; thread_n++)

    //
    // Wait for all threads to finish.
    //

    if (nthreads > 1) {
      for (dream_idx_type thread_n=0; thread_n < nthreads; thread_n++) {
        threads[thread_n].join();
      }
    }

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
      the_str = fftw_export_wisdom_to_string();
      buflen = strlen(the_str);

      std::string cmout( buflen, ' ' );
      cmout.insert( buflen, (const char*) the_str);

      // Add to output args.
      oct_retval.append( cmout);

      fftw_free(the_str);
    }

    return oct_retval;
  }


  if ( (nrhs == 4 && !load_wisdom) || nrhs == 5 || nrhs == 6) { // In-place mode.

    in_place = true;

    if (args(3).matrix_value().rows() != A_M+B_M-1) {
      error("Wrong number of rows in argument 5!");
      return oct_retval;
    }

    if (args(3).matrix_value().cols() != A_N) {
      error("Wrong number of columns in argument 5!");
      return oct_retval;
    }

    const Matrix Ytmp = args(3).matrix_value();
    double *Y = (double*) Ytmp.data();

    SIRData ymat(Y, Ytmp.rows(), Ytmp.cols());
    ymat.clear();

    //
    // Call the CONV subroutine.
    //

    running = true;

    // Allocate mem for the threads.
    threads = new std::thread[nthreads]; // Init thread data.
    if (!threads) {
      error("Failed to allocate memory for threads!");
      return oct_retval;
    }

    // Allocate local data.
    DATA *D = (DATA*) malloc(nthreads*sizeof(DATA));
    if (!D) {
      error("Failed to allocate memory for thread data!");
      return oct_retval;
      }

    for (dream_idx_type thread_n=0; thread_n<nthreads; thread_n++) {

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
      D[thread_n].block_len = block_len;
      D[thread_n].Y = Y;

      // Start the threads.
      threads[thread_n] = std::thread(smp_dream_fftconv_ola, &D[thread_n]);

    } // for (thread_n = 0; thread_n < nthreads; thread_n++)

      // Wait for all threads to finish.
    for (dream_idx_type thread_n = 0; thread_n < nthreads; thread_n++)
      threads[thread_n].join();

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
      the_str = fftw_export_wisdom_to_string();
      buflen = strlen(the_str);

      std::string cmout( buflen, ' ' );
      cmout.insert( buflen, (const char*) the_str);

      // Add to output args.
      oct_retval.append( cmout);

      fftw_free(the_str);
    }

    // Cleanup the plans.
    fftw_destroy_plan(p_forward);
    fftw_destroy_plan(p_backward);

    // Clean up FFTW
    //fftw_cleanup(); This seems to put Octave in an unstable state. Calling fftconv will crash Octave.

    return oct_retval;
  }

  return oct_retval;
}
