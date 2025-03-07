/***
*
* Copyright (C) 2006,2007,2008,2009,2014,2015,2016,2021,2022,2023,2025 Fredrik Lingvall
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

#include <octave/oct.h>

#include "dream.h"
#include "affinity.h"
#include "fftconv.h"

/***
 *
 *  Parallel (threaded) FFTW based convolution and summation.
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
  int L;
  double **H;
  dream_idx_type H_M;
  dream_idx_type H_N;
  double *U;
  dream_idx_type U_M;
  dream_idx_type U_N;
  double *Y;
  FFT *fft;
  ConvMode conv_mode;
} DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//

void* smp_dream_sum_fftconv(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

/***
 *
 * Thread function.
 *
 ***/

void* smp_dream_sum_fftconv(void *arg)
{
  DATA D = *(DATA *)arg;
  dream_idx_type col_start=D.col_start, col_stop=D.col_stop, n;
  double **H = D.H, *U = D.U, *Y = D.Y;
  dream_idx_type H_M = D.H_M, U_M = D.U_M; //, U_N = D.U_N;
  dream_idx_type L = D.L;
  FFT fft = *D.fft;
  ConvMode conv_mode = D.conv_mode;

  dream_idx_type fft_len = H_M + U_M - 1;

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

  // NB. U is always a matrix here

  for (n=col_start; n<col_stop; n++) {

    add_fftconv(fft,
                H, L, H_M,
                n,
                U, U_M, // U must be U_M x L
                &Y[0+n*fft_len],
                a, b, c, af, bf, cf,
                conv_mode);

    if (running==false) {
      octave_stdout << "sum_fftconv: thread for column " << col_start+1 << " -> " << col_stop << " bailing out!\n";
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
 * Octave (oct) gateway function for SUM_FFTCONV.
 *
 ***/

DEFUN_DLD (sum_fftconv, args, nlhs,
           "-*- texinfo -*-\n\
@unnumberedsec Normal mode\n\
\n\
@deftypefn {Loadable Function} {} Y = sum_fftconv(H, U, wisdom_str);\n\
\n\
SUM_FFTCONV - Computes (using parallel threaded processing) the sum of one dimensional\n\
convolutions of the columns in each 2D matrix in the 3D matrix H\n\
with the corresponding columns in the matrix U.\n\
\n\
In normal mode sum_fftconv performs an operation similar to:\n\
@verbatim\n\
\n\
fft_len = M+K-1;\n\
YF = zeros(fft_len,N);\n\
for l=1:L\n\
  for n=1:N\n\
    YF(:,n) = YF(:,n) + fft(H(:,n,l),fft_len).* fft(U(:,l),fft_len);\n\
  end\n\
end\n\
Y = real(ifft(YF))/fft_len;\n\
\n\
@end verbatim\n\
using threaded processing. The computations are performed using FFT:s.\n\
\n\
Input parameters:\n\
\n\
@table @code\n\
@item H\n\
  An MxNxL 3D matrix.\n\
@item U\n\
  H KxL matrix.\n\
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
[tmp,wisdom_str] = fftconv_p(H(:,1,1),U(:,1));\n\
\n\
for i=1:N\n\
\n\
  % Do some stuff here.\n\
\n\
  Y = sum_fftconv(H, U, wisdom_str);\n\
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
@end deftypefn\n\
\n\
@unnumberedsec In-place mode\n\
\n\
In in-place mode sum_fftconv performs the operations in-place on a pre-allocated matrix:\n\
\n\
@verbatim\n\
 sum_fftconv(H, U, Y, wisdom_str);\n\
\n\
@end verbatim\n\
@noindent Here sum_fftconv do not have any output arguments and the\n \
results are instead stored directly in the pre-allocated (M+K-1)xN input matrix Y. A typical usage is:\n\
@verbatim\n\
\n\
Y = zeros(M+K-1,N); % Allocate space for Y.\n\
\n\
for i=1:N\n\
\n\
  % Do some stuff here.\n\
\n\
   sum_fftconv(H, U, Y, wisdom_str);\n\
end\n\
\n\
@end verbatim\n\
@noindent where memory allocation of the (possible large) matrix Y now is avoided inside the for loop.\n\
\n\
NOTE: A side-effect of using in-place mode is that if a copy Y2 of Y is made\n\
then both Y2 and Y will be altered by sum_fftconv. That is, by performing,\n\
@verbatim\n\
\n\
Y  = zeros(M+K-1,N);\n\
Y2 = Y; % No actual copy of data here.\n\
sum_fftconv(H,U,Y,wisdom_str);\n\
\n\
@end verbatim\n\
@noindent then both Y and Y2 will be changed (since Octave do not make a new copy of the data\n\
in Y2 unless Y2 is changed before the sum_fftconv call).\n\
\n\
sum_fftconv is a part of the DREAM Toolbox available at\n\
@url{https://github.com/frli8848/DREAM}.\n\
\n\
Copyright @copyright{} 2006-2023 Fredrik Lingvall.\n\
@seealso {conv, conv_p, fftconv, fftconv_p, fftw_wisdom}")
{
  sighandler_t old_handler, old_handler_abrt, old_handler_keyint;
  dream_idx_type col_start, col_stop, H_M, H_N, H_L, U_M, U_N;
  double *Y = nullptr;
  std::thread *threads;
  dream_idx_type thread_n, nthreads;
  DATA *D = nullptr;
  int plan_method = 4; // Default to FFTW_ESTIMATE
  dream_idx_type fft_len;
  bool return_wisdom = false, load_wisdom = false;
  ConvMode conv_mode=ConvMode::equ;

  octave_value_list oct_retval;

  int nrhs = args.length();

  //
  // Set the method which fftw computes plans
  //

  // If we want to save a plan (in the second output arg)
  // then use the more time-consuming MEAUSURE method.
  if (nlhs == 2) {
    plan_method = 3; // 3 = MEASURE.
  }

  //
  // Check for proper inputs arguments.
  //

  // Num inputs
  if ( (nrhs < 2) ||  (nrhs > 4) ) {
    dream_err_msg("sum_fftconv requires 2 to 4 input arguments!");
    return oct_retval;
  }

  // Num outputs
  if (nlhs > 2) {
    dream_err_msg("Too many output arguments for sum_fftconv!");
    return oct_retval;
  }

  // 2nd output is wisdom string
  if (nlhs == 2) {
    return_wisdom = true;
  }

  const NDArray tmp0 = args(0).array_value();
  dream_idx_type dims = args(0).ndims();
  if (dims != 3) {
    dream_err_msg("Argument 1 should be a 3D Matrix\n");
    return oct_retval;
  }

  dim_vector dv = args(0).dims();
  H_M = dv(0);
  H_N = dv(1);
  H_L = dv(2);

  // Store pointers to the L A-matrices in a vector.
  double **H = (double**) malloc(H_L*sizeof(double*));
  for (dream_idx_type k=0; k<H_L; k++) {
    H[k] = (double*) &(tmp0.data()[H_M*H_N*k]);
  }

  const Matrix tmp1 = args(1).matrix_value();
  U_M = tmp1.rows();
  U_N = tmp1.cols();
  double *U = (double*) tmp1.data();

  if (H_L != U_N) {
    dream_err_msg("3rd dimension of arg 1 must match the number of columns in arg 2\n");
    return oct_retval;
  }

  if (U_M == 1 || U_N == 1 ) { // U is a vector.
    U_M = U_M*U_N;
    U_N = 1;
  }

  fft_len = H_M+U_M-1;

  // NB We call the FFTW planners only from the main thread once
  // so no need to lock it with a mutex
  FFT fft(fft_len, nullptr, plan_method);

  std::string wisdom_str ="";

  switch (nrhs) {

  case 2:
    break;

  case 3:
    if ( args(2).is_string() ) { // 3rd arg is a fftw wisdom string.

      wisdom_str = args(2).string_value();

      //
      // If 3rd arg is a string then only a wisdom string is valid.
      //

      if (!fft.is_wisdom(wisdom_str)) {
        dream_err_msg("The string in arg 3 do not seem to be in a FFTW wisdom format!");
        return oct_retval;
      }
      else {
        load_wisdom = true;
      }
    } else { // 3rd arg not a string then assume in-place mode.
      fft.forget_wisdom(); // Clear wisdom history (a new wisdom will be created below).
      if (nlhs > 0) {
        dream_err_msg("3rd arg is not a FFTW wisdom string and in-place mode is assumed. But then there should be no output args!");
        return oct_retval;
      }
    }
    break;

  case 4:
    if ( args(3).is_string() ) { // 4th arg is a fftw wisdom string.

      wisdom_str = args(3).string_value();

      if (!fft.is_wisdom(wisdom_str)) {
        dream_err_msg("The string in arg 4 do not seem to be in fftw wisdom format!");
        return oct_retval;
      }
      else {
        load_wisdom = true;
      }
    }
    else { // 4th arg not a string
      dream_err_msg("Argument 4 do not seem to be in fftw wisdom string format!");
      return oct_retval;
    }
    break;

  default:
    dream_err_msg("sum_fftconv requires 2 to 4 input arguments!");
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

  // We cannot have more threads than the number of observation points.
  if (nthreads > H_N) {
    nthreads = H_N;
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
      dream_err_msg("Failed to load FFTW wisdom!");
      return oct_retval;
    }
  }

  //
  // Normal (non in-place) mode.
  //

  if (nrhs == 2 ||  (nrhs == 3 && load_wisdom)) {

    Matrix Ymat(fft_len, H_N);
    Y = (double*) Ymat.data();

    // Clear output in normal mode
    SIRData ymat(Y, fft_len, H_N);
    ymat.clear();

    oct_retval.append(Ymat);
  }

  //
  // In-place mode.
  //

  if ( (nrhs == 3 && !load_wisdom) || nrhs == 4 ) {

    conv_mode = ConvMode::sum;

    if (args(2).matrix_value().rows() != fft_len) {
      dream_err_msg("Wrong number of rows in argument 3!");
      return oct_retval;
    }

    if (args(2).matrix_value().cols() != H_N) {
      dream_err_msg("Wrong number of columns in argument 3!");
      return oct_retval;
    }

    const Matrix Ytmp = args(2).matrix_value(); // This will make it stay in scope.
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
  threads = new std::thread[nthreads]; // Init thread data.
  if (!threads) {
    dream_err_msg("Failed to allocate memory for threads!");
    return oct_retval;
  }

  for (thread_n=0; thread_n < nthreads; thread_n++) {

    col_start = thread_n * H_N/nthreads;
    col_stop = (thread_n+1) * H_N/nthreads;

    // Init local data.
    D[thread_n].col_start = col_start; // Local start index;
    D[thread_n].col_stop = col_stop; // Local stop index;
    D[thread_n].H = H;
    D[thread_n].H_M = H_M;
    D[thread_n].H_N = H_N;
    D[thread_n].U = U;
    D[thread_n].U_M = U_M;
    D[thread_n].U_N = U_N;
    D[thread_n].Y = Y;
    D[thread_n].L = H_L;
    D[thread_n].fft = &fft;
    D[thread_n].conv_mode = conv_mode;

    if (nthreads > 1) {
      // Start the threads.
      threads[thread_n] = std::thread(smp_dream_sum_fftconv, &D[thread_n]);
    } else {
      smp_dream_sum_fftconv(&D[0]);
    }
  }

  if (nthreads > 1) {
    // Wait for all threads to finish.
    for (dream_idx_type thread_n = 0; thread_n < nthreads; thread_n++) {
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

  if (H) {
    free(H);
  }

  // Return the FFTW Wisdom so that the plans can be re-used.
  if (return_wisdom) {
    std::string cmout = fft.get_wisdom();
    oct_retval.append(cmout); // Add to output args.
  }

  return oct_retval; // Just to fix compiler warnings.
}
