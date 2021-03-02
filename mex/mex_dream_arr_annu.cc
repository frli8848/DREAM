/***
*
* Copyright (C) 2003,2005,2006,2007,2008,2009,2014,2015,2019 Fredrik Lingvall
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
#include <signal.h>

#include <iostream>
#include <thread>
#include <mutex>

#include "dream_arr_annu.h"
#include "affinity.h"
#include "dream_error.h"

#include "mex.h"

//
// Globals
//

volatile int out_err = NONE;
std::mutex err_lock;
int running;

//
// typedef:s
//

typedef struct
{
  size_t no;
  size_t start;
  size_t stop;
  double *ro;
  double dx;
  double dy;
  double dt;
  size_t nt;
  int delay_method;
  double *delay;
  double v;
  double cp;
  Attenuation *att;
  size_t num_radii;
  double *gr;
  int foc_met;
  bool do_apod;
  int apod_met;
  double *focal;
  double *apod;
  double param;
  double *h;
  int err_level;
} DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//

void* smp_dream_arr_annu(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

/***
 *
 * Thread function.
 *
 ***/

void* smp_dream_arr_annu(void *arg)
{
  int tmp_err = NONE, err = NONE;
  size_t n;
  DATA D = *(DATA *)arg;
  double xo, yo, zo;
  double *h = D.h;
  double dx=D.dx, dy=D.dy, dt=D.dt;
  size_t no=D.no, nt=D.nt;
  int    tmp_lev, err_level=D.err_level;
  double *delay=D.delay, *ro=D.ro, v=D.v, cp=D.cp;
  Attenuation *att = D.att;
  size_t start=D.start, stop=D.stop;
  int    foc_met=D.foc_met, do_apod = D.do_apod,apod_met=D.apod_met;
  double *focal=D.focal, *apod=D.apod, param=D.param;
  size_t  num_radii = D.num_radii;

  double *gr=D.gr;

  // Buffers for the FFTs in the Attenuation
  std::unique_ptr<FFTCVec> xc_vec;
  std::unique_ptr<FFTVec> x_vec;
  if (att) {
    xc_vec = std::make_unique<FFTCVec>(nt);
    x_vec = std::make_unique<FFTVec>(nt);
  }

  // Let the thread finish and then catch the error.
  if (err_level == STOP) {
    tmp_lev = PARALLEL_STOP;
  } else {
    tmp_lev = err_level;
  }

  for (n=start; n<stop; n++) {
    xo = ro[n];
    yo = ro[n+1*no];
    zo = ro[n+2*no];

    double dlay = 0.0;
    if (D.delay_method == SINGLE_DELAY) {
      dlay = delay[0];
    } else { // MULTIPLE_DELAYS delays.
      dlay = delay[n];
    }

    if (att == nullptr) {
      err = dream_arr_annu(xo, yo, zo,
                           dx, dy, dt, nt,
                           dlay, v,cp,
                           num_radii,gr,
                           foc_met, focal,
                           apod, do_apod, apod_met, param,
                           &h[n*nt], tmp_lev);
    } else {
      err = dream_arr_annu(*att, *xc_vec, *x_vec,
                           xo, yo, zo,
                           dx, dy, dt, nt,
                           dlay, v,cp,
                           num_radii,gr,
                           foc_met, focal,
                           apod, do_apod, apod_met, param,
                           &h[n*nt], tmp_lev);
    }

    if (err != NONE || out_err ==  PARALLEL_STOP) {
      tmp_err = err;
      if (err == PARALLEL_STOP || out_err ==  PARALLEL_STOP) {
        break; // Jump out when a STOP error occurs.
      }
    }

    if (!running) {
      std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!" << std::endl;
      return(NULL);
    }

  }

  // Lock out_err for update, update it, and unlock.
  err_lock.lock();

  if ((tmp_err != NONE) && (out_err == NONE)) {
    out_err = tmp_err;
  }

  err_lock.unlock();

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
 * Matlab (MEX) gateway function for dream_arr_annu
 *
 ***/

extern void _main();

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *ro, *s_par,*m_par;
  char   apod_str[50], foc_str[50];
  int    buflen;
  double dx,dy,dt;
  size_t nt, no;
  double param=0,*delay,v,cp,alpha;
  size_t num_radii=0;
  double *gr;
  int    foc_met=0;
  double *focal=nullptr;
  double *apod=NULL;
  bool   do_apod=false;
  int    apod_met=0;
  double *h, *err_p;
  int    err_level=STOP, is_set = false;
  char   err_str[50];
  DATA   *D;
  size_t start, stop;
  std::thread *threads;
  unsigned int thread_n, nthreads;
  sighandler_t old_handler, old_handler_abrt, old_handler_keyint;

  // Check for proper number of arguments.
  if (!((nrhs == 10) || (nrhs == 11))) {
    dream_err_msg("dream_arr_annu requires 10 or 11 input arguments!");
  } else {
    if (nlhs > 2) {
      dream_err_msg("Too many output arguments for dream_arr_annu!");
    }
  }

  //
  // Observation point.
  //

  // Check that arg (number of observation points) x 3 matrix
  if (mxGetN(prhs[0]) != 3) {
    dream_err_msg("Argument 1 to dream_arr_annu must be a (number of observation points) x 3 matrix!");
  }

  no = mxGetM(prhs[0]); // Number of observation points.
  ro = mxGetPr(prhs[0]);

  //
  // Grid function (position vectors of the elements).
  //

  if ((mxGetM(prhs[1]) > 1) & (mxGetN(prhs[1]) > 1)) {
    dream_err_msg("Argument 2 must a vector (number of array elements)");
  }

  num_radii = (int) mxGetM(prhs[1])*mxGetN(prhs[1]); // Number of elementents in the array.
  gr = mxGetPr(prhs[1]);	// Vector of annular radi,

  //
  // Temporal and spatial sampling parameters.
  //


  // Check that arg 2 is a 4 element vector
  if (!((mxGetM(prhs[2])==4 && mxGetN(prhs[2])==1) || (mxGetM(prhs[2])==1 && mxGetN(prhs[2])==4))) {
    dream_err_msg("Argument 3 must be a vector of length 4!");
  }

  s_par = mxGetPr(prhs[2]);
  dx    = s_par[0];		// Spatial x discretization size.
  dy    = s_par[1];		// Spatial dy iscretization size.
  dt    = s_par[2];		// Temporal discretization size (= 1/sampling freq).
  nt    = (size_t) s_par[3];	// Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //

  // Check that arg 4 is a scalar.
  if ( (mxGetM(prhs[3]) * mxGetN(prhs[3]) !=1) && ((mxGetM(prhs[3]) * mxGetN(prhs[3])) != no)) {
    dream_err_msg("Argument 4 must be a scalar or a vector with a length equal to the number of observation points!");
  }

  delay = mxGetPr(prhs[3]);

  //
  // Material parameters
  //

  // Check that arg 5 is a 3 element vector.
  if (!((mxGetM(prhs[4])==3 && mxGetN(prhs[4])==1) || (mxGetM(prhs[4])==1 && mxGetN(prhs[4])==3))) {
    dream_err_msg("Argument 5 must be a vector of length 3!");
  }

  m_par = mxGetPr(prhs[4]);
  v     = m_par[0]; // Normal velocity of transducer surface.
  cp    = m_par[1]; // Sound speed.
  alpha = m_par[2]; // Attenuation coefficient [dB/(cm MHz)],

  //
  // Focusing parameters.
  //

  if (nrhs >= 6) {

    if (!mxIsChar(prhs[5])) {
      dream_err_msg("Argument 6 must be a string");
    }

    buflen = (mxGetM(prhs[5]) * mxGetN(prhs[5]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[5],foc_str,buflen);

    if (!strcmp(foc_str,"off")) {
      foc_met = NO_FOCUS;
    }

    if (!strcmp(foc_str,"on")) {
      foc_met = FOCUS_XY;
    }

    if (!strcmp(foc_str,"ud")) {
      foc_met = FOCUS_UD;
      is_set = true;

      if (mxGetM(prhs[6]) * mxGetN(prhs[6]) != num_radii ) {
        printf("The time delay vector (argument 7) for user defined ('ud') focusing\n") ;
        dream_err_msg("delays must have the same length as the number of array elements.!");

      }
      focal = mxGetPr(prhs[6]);
    } else {

      // Check that arg 7 is a scalar.
      if (mxGetM(prhs[6]) * mxGetN(prhs[6]) !=1 ) {
        dream_err_msg("Argument 7 must be a scalar!");
      }

      // Focal point (in mm).
      focal = mxGetPr(prhs[6]);
    }

  } else {
    foc_met = NO_FOCUS;
  }

  //
  // Apodization.
  //

  if (nrhs >= 8) {

    if (!mxIsChar(prhs[7])) {
      dream_err_msg("Argument 8 must be a string");
    }

    buflen = (mxGetM(prhs[7]) * mxGetN(prhs[7]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[7],apod_str,buflen);

    do_apod = false;			// default off.
    is_set = false;

    if (!strcmp(apod_str,"off")) {
      do_apod = false;
      is_set = true;
    }

    if (!strcmp(apod_str,"ud")) {
      do_apod = true;
      apod_met = 0;
      is_set = true;

      // Vector of apodization weights.
      if (mxGetM(prhs[8]) * mxGetN(prhs[8]) != num_radii) {
        dream_err_msg("The length of argument 9 (apodization vector) must be the same as the number of array elements!");
      }

      apod = mxGetPr(prhs[8]);
    }

    if (!strcmp(apod_str,"triangle")) {
      do_apod = true;
      apod_met = APOD_TRIANGLE;
      is_set = true;
    }

    if (!strcmp(apod_str,"gauss")) {
      do_apod = true;
      apod_met = APOD_GAUSS;
      is_set = true;
    }

    if (!strcmp(apod_str,"raised")) {
      do_apod = true;
      apod_met = APOD_RISED_COSINE;
      is_set = true;
    }

    if (!strcmp(apod_str,"simply")) {
      do_apod = true;
      apod_met = APOD_SIMPLY_SUPPORTED;
      is_set = true;
    }

    if (!strcmp(apod_str,"clamped")) {
      do_apod = true;
      apod_met = APOD_CLAMPED;
      is_set = true;
    }

    if (is_set == false) {
      dream_err_msg("Unknown apodization method!");
    }

    // Parameter for raised cos and Gaussian apodization functions.
    if (mxGetM(prhs[9]) * mxGetN(prhs[9]) !=1 ) {
      dream_err_msg("Argument 10 must be a scalar");
    }

    param = mxGetScalar(prhs[9]);
  } else {
    do_apod = false;
  }

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11)
  nthreads = std::thread::hardware_concurrency();

  // Read DREAM_NUM_THREADS env var
  if(const char* env_p = std::getenv("DREAM_NUM_THREADS")) {
    unsigned int dream_threads = std::stoul(env_p);
    if (dream_threads < nthreads) {
      nthreads = dream_threads;
    }
  }

  // nthreads can't be larger then the number of observation points.
  if (nthreads > no) {
    nthreads = no;
  }

  //
  // Error Reporting.
  //

  if (nrhs == 11) {

    if (!mxIsChar(prhs[10])) {
      dream_err_msg("Argument 11 must be a string");
    }

    buflen = (mxGetM(prhs[10]) * mxGetN(prhs[10]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[10],err_str,buflen);

    is_set = false;

    if (!strcmp(err_str,"ignore")) {
      err_level = IGNORE;
      is_set = true;
    }

    if (!strcmp(err_str,"warn")) {
      err_level = WARN;
      is_set = true;
    }

    if (!strcmp(err_str,"stop")) {
      err_level = STOP;
      is_set = true;
    }

    if (is_set == false) {
      dream_err_msg("Unknown error level!");
    }
  } else {
    err_level = STOP; // Default.
  }

  // Create an output matrix for the impulse response
  plhs[0] = mxCreateDoubleMatrix(nt,no,mxREAL);
  h = mxGetPr(plhs[0]);

  //
  // Register signal handlers.
  //

  if ((old_handler = signal(SIGTERM, &sighandler)) == SIG_ERR) {
    printf("Couldn't register SIGTERM signal handler.\n");
  }

  if (( old_handler_abrt=signal(SIGABRT, &sighandler)) == SIG_ERR) {
    printf("Couldn't register SIGABRT signal handler.\n");
  }

  if (( old_handler_keyint=signal(SIGINT, &sighandler)) == SIG_ERR) {
    printf("Couldn't register SIGINT signal handler.\n");
  }

  //
  // Call the DREAM subroutine.
  //

  out_err = NONE;
  running = true;

  // Do we have attenuation?
  Attenuation att(nt, dt, alpha);
  Attenuation *att_ptr = nullptr;
  if (alpha > std::numeric_limits<double>::epsilon() ) {
    att_ptr = &att;

    // FIXME: Force to run in the main thread when we have nonzero attenuation
    // due to Matlab FFT threading issues.
    nthreads = 1;
  }

  // Allocate local data.
  D = (DATA*) malloc(nthreads*sizeof(DATA));

  // Allocate mem for the threads.
  threads = new std::thread[nthreads]; // Init thread data.

  for (thread_n = 0; thread_n < nthreads; thread_n++) {

    start = thread_n * no/nthreads;
    stop =  (thread_n+1) * no/nthreads;

    // Init local data.
    D[thread_n].start = start; // Local start index;
    D[thread_n].stop = stop; // Local stop index;
    D[thread_n].no = no;
    D[thread_n].ro = ro;
    D[thread_n].dx = dx;
    D[thread_n].dy = dy;
    D[thread_n].dt = dt;
    D[thread_n].nt = nt;

    if (mxGetM(prhs[3]) * mxGetN(prhs[3]) == 1)
      D[thread_n].delay_method = SINGLE_DELAY; // delay is a scalar.
    else
      D[thread_n].delay_method = MULTIPLE_DELAYS; // delay is a vector.

    D[thread_n].delay = delay;
    D[thread_n].v = v;
    D[thread_n].cp = cp;
    D[thread_n].att = att_ptr;
    D[thread_n].num_radii = num_radii;
    D[thread_n].gr = gr;
    D[thread_n].foc_met = foc_met;
    D[thread_n].do_apod = do_apod;
    D[thread_n].apod_met = apod_met;
    D[thread_n].focal = focal;
    D[thread_n].apod = apod;
    D[thread_n].param = param;
    D[thread_n].h = h;
    D[thread_n].err_level = err_level;

    if (nthreads > 1) {
      // Starts the threads.
      threads[thread_n] = std::thread(smp_dream_arr_annu, &D[thread_n]); // Start the threads.
      set_dream_thread_affinity(thread_n, nthreads, threads);
    } else {
      smp_dream_arr_annu(&D[0]);
    }

  }

  // Wait for all threads to finish.
  if (nthreads > 1) {
    for (thread_n = 0; thread_n < nthreads; thread_n++) {
      threads[thread_n].join();
    }
  }

  // Free memory.
  free((void*) D);

  //
  // Restore old signal handlers.
  //

  if (signal(SIGTERM, old_handler) == SIG_ERR) {
    printf("Couldn't register old SIGTERM signal handler.\n");
  }

  if (signal(SIGABRT,  old_handler_abrt) == SIG_ERR) {
    printf("Couldn't register old SIGABRT signal handler.\n");
  }

  if (signal(SIGINT, old_handler_keyint) == SIG_ERR) {
    printf("Couldn't register old SIGINT signal handler.\n");
  }

  if (!running) {
    dream_err_msg("CTRL-C pressed!\n"); // Bail out.
  }

  //
  // Check for Error.
  //

  if ( (err_level == STOP) && (out_err != NONE)) {
    dream_err_msg(""); // Bail out if error.
  }

  //
  // Return error.
  //
  if (nlhs == 2) {
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    err_p =  mxGetPr(plhs[1]);
    err_p[0] = (double) out_err;
  }

  return;
}
