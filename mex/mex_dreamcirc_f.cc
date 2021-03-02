/***
 *
 * Copyright (C) 2003,2006,2007,2008,2009,2014,2015,2019,2021 Fredrik Lingvall
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

#include "affinity.h"
#include "dreamcirc_f.h"
#include "dream_error.h"

#include "mex.h"

//
// Globals
//

volatile ErrorLevel out_err=ErrorLevel::none;
std::mutex err_lock;
int running;

//
// typedef:s
//

typedef struct
{
  dream_idx_type no;
  dream_idx_type start;
  dream_idx_type stop;
  double *ro;
  double R;
  FocusMet foc_met;
  double focal;
  double dx;
  double dy;
  double dt;
  dream_idx_type nt;
  DelayType delay_type;
  double *delay;
  double v;
  double cp;
  Attenuation *att;
  double *h;
  ErrorLevel err_level;
} DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//

void* smp_dream_circ_f(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

/***
 *
 * Thread function.
 *
 ***/

void* smp_dream_circ_f(void *arg)
{
  ErrorLevel tmp_err=ErrorLevel::none, err=ErrorLevel::none;
  DATA D = *(DATA *)arg;
  double xo, yo, zo;
  double *h = D.h;
  double R=D.R, dx=D.dx, dy=D.dy, dt=D.dt;
  dream_idx_type n, no=D.no, nt=D.nt;
  ErrorLevel tmp_lev=ErrorLevel::none, err_level=D.err_level;
  double *delay=D.delay, *ro=D.ro, v=D.v, cp=D.cp, focal=D.focal;
  Attenuation *att = D.att;
  dream_idx_type start=D.start, stop=D.stop;
  FocusMet foc_met=D.foc_met;

  // Buffers for the FFTs in the Attenuation
  std::unique_ptr<FFTCVec> xc_vec;
  std::unique_ptr<FFTVec> x_vec;
  if (att) {
    xc_vec = std::make_unique<FFTCVec>(nt);
    x_vec = std::make_unique<FFTVec>(nt);
  }

  // Let the thread finish and then catch the error.
  if (err_level == ErrorLevel::stop) {
    tmp_lev = ErrorLevel::parallel_stop;
  } else {
    tmp_lev = err_level;
  }

  for (n=start; n<stop; n++) {
    xo = ro[n];
    yo = ro[n+1*no];
    zo = ro[n+2*no];

    double dlay = 0.0;
    if (D.delay_type == DelayType::single) {
      dlay = delay[0];
    } else { // DelayType::multiple.
      dlay = delay[n];
    }

    if (att == nullptr) {
      err = dreamcirc_f(xo, yo, zo,
                        R, foc_met, focal,
                        dx, dy, dt,
                        nt,dlay, v, cp,
                        &h[n*nt], tmp_lev);
    } else {
      err = dreamcirc_f(*att, *xc_vec, *x_vec,
                        xo, yo, zo,
                        R, foc_met, focal,
                        dx, dy, dt,
                        nt, dlay, v, cp,
                        &h[n*nt], tmp_lev);
    }

    if (err != ErrorLevel::none || out_err ==  ErrorLevel::parallel_stop) {
      tmp_err = err;
      if (err == ErrorLevel::parallel_stop || out_err ==  ErrorLevel::parallel_stop) {
        break; // Jump out when a ErrorLevel::stop error occurs.
      }
    }

    if (!running) {
      std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!" << std::endl;
      return(NULL);
    }

  }

  // Lock out_err for update, update it, and unlock.
  err_lock.lock();

  if ((tmp_err != ErrorLevel::none) && (out_err == ErrorLevel::none)) {
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
 * Matlab (MEX) gateway function for dreamcirc_f.
 *
 ***/

extern void _main();

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *ro, *geom_par, *s_par, *m_par;
  dream_idx_type nt, no;
  dream_idx_type buflen;
  double R, dx, dy, dt;
  double *delay=nullptr, v, cp, alpha;
  FocusMet foc_met=FocusMet::none;
  char   foc_str[50];
  double focal=0.0;
  double *h, *err_p;
  ErrorLevel err_level=ErrorLevel::stop;
  bool is_set=false;
  char err_str[50];
  DATA *D;
  dream_idx_type start, stop;
  std::thread *threads;
  unsigned int thread_n, nthreads;
  sighandler_t old_handler, old_handler_abrt, old_handler_keyint;

  // Check for proper number of arguments

  if (!((nrhs == 7) || (nrhs == 8))) {
    dream_err_msg("dreamcirc_f requires 7 or 8 input arguments!");
  } else {
    if (nlhs > 2) {
      dream_err_msg("Too many output arguments for dreamcirc_f!");
    }
  }

  //
  // Observation point.
  //

  // Check that arg (number of observation points) x 3 matrix
  if (mxGetN(prhs[0]) != 3) {
    dream_err_msg("Argument 1 must be a (number of observation points) x 3 matrix!");
  }

  no = mxGetM(prhs[0]); // Number of observation points.
  ro = mxGetPr(prhs[0]);

  //
  // Transducer geometry
  //

  // Check that arg 2 is a scalar.
  if (!((mxGetM(prhs[1])==1 && mxGetN(prhs[1])==1))) {
    dream_err_msg("Argument 2 must be a scalar!");
  }

  geom_par = mxGetPr(prhs[1]);
  R = geom_par[0];              // Radius of the transducer.

  //
  // Temporal and spatial sampling parameters.
  //

  // Check that arg 3 is a 4 element vector
  if (!((mxGetM(prhs[2])==4 && mxGetN(prhs[2])==1) || (mxGetM(prhs[2])==1 && mxGetN(prhs[2])==4))) {
    dream_err_msg("Argument 3 must be a vector of length 4!");
  }

  s_par = mxGetPr(prhs[2]);
  dx = s_par[0];  // Spatial x discretization size.
  dy = s_par[1];  // Spatial dy iscretization size.
  dt = s_par[2];  // Temporal discretization size (= 1/sampling freq).
  nt = (dream_idx_type) s_par[3]; // Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //

  // Check that arg 4 is a scalar or vector
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
  alpha  = m_par[2]; // Attenuation coefficient [dB/(cm MHz)].

  //
  // Focusing parameters.
  //

  if (nrhs >= 6) {

    if (!mxIsChar(prhs[5])) {
      dream_err_msg("Argument 6 must be a string");
    }

    buflen = (mxGetM(prhs[5]) * mxGetN(prhs[5]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[5],foc_str,buflen);

    is_set = false;

    if (!strcmp(foc_str,"off")) {
      foc_met = FocusMet::none;
      is_set = true;
    }

    if (!strcmp(foc_str,"x")) {
      foc_met = FocusMet::x;
      is_set = true;
    }

    if (!strcmp(foc_str,"y")) {
      foc_met = FocusMet::y;
      is_set = true;
    }

    if (!strcmp(foc_str,"xy")) {
      foc_met = FocusMet::xy;
      is_set = true;
    }

    if (!strcmp(foc_str,"x+y")) {
      foc_met = FocusMet::x_y;
      is_set = true;
    }

    if (is_set == false) {
      dream_err_msg("Unknown focusing method!");
    }

    // Check that arg 7 is a scalar.
    if (mxGetM(prhs[6]) * mxGetN(prhs[6]) != 1 ) {
      dream_err_msg("Argument 7 must be a scalar!");
    }

    // Focal point (in mm).
    focal = mxGetScalar(prhs[6]);

  } else {
    foc_met = FocusMet::none;
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
  // Error reporting.
  //

  if (nrhs == 8) {

    if (!mxIsChar(prhs[7])) {
      dream_err_msg("Argument 8 must be a string");
    }

    buflen = (mxGetM(prhs[7]) * mxGetN(prhs[7]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[7],err_str,buflen);

    is_set = false;

    if (!strcmp(err_str,"ignore")) {
      err_level = ErrorLevel::ignore;
      is_set = true;
    }

    if (!strcmp(err_str,"warn")) {
      err_level = ErrorLevel::warn;
      is_set = true;
    }

    if (!strcmp(err_str,"stop")) {
      err_level = ErrorLevel::stop;
      is_set = true;
    }

    if (is_set == false) {
      dream_err_msg("Unknown error level!");
    }
  } else {
    err_level = ErrorLevel::stop; // Default.
  }

  // Create an output matrix for the impulse response
  plhs[0] = mxCreateDoubleMatrix(nt,no,mxREAL);
  h = mxGetPr(plhs[0]);

  //
  // Register signal handlers.
  //

  if ((old_handler = signal(SIGTERM, &sighandler)) == SIG_ERR) {
    std::cerr << "Couldn't register SIGTERM signal handler!" << std::endl;
  }

  if ((old_handler_abrt=signal(SIGABRT, &sighandler)) == SIG_ERR) {
    std::cerr << "Couldn't register SIGABRT signal handler!" << std::endl;
  }

  if ((old_handler_keyint=signal(SIGINT, &sighandler)) == SIG_ERR) {
    std::cerr << "Couldn't register SIGINT signal handler!" << std::endl;
  }

  //
  // Call the DREAM subroutine.
  //

  out_err = ErrorLevel::none;
  running=true;

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
    D[thread_n].R = R;
    D[thread_n].foc_met = foc_met;
    D[thread_n].focal = focal;
    D[thread_n].dx = dx;
    D[thread_n].dy = dy;
    D[thread_n].dt = dt;
    D[thread_n].nt = nt;

    if (mxGetM(prhs[3]) * mxGetN(prhs[3]) == 1) {
      D[thread_n].delay_type = DelayType::single; // delay is a scalar.
    } else {
      D[thread_n].delay_type = DelayType::multiple; // delay is a vector.
    }

    D[thread_n].delay = delay;
    D[thread_n].v = v;
    D[thread_n].cp = cp;
    D[thread_n].att = att_ptr;
    D[thread_n].h = h;
    D[thread_n].err_level = err_level;

    if (nthreads > 1) {
      // Starts the threads.
      threads[thread_n] = std::thread(smp_dream_circ_f, &D[thread_n]); // Start the threads.
      set_dream_thread_affinity(thread_n, nthreads, threads);
    } else {
      smp_dream_circ_f(&D[0]);
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
    std::cerr << "Couldn't register old SIGTERM signal handler!" << std::endl;
  }

  if (signal(SIGABRT,  old_handler_abrt) == SIG_ERR) {
    std::cerr << "Couldn't register old SIGABRT signal handler!" << std::endl;
  }

  if (signal(SIGINT, old_handler_keyint) == SIG_ERR) {
    std::cerr << "Couldn't register old SIGINT signal handler!" << std::endl;
  }

  if (!running) {
    dream_err_msg("CTRL-C pressed!\n"); // Bail out.
  }

  //
  // Check for Error.
  //

  if ( (err_level == ErrorLevel::stop) && (out_err != ErrorLevel::none)) {
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
