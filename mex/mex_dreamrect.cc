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

#include <csignal>
#include <thread>
#include <mutex>

#include "dreamrect.h"
#include "affinity.h"
#include "arg_parser.h"

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
  size_t no;
  size_t start;
  size_t stop;
  double *ro;
  double a;
  double b;
  double dx;
  double dy;
  double dt;
  size_t nt;
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

void* smp_dream_rect(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

/***
 *
 * Thread function.
 *
 ***/

void* smp_dream_rect(void *arg)
{
  ErrorLevel tmp_err=ErrorLevel::none, err=ErrorLevel::none;
  DATA D = *(DATA *)arg;
  double xo, yo, zo;
  double *h = D.h;
  double a=D.a, b=D.b, dx=D.dx, dy=D.dy, dt=D.dt;
  size_t n, no=D.no, nt=D.nt;
  ErrorLevel tmp_lev=ErrorLevel::none, err_level=D.err_level;
  double *delay=D.delay, *ro=D.ro, v=D.v, cp=D.cp;
  Attenuation *att = D.att;
  size_t start=D.start, stop=D.stop;

  // Buffers for the FFTs in the Attenuation
  std::unique_ptr<FFTCVec> xc_vec;
  std::unique_ptr<FFTVec> x_vec;
  if (att) {
    xc_vec = std::make_unique<FFTCVec>(nt);
    x_vec = std::make_unique<FFTVec>(nt);
  }

  // Let the thread finish and then catch the error.
  if (err_level == ErrorLevel::stop)
    tmp_lev = ErrorLevel::parallel_stop;
  else
    tmp_lev = err_level;

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
      err = dreamrect(xo, yo, zo,
                      a, b,
                      dx, dy,
                      dt, nt, dlay, v, cp,
                      &h[n*nt], tmp_lev);
    } else {
      err = dreamrect(*att, *xc_vec, *x_vec,
                      xo, yo, zo,
                      a, b,
                      dx, dy,
                      dt, nt, dlay, v, cp,
                      &h[n*nt], tmp_lev);
    }

    if (err != ErrorLevel::none || out_err ==  ErrorLevel::parallel_stop) {
      tmp_err = err;
      if (err == ErrorLevel::parallel_stop || out_err ==  ErrorLevel::parallel_stop)
        break; // Jump out when a ErrorLevel::stop error occurs.
    }

    if (!running) {
      std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!" << std::endl;
      return(NULL);
    }

  }

  // Lock out_err for update, update it, and unlock.
  err_lock.lock();

  if ((tmp_err != ErrorLevel::none) && (out_err == ErrorLevel::none))
    out_err = tmp_err;

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
 * mex_dreamrect.cc -  Matlab (MEX) gateway function for (parallel) dreamrect.
 *
 ***/

extern void _main();

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  ErrorLevel err_level=ErrorLevel::stop;
  DATA *D;
  std::thread *threads;
  sighandler_t old_handler, old_handler_abrt, old_handler_keyint;

  ArgParser ap;

  // Check for proper number of arguments

  ap.check_arg_in("dreamrect", nrhs, 5, 6);
  ap.check_arg_out("dreamrect", nlhs, 0, 2);

  //
  // Observation point.
  //

  ap.check_obs_points("dreamrect", prhs, 0);
  dream_idx_type no = mxGetM(prhs[0]); // Number of observation points.
  double *ro = mxGetPr(prhs[0]);

  //
  // Transducer geometry
  //

  double a=0.0, b=0.0, dummy=0.0;
  ap.parse_geometry("dreamrect", prhs, 1, 2, a, b, dummy);

  //
  // Temporal and spatial sampling parameters.
  //

  double dx=0.0, dy=0.0, dt=0.0;
  dream_idx_type nt=0;
  ap.parse_sampling("dreamect", prhs, 2, 4, dx, dy, dt, nt);

  //
  // Start point of impulse response vector ([us]).
  //

  ap.check_delay("dreamrect", prhs, 3, no);
  double *delay = mxGetPr(prhs[3]);

  //
  // Material parameters
  //

  double v=1.0, cp=1000.0, alpha=0.0;
  ap.parse_material("dreamrect", prhs, 4, v, cp, alpha);

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

  // nthreads can't be larger then the number of observation points.
  if (nthreads > no) {
    nthreads = no;
  }

  //
  // Error reporting.
  //

  if (nrhs == 6) {
    ap.parse_error_arg("dreamline", prhs, 5, err_level);
  } else {
    err_level = ErrorLevel::stop; // Default.
  }

  // Create an output matrix for the impulse response
  plhs[0] = mxCreateDoubleMatrix(nt,no,mxREAL);
  double *h = mxGetPr(plhs[0]);

  SIRData hsir(h, nt, no);
  hsir.clear();

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
  // Call the DREAM subroutine.
  //

  out_err = ErrorLevel::none;
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

  for (dream_idx_type thread_n=0; thread_n < nthreads; thread_n++) {

    dream_idx_type start = thread_n * no/nthreads;
    dream_idx_type stop =  (thread_n+1) * no/nthreads;

    // Init local data.
    D[thread_n].start = start; // Local start index;
    D[thread_n].stop = stop; // Local stop index;
    D[thread_n].no = no;
    D[thread_n].ro = ro;
    D[thread_n].a = a;
    D[thread_n].b = b;
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
      threads[thread_n] = std::thread(smp_dream_rect, &D[thread_n]);
      set_dream_thread_affinity(thread_n, nthreads, threads);
    } else {
      smp_dream_rect(&D[0]);
    }
  }

  // Wait for all threads to finish.
  if (nthreads > 1) {
    for (dream_idx_type thread_n=0; thread_n < nthreads; thread_n++) {
      threads[thread_n].join();
    }
  }

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
    double *err_p =  mxGetPr(plhs[1]);
    err_p[0] = (double) out_err;
  }

  return;
}
