/***
 *
 * Copyright (C) 2003,2006,2007,2008,2009,2014,2015,2016,2019,2021,2023 Fredrik Lingvall
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
#include <iostream>
#include <string>
#include <thread>

#include "affinity.h"
#include "circ_sir.h"
#include "arg_parser.h"

#include "mex.h"

//
// Globals
//

//volatile ErrorLevel out_err=ErrorLevel::none;
//std::mutex err_lock;
int running;

//
// typedef:s
//

typedef struct
{
  size_t No;
  size_t start;
  size_t stop;
  double *Ro;
  double r;
  double dt;
  size_t nt;
  DelayType delay_type;
  double *delay;
  double v;
  double cp;
  //double alpha;
  double *h;
  //ErrorLevel err_level;
} DATA;

typedef void (*sighandler_t)(int);

/***
 *
 * Thread function.
 *
 ***/

void* smp_dream_circ_sir(void *arg)
{
  //ErrorLevel tmp_err=ErrorLevel::none, err=ErrorLevel::none;
  DATA D = *(DATA *)arg;
  double xo, yo, zo;
  double *h = D.h;
  double r=D.r, dt=D.dt;
  size_t n, No=D.No, nt=D.nt;
  //int    tmp_lev, err_level=D.err_level;
  double *delay=D.delay, *Ro=D.Ro, v=D.v, cp=D.cp;// alpha=D.alpha;
  size_t start=D.start, stop=D.stop;

  // Let the thread finish and then catch the error.
  /*
  if (err_level == ErrorLevel::stop)
    tmp_lev = ErrorLevel::parallel_stop;
  else
    tmp_lev = err_level;
  */

  if (D.delay_type == DelayType::single) {
    for (n=start; n<stop; n++) {
      xo = Ro[n];
      yo = Ro[n+1*No];
      zo = Ro[n+2*No];

      circ_sir(xo,yo,zo,r,dt,nt,delay[0],v,cp,&h[n*nt]); // TODO: Add attenuation.

      if (!running) {
        std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!\n";
        return(NULL);
      }

    }
  } else { // DelayType::multiple.
    for (n=start; n<stop; n++) {
      xo = Ro[n];
      yo = Ro[n+1*No];
      zo = Ro[n+2*No];

      circ_sir(xo,yo,zo,r,dt,nt,delay[n],v,cp,&h[n*nt]); // TODO: Add attenuation.

      if (!running) {
        std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!\n";
        return(NULL);
      }

    }
  }

  // Lock out_err for update, update it, and unlock.
  /*
  err_lock.lock();

  if ((tmp_err != ErrorLevel::none) && (out_err == ErrorLevel::none))
    out_err = tmp_err;

  err_lock.unlock();
  */

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
 * Matlab (MEX) gateway function for circ_sir.
 *
 ***/

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *Ro,*geom_par, *s_par, *m_par;
  dream_idx_type nt, No;
  double r, dt;
  double *delay, v, cp;
  double *h;
  DATA   *D;
  size_t start, stop;
  std::thread *threads;
  dream_idx_type thread_n, nthreads;
  sighandler_t  old_handler, old_handler_abrt, old_handler_keyint;

  ArgParser ap;

  // Check for proper number of arguments

  ap.check_arg_in("circ_sir", nrhs, 5, 5);
  ap.check_arg_out("circ_sir", nlhs, 0, 1);

  //
  // Observation point.
  //

  ap.check_obs_points("circ_sir", prhs, 0);
  No = mxGetM(prhs[0]); // Number of observation points.
  Ro = mxGetPr(prhs[0]);

  //
  // Transducer geometry
  //

  ap.check_geometry("circ_sir", prhs, 1, 1);
  geom_par = mxGetPr(prhs[1]);
  r  = geom_par[0];		// Radius of the transducer.

  //
  // Temporal and spatial sampling parameters.
  //

  ap.check_sampling("circ_sir", prhs, 2, 2);
  s_par = mxGetPr(prhs[2]);
  dt    = s_par[0]; // Temporal discretization size (= 1/sampling freq).
  nt    = (size_t) s_par[1]; // Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //

  ap.check_delay("circ_sir", prhs, 3, No);
  delay = mxGetPr(prhs[3]);

  //
  // Material parameters
  //

  ap.check_material("circ_sir", prhs, 4, 2);
  m_par = mxGetPr(prhs[4]);
  v     = m_par[0]; // Normal velocity of transducer surface.
  cp    = m_par[1]; // Sound speed.

  //
  // Number of threads
  //

  // Get number of CPU cores (including hypethreading, C++11)
  nthreads = std::thread::hardware_concurrency();

  // Read DREAM_NUM_THREADS env var
  if(const char* env_p = std::getenv("DREAM_NUM_THREADS")) {
    dream_idx_type dream_threads = (dream_idx_type) std::stoul(env_p);
    if (dream_threads < nthreads) {
      nthreads = dream_threads;
    }
  }

  // nthreads can't be larger then the number of observation points.
  if (nthreads > No) {
    nthreads = No;
  }

  //
  // Create an output matrix for the impulse response(s).
  //

  plhs[0] = mxCreateDoubleMatrix(nt,No,mxREAL);
  h = mxGetPr(plhs[0]);

  SIRData hsir(h, nt, No);
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
  // Call the analytic SIR subroutine.
  //

  //out_err = ErrorLevel::none;
  running = true;

  // Allocate local data.
  D = (DATA*) malloc(nthreads*sizeof(DATA));

  // Allocate mem for the threads.
  threads = new std::thread[nthreads]; // Init thread data.

  for (thread_n = 0; thread_n < nthreads; thread_n++) {

    start = thread_n * No/nthreads;
    stop =  (thread_n+1) * No/nthreads;

    // Init local data.
    D[thread_n].start = start; // Local start index;
    D[thread_n].stop = stop; // Local stop index;
    D[thread_n].No = No;
    D[thread_n].Ro = Ro;
    D[thread_n].r = r;
    D[thread_n].dt = dt;
    D[thread_n].nt = nt;

    if (mxGetM(prhs[3]) * mxGetN(prhs[3]) == 1)
      D[thread_n].delay_type = DelayType::single; // delay is a scalar.
    else
      D[thread_n].delay_type = DelayType::multiple; // delay is a vector.

    D[thread_n].delay = delay;
    D[thread_n].v = v;
    D[thread_n].cp = cp;
    //D[thread_n].alpha = alpha;
    D[thread_n].h = h;
    //D[thread_n].err_level = err_level;

    // Start the threads.
    threads[thread_n] = std::thread(smp_dream_circ_sir, &D[thread_n]); // Start the threads.
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
    dream_err_msg("CTRL-C pressed!\n"); // Bail out.
  }

  return;
}
