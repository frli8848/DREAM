/***
*
* Copyright (C) 2003,2006,2007,2008,2009,2014,2015 Fredrik Lingvall
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
#include <iostream>
#include <thread>
#include <mutex>
#include <signal.h>
#include "mex.h"
#include "dreamcirc_f.h"
#include "dream_error.h"

#define SINGLE 0
#define MULTIPLE 1

#ifdef USE_FFTW
#include "att.h"
#endif

// Parallel implementation.

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
  dream_idx_type no;
  dream_idx_type start;
  dream_idx_type stop;
  double *RESTRICT ro;
  double r;
  double dx;
  double dy;
  double dt;
  dream_idx_type nt;
  int delay_method;
  double *RESTRICT delay;
  double v;
  double cp;
  double alfa;
  int ifoc;
  double focal;
  double *RESTRICT h;
  int err_level;
} DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//
void* smp_process(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

/***
 *
 * Thread function.
 *
 ***/

void* smp_process(void *arg)
{
  int tmp_err = NONE, err = NONE;
  DATA D = *(DATA *)arg;
  double xo, yo, zo;
  double *RESTRICT h = D.h;
  double r=D.r, dx=D.dx, dy=D.dy, dt=D.dt;
  size_t n, no=D.no, nt=D.nt;
  int    tmp_lev, err_level=D.err_level;
  double *RESTRICT delay=D.delay, *ro=D.ro, v=D.v, cp=D.cp, alfa=D.alfa, focal=D.focal;
  size_t  start=D.start, stop=D.stop;
  int ifoc = D.ifoc;

  // Let the thread finish and then catch the error.
  if (err_level == STOP)
    tmp_lev = PARALLEL_STOP;
  else
    tmp_lev = err_level;

  if (D.delay_method == SINGLE) {
    for (n=start; n<stop; n++) {
      xo = ro[n];
      yo = ro[n+1*no];
      zo = ro[n+2*no];
      err = dreamcirc_f(xo,yo,zo,r,dx,dy,dt,nt,delay[0],v,cp,alfa,
                      ifoc,focal,&h[n*nt],tmp_lev);

      if (err != NONE || out_err ==  PARALLEL_STOP) {
        tmp_err = err;
        if (err == PARALLEL_STOP || out_err ==  PARALLEL_STOP)
          break; // Jump out when a STOP error occurs.
      }

      if (!running) {
        std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!" << std::endl;
        return(NULL);
      }

    }
  } else { // MULTIPLE delays.
    for (n=start; n<stop; n++) {
      xo = ro[n];
      yo = ro[n+1*no];
      zo = ro[n+2*no];
      err = dreamcirc_f(xo,yo,zo,r,dx,dy,dt,nt,delay[n],v,cp,alfa,
                      ifoc,focal,&h[n*nt],tmp_lev);

      if (err != NONE || out_err ==  PARALLEL_STOP) {
        tmp_err = err;
        if (err == PARALLEL_STOP || out_err ==  PARALLEL_STOP)
          break; // Jump out when a STOP error occurs.
      }

      if (!running) {
        std::cout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!" << std::endl;
        return(NULL);
      }

    }
  }

  // Lock out_err for update, update it, and unlock.
  err_lock.lock();

  if ((tmp_err != NONE) && (out_err == NONE))
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
 * Matlab (MEX) gateway function for parallel dreamcirc_f.
 *
 ***/

extern void _main();

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *RESTRICT ro, *RESTRICT geom_par, *RESTRICT s_par, *RESTRICT m_par;
  size_t nt, no;
  int    ifoc=0;
  char   foc_met[50];
  int    buflen;
  double r, dx, dy, dt;
  double *RESTRICT delay, v, cp, alfa,focal=0;
  double *RESTRICT h, *err_p;
  int    err_level=STOP, set = false;
  char   err_str[50];
  DATA   *D;
  size_t start, stop;
  std::thread *threads;
  unsigned int thread_n, nthreads;
  sighandler_t old_handler, old_handler_abrt, old_handler_keyint;

  // Check for proper number of arguments

  if (!((nrhs == 7) || (nrhs == 8))) {
    dream_err_msg("dreamcirc_f requires 7 or 8 input arguments!");
  }
  else
    if (nlhs > 2) {
      dream_err_msg("Too many output arguments for dreamcirc_f!");
    }

  //
  // Observation point.
  //

  // Check that arg (number of observation points) x 3 matrix
  if (mxGetN(prhs[0]) != 3)
    dream_err_msg("Argument 1 must be a (number of observation points) x 3 matrix!");

  no = mxGetM(prhs[0]); // Number of observation points.
  ro = mxGetPr(prhs[0]);

  //
  // Transducer geometry
  //

  // Check that arg 2 is a scalar.
  if (!((mxGetM(prhs[1])==1 && mxGetN(prhs[1])==1)))
    dream_err_msg("Argument 2 must be a scalar!");

  geom_par = mxGetPr(prhs[1]);
  r  = geom_par[0];		// Radius of the transducer.


  //
  // Temporal and spatial sampling parameters.
  //

  // Check that arg 3 is a 4 element vector
  if (!((mxGetM(prhs[2])==4 && mxGetN(prhs[2])==1) || (mxGetM(prhs[2])==1 && mxGetN(prhs[2])==4)))
    dream_err_msg("Argument 3 must be a vector of length 4!");

  s_par = mxGetPr(prhs[2]);
  dx    = s_par[0];		// Spatial x discretization size.
  dy    = s_par[1];		// Spatial dy iscretization size.
  dt    = s_par[2];		// Temporal discretization size (= 1/sampling freq).
  nt    = (dream_idx_type) s_par[3];	// Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //

  // Check that arg 4 is a scalar or vector
  if ( (mxGetM(prhs[3]) * mxGetN(prhs[3]) !=1) && ((mxGetM(prhs[3]) * mxGetN(prhs[3])) != no))
    dream_err_msg("Argument 4 must be a scalar or a vector with a length equal to the number of observation points!");

  delay = mxGetPr(prhs[3]);

  //
  // Material parameters
  //

  // Check that arg 5 is a 3 element vectora
  if (!((mxGetM(prhs[4])==3 && mxGetN(prhs[4])==1) || (mxGetM(prhs[4])==1 && mxGetN(prhs[4])==3)))
    dream_err_msg("Argument 5 must be a vector of length 3!");

  m_par = mxGetPr(prhs[4]);
  v     = m_par[0]; // Normal velocity of transducer surface.
  cp    = m_par[1]; // Sound speed.
  alfa  = m_par[2]; // Attenuation coefficient [dB/(cm MHz)].

  //
  // Focusing parameters.
  //

  //  ifoc = 1 - no foc, 2 foc x ,3 foc y, 4 foc xy (del=fsqrt(x*x+y*y)), 5 focx+focy.

  if (nrhs >= 6) {

    if (!mxIsChar(prhs[5]))
      dream_err_msg("Argument 6 must be a string");

    buflen = (mxGetM(prhs[5]) * mxGetN(prhs[5]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[5],foc_met,buflen);

    set = false;

    if (!strcmp(foc_met,"off")) {
      ifoc = 1;
      set = true;
    }

    if (!strcmp(foc_met,"x")) {
      ifoc = 2;
      set = true;
    }

    if (!strcmp(foc_met,"y")) {
      ifoc = 3;
      set = true;
    }

    if (!strcmp(foc_met,"xy")) {
      ifoc = 4;
      set = true;
    }

    if (!strcmp(foc_met,"x+y")) {
      ifoc = 5;
      set = true;
    }

   if (set == false)
      dream_err_msg("Unknown focusing method!");

    // Check that arg 7 is a scalar.
    if (mxGetM(prhs[6]) * mxGetN(prhs[6]) !=1 )
      dream_err_msg("Argument 7 must be a scalar!");

    // Focal point (in mm).
    focal = mxGetScalar(prhs[6]);

  } else
    ifoc = 1;

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11)
  nthreads = std::thread::hardware_concurrency();

  // Read OMP_NUM_THREADS env var
  if(const char* env_p = std::getenv("OMP_NUM_THREADS")) {
    unsigned int omp_threads = std::stoul(env_p);
    if (omp_threads < nthreads) {
      nthreads = omp_threads;
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

    if (!mxIsChar(prhs[7]))
      dream_err_msg("Argument 8 must be a string");

    buflen = (mxGetM(prhs[7]) * mxGetN(prhs[7]) * sizeof(mxChar)) + 1;
    mxGetString(prhs[7],err_str,buflen);

    set = false;

    if (!strcmp(err_str,"ignore")) {
      err_level = IGNORE;
      set = true;
    }

    if (!strcmp(err_str,"warn")) {
      err_level = WARN;
      set = true;
    }

    if (!strcmp(err_str,"stop")) {
      err_level = STOP;
      set = true;
    }

    if (set == false)
      dream_err_msg("Unknown error level!");

  }
  else
    err_level = STOP; // Default.

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
  running=true;

#ifdef USE_FFTW
  if (alfa != (double) 0.0)
    att_init(nt,nthreads);
#endif

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
    D[thread_n].r = r;
    D[thread_n].dx = dx;
    D[thread_n].dy = dy;
    D[thread_n].dt = dt;
    D[thread_n].nt = nt;

    if (mxGetM(prhs[3]) * mxGetN(prhs[3]) == 1)
      D[thread_n].delay_method = SINGLE; // delay is a scalar.
    else
      D[thread_n].delay_method = MULTIPLE; // delay is a vector.

    D[thread_n].delay = delay;
    D[thread_n].v = v;
    D[thread_n].cp = cp;
    D[thread_n].alfa = alfa;
    D[thread_n].ifoc = ifoc;
    D[thread_n].focal = focal;
    D[thread_n].h = h;
    D[thread_n].err_level = err_level;

    // Starts the threads.
    threads[thread_n] = std::thread(smp_process, &D[thread_n]); // Start the threads.
  }

  // Wait for all threads to finish.
  for (thread_n = 0; thread_n < nthreads; thread_n++)
    threads[thread_n].join();

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
    printf("Couldn't register old  SIGINT signal handler.\n");
  }

#ifdef USE_FFTW
  if (alfa != (double) 0.0)
    att_close();
#endif

  if (!running) {
    dream_err_msg("CTRL-C pressed!\n"); // Bail out.
  }

  //
  // Check for Error.
  //

  if ( (err_level == STOP) && (out_err != NONE))
    dream_err_msg(""); // Bail out if error.

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
