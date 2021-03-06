/***
*
* Copyright (C) 2008,2009,2012,2015,2016,2021,2020 Fredrik Lingvall
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

#include <thread>
#include <mutex>

#include "rect_sir.h"
#include "affinity.h"
#include "dream_error.h"

//
// Octave headers.
//

#include <octave/oct.h>

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
  octave_idx_type no;
  octave_idx_type start;
  octave_idx_type stop;
  double *ro;
  double a;
  double b;
  double dt;
  octave_idx_type nt;
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

void* smp_dream_rect_sir(void *arg)
{
  //ErrorLevel tmp_err=ErrorLevel::none, err=ErrorLevel::none;
  DATA D = *(DATA *)arg;
  double xo, yo, zo;
  double *h = D.h;
  double a=D.a, b=D.b, dt=D.dt;
  octave_idx_type n, no=D.no, nt=D.nt;
  //int    tmp_lev, err_level=D.err_level;
  double *delay=D.delay, *ro=D.ro, v=D.v, cp=D.cp; // alfa=D.alfa;
  octave_idx_type start=D.start, stop=D.stop;

  // Let the thread finish and then catch the error.
  /*
  if (err_level == ErrorLevel::stop)
    tmp_lev = ErrorLevel::parallel_stop;
  else
    tmp_lev = err_level;
  */

  if (D.delay_type == DelayType::single) {
    for (n=start; n<stop; n++) {
      xo = ro[n];
      yo = ro[n+1*no];
      zo = ro[n+2*no];

      rect_sir(xo,yo,zo,a,b,dt,nt,delay[0],v,cp,&h[n*nt]); // TODO: Add attenuation.

      if (!running) {
        octave_stdout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!\n";
        return(NULL);
      }

    }
  } else { // DelayType::multiple.
    for (n=start; n<stop; n++) {
      xo = ro[n];
      yo = ro[n+1*no];
      zo = ro[n+2*no];

      rect_sir(xo,yo,zo,a,b,dt,nt,delay[n],v,cp,&h[n*nt]); // TODO: Add attenuation.

      if (!running) {
        octave_stdout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!\n";
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
 *  rect_sir - Octave (oct) gateway function for RECT_SIR.
 *
 ***/

DEFUN_DLD (rect_sir, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {}  [Y] = rect_sir(Ro,geom_par,delay,s_par,m_par,device).\n \
\n\
RECT_SIR Computes the time-continous (analytic) spatial impulse response(s) for a\n\
rectangular transducer.\n\
\n\
Observation point(s) ([mm]):\n\
\n\
@table @code\n\
@item Ro\n\
An N x 3 matrix, Ro = [xo1 yo1 zo2; xo2 yo2 zo2; ... xoN yoN zoN]; where N is the number of observation points.\n\
@end table\n\
\n\
Geometrical parameters: geom_par = [a b];\n\
\n\
@table @code\n\
@item a\n\
@item b\n\
Dimensions of the transducer [mm].\n\
@end table\n\
\n\
Sampling parameters: s_par = [dt nt]; \n\
\n\
@table @code\n\
@item dt\n\
Temporal discretization period (= 1/sampling freq) [us].\n\
@item  nt\n\
Length of impulse response vector.\n\
@end table\n\
\n\
Start point of SIR:\n\
\n\
@table @code\n\
@item  delay\n\
Scalar delay for all observation points or a vector with individual delays for each observation point [us].\n\
@end table\n\
\n\
Material parameters: m_par = [v cp];\n\
\n\
@table @code\n\
@item v\n\
Normal velocity [m/s].\n\
@item cp\n\
Sound velocity [m/s].\n\
@end table\n\
\n\
Compute device:\n\
\n\
@table @code\n\
@item 'device'\n\
A string which can be one of 'cpu' or 'gpu'.\n\
@end table\n\
\n\
rect_sir is an oct-function that is a part of the DREAM Toolbox available at\n\
@url{http://www.signal.uu.se/Toolbox/dream/}.\n\
\n\
Copyright @copyright{} 2008-2020 Fredrik Lingvall.\n\
@seealso {dreamrect,circ_sir}\n\
@end deftypefn")
{
  double *ro,*geom_par,*s_par,*m_par;
  dream_idx_type nt, no;
  double  a, b, dt;
  double *delay,v,cp;
  double *h;
  DATA   *D;
  octave_idx_type start, stop;
  std::thread *threads;
  unsigned int thread_n, nthreads;
  sighandler_t  old_handler, old_handler_abrt, old_handler_keyint;
  std::string device;
  octave_value_list oct_retval;

  int nrhs = args.length ();

  // Check for proper number of arguments

  if (nrhs != 5 && nrhs != 6) {
    dream_err_msg("rect_sir requires 5 or 6 input arguments!");
  }
  else
    if (nlhs > 1) {
      dream_err_msg("Too many output arguments for rect_sir!");
    }

  //
  // Observation point.
  //

  // Check that arg (number of observation points) x 3 matrix
  if (mxGetN(0) != 3) {
    error("Argument 1 must be a (number of observation points) x 3 matrix!");
    return oct_retval;
  }

  no = mxGetM(0); // Number of observation points.
  const Matrix tmp0 = args(0).matrix_value();
  ro = (double*) tmp0.fortran_vec();

  //
  // Transducer geometry
  //

  // Check that arg 2 is a two element vector.
  if (!((mxGetM(1)==2 && mxGetN(1)==1) || (mxGetM(1)==1 && mxGetN(1)==2))) {
    error("Argument 2 must be a two element vector!");
    return oct_retval;
  }

  const Matrix tmp1 = args(1).matrix_value();
  geom_par = (double*) tmp1.fortran_vec();
  a  = geom_par[0];		// x-dim of the transducer.
  b  = geom_par[1];		// y-dim of the transducer.

  //
  // Temporal and spatial sampling parameters.
  //

  // Check that arg 3 is a 2 element vector
  if (!((mxGetM(2)==2 && mxGetN(2)==1) || (mxGetM(2)==1 && mxGetN(2)==2))) {
    error("Argument 3 must be a vector of length 2!");
    return oct_retval;
  }

  const Matrix tmp2 = args(2).matrix_value();
  s_par = (double*) tmp2.fortran_vec();
  dt    = s_par[0]; // Temporal discretization size (= 1/sampling freq).
  nt    = (dream_idx_type) s_par[1];	// Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //

  // Check that arg 4 is a scalar (or vector).
  if ( (mxGetM(3) * mxGetN(3) !=1) && ((mxGetM(3) * mxGetN(3)) != no)) {
    error("Argument 4 must be a scalar or a vector with a length equal to the number of observation points!");
    return oct_retval;
  }

  const Matrix tmp3 = args(3).matrix_value();
  delay = (double*) tmp3.fortran_vec();

  //
  // Material parameters
  //

  // Check that arg 5 is a 2 element vector.
  if (!((mxGetM(4)==2 && mxGetN(4)==1) || (mxGetM(4)==1 && mxGetN(4)==2))) {
    error("Argument 5 must be a vector of length 2!");
    return oct_retval;
  }

  const Matrix tmp4 = args(4).matrix_value();
  m_par = (double*) tmp4.fortran_vec();
  v     = m_par[0]; // Normal velocity of transducer surface.
  cp    = m_par[1]; // Sound speed.

  //
  // Compute device
  //

  if (nrhs == 6) {

    if (!mxIsChar(5)) {
      error("Argument 6 must be a string");
      return oct_retval;
    }

    device = args(5).string_value();
    octave_stdout << "Compute device: " <<   device  << "\n";
  }

  //
  // Create an output matrix for the impulse response(s).
  //

  Matrix h_mat(nt, no);
  h = h_mat.fortran_vec();

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
  // Call the rect_sir subroutine.
  //

  running = true;

  // Check if we should use the GPU

#ifdef USE_OPENCL
  if (device == "gpu") {
    cl_rect_sir(ro, no, a, b, dt, nt,  delay[0], v, cp, h);

  } else { // Otherwise use the cpu
#endif

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
    if (nthreads > (unsigned int) no) {
      nthreads = no;
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
    D[thread_n].a = a;
    D[thread_n].b = b;
    D[thread_n].dt = dt;
    D[thread_n].nt = nt;

    if (mxGetM(3) * mxGetN(3) == 1)
      D[thread_n].delay_type = DelayType::single; // delay is a scalar.
    else
      D[thread_n].delay_type = DelayType::multiple; // delay is a vector.

    D[thread_n].delay = delay;
    D[thread_n].v = v;
    D[thread_n].cp = cp;
    //D[thread_n].alpha = alpha;
    D[thread_n].h = h;
    //D[thread_n].err_level = err_level;

    // Starts the threads.
    threads[thread_n] = std::thread(smp_dream_rect_sir, &D[thread_n]); // Start the threads.
    set_dream_thread_affinity(thread_n, nthreads, threads);
  }

  // Wait for all threads to finish.
  for (thread_n = 0; thread_n < nthreads; thread_n++)
    threads[thread_n].join();

  // Free memory.
  free((void*) D);

#ifdef USE_OPENCL
  }
#endif

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

  oct_retval.append(h_mat);

  return oct_retval;
}
