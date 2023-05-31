/***
*
* Copyright (C) 2008,2009,2015,2016,2021,2023 Fredrik Lingvall
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

#include <iostream>
#include <csignal>
#include <thread>
#include <mutex>

#include "affinity.h"
#include "circ_sir.h"
#include "arg_parser.h"

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
  double r;
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

void* smp_dream_circ_sir(void *arg)
{
  //ErrorLevel tmp_err=ErrorLevel::none, err=ErrorLevel::none;
  DATA D = *(DATA *)arg;
  double xo, yo, zo;
  double *h = D.h;
  double r=D.r, dt=D.dt;
  octave_idx_type n, no=D.no, nt=D.nt;
  //int    tmp_lev, err_level=D.err_level;
  double *delay=D.delay, *ro=D.ro, v=D.v, cp=D.cp;// alpha=D.alpha;
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

      circ_sir(xo,yo,zo,r,dt,nt,delay[0],v,cp,&h[n*nt]); // TODO: Add attenuation.

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

      circ_sir(xo,yo,zo,r,dt,nt,delay[n],v,cp,&h[n*nt]); // TODO: Add attenuation.

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
 * circ_sir - Octave (oct) gateway function for CIRC_SIR.
 *
 ***/

DEFUN_DLD (circ_sir, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {}  [Y] = circ_sir(Ro,geom_par,delay,s_par,m_par).\n\
\n\
CIRC_SIR Computes the time-continous (analytic) spatial impulse response(s) for a\n\
circular transducer.\n\
\n\
Observation point(s) ([mm]):\n\
\n\
@table @code\n\
@item Ro\n\
An N x 3 matrix, Ro = [xo1 yo1 zo2; xo2 yo2 zo2; ... xoN yoN zoN]; where N is the number of observation points.\n\
@end table\n\
\n\
Geometrical parameters: geom_par = [r];\n\
\n\
@table @code\n\
@item r\n\
Radius of the transducer [mm].\n\
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
\n\
@end table\n\
\n\
circ_sir is an oct-function that is a part of the DREAM Toolbox available at\n\
@url{https://github.com/frli8848/DREAM}.\n\
\n\
Copyright @copyright{} 2008-2023 Fredrik Lingvall.\n\
@seealso {dreamcirc, scirc_sir}\n\
@end deftypefn")
{
  double *ro,*geom_par,*s_par,*m_par;
  dream_idx_type nt, no;
  double r, dt;
  double *delay,v,cp;
  double *h;
  DATA   *D;
  octave_idx_type start, stop;
  std::thread *threads;
  unsigned int thread_n, nthreads;
  sighandler_t  old_handler, old_handler_abrt, old_handler_keyint;
  octave_value_list oct_retval;

  int nrhs = args.length ();

  ArgParser ap;

  // Check for proper number of arguments

  if (!ap.check_arg_in("circ_sir", nrhs, 5, 5)) {
    return oct_retval;
  }

  if (!ap.check_arg_out("circ_sir", nlhs, 0, 1)) {
    return oct_retval;
  }

  //
  // Observation point.
  //

  if (!ap.check_obs_points("circ_sir", args, 0)) {
    return oct_retval;
  }

  no = mxGetM(0); // Number of observation points.
  const Matrix tmp0 = args(0).matrix_value();
  ro = (double*) tmp0.data();

  //
  // Transducer geometry
  //

  if (!ap.check_geometry("circ_sir", args, 1, 1)) {
    return oct_retval;
  }

  const Matrix tmp1 = args(1).matrix_value();
  geom_par = (double*) tmp1.data();
  r = geom_par[0];		// Radius of the transducer.

  //
  // Temporal and spatial sampling parameters.
  //

  if (!ap.check_sampling("circ_sir", args, 2, 2)) {
    return oct_retval;
  }

  const Matrix tmp2 = args(2).matrix_value();
  s_par = (double*) tmp2.data();
  dt    = s_par[0]; // Temporal discretization size (= 1/sampling freq).
  nt    = (dream_idx_type) s_par[1];	// Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //

  if (!ap.check_delay("circ_sir", args, 3, no)) {
    return oct_retval;
  }

  const Matrix tmp3 = args(3).matrix_value();
  delay = (double*) tmp3.data();

  //
  // Material parameters
  //

  if (!ap.check_material("circ_sir", args, 4, 2)) {
    return oct_retval;
  }

  const Matrix tmp4 = args(4).matrix_value();
  m_par = (double*) tmp4.data();
  v     = m_par[0]; // Normal velocity of transducer surface.
  cp    = m_par[1]; // Sound speed.

  //
  // Number of threads
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

  //
  // Create an output matrix for the impulse response(s).
  //

  Matrix h_mat(nt, no);
  h = (double*) h_mat.data();

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
  // Call the analytic SIR subroutine.
  //

  //out_err = ErrorLevel::none;
  running = true;

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

    // Start the threads.
    threads[thread_n] = std::thread(smp_dream_circ_sir, &D[thread_n]);
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

  oct_retval.append(h_mat);

  return oct_retval;
}
