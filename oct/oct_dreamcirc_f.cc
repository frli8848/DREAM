/***
 *
 * Copyright (C) 2006,2007,2008,2009,2012,2014,2016,2019,2021 Fredrik Lingvall
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

#include <signal.h>

#include <thread>
#include <mutex>

#include "dreamcirc_f.h"
#include "affinity.h"
#include "arg_parser.h"

//
// Octave headers.
//

#include <octave/oct.h>

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
  octave_idx_type no;
  octave_idx_type start;
  octave_idx_type stop;
  double *ro;
  double R;
  FocusMet foc_met;
  double focal;
  double dx;
  double dy;
  double dt;
  octave_idx_type nt;
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
  octave_idx_type n, no=D.no, nt=D.nt;
  ErrorLevel tmp_lev=ErrorLevel::none, err_level=D.err_level;
  double *delay=D.delay, *ro=D.ro, v=D.v, cp=D.cp, focal=D.focal;
  Attenuation *att = D.att;
  octave_idx_type start=D.start, stop=D.stop;
  FocusMet foc_met = D.foc_met;

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
                        R,
                        foc_met, focal,
                        dx, dy, dt,
                        nt, dlay, v, cp,
                        &h[n*nt], tmp_lev);
    } else {
      err = dreamcirc_f(*att, *xc_vec.get(),*x_vec.get(),
                        xo, yo, zo,
                        R,
                        foc_met, focal,
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
      octave_stdout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!\n";
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
 * Octave (oct) gateway function for parallel dreamcirc_f.
 *
 ***/
DEFUN_DLD (dreamcirc_f, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} [H,err] = dreamcirc_f(Ro,geom_par,s_par,delay,m_par,foc_met,focal,err_level)\n\
\n\
DREAMCIRC_F - Computes the spatial impulse response\n\
for a circular focused transducer using parallel processing\n\
(using threads).\n\
\n\
Observation point(s) ([mm]):\n\
\n\
@table @code\n\
@item Ro\n\
An N x 3 matrix, Ro = [xo1 yo1 zo2; xo2 yo2 zo2; ... xoN yoN zoN]; where N is the number of observation points.\n\
@end table\n\
\n\
Geometrical parameters: geom_par = [R];\n\
\n\
@table @code\n\
@item R\n\
Radius of the transducer elements [mm].\n\
@end table\n\
\n\
Sampling parameters: s_par = [dx dy dt nt]; \n\
\n\
@table @code\n\
@item dx\n\
Spatial x-direction discretization size [mm].\n\
@item dy\n\
Spatial y-direction discretization size [mm].\n\
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
Material parameters: m_par = [v cp alpha];\n\
\n\
@table @code\n\
@item v\n\
Normal velocity [m/s].\n\
@item cp\n\
Sound velocity [m/s].\n\
@item alpha\n\
Attenuation coefficient [dB/(cm MHz)].\n\
\n\
@end table\n\
\n\
Focusing parameters: foc_met and focal:\n\
\n\
@table @code\n\
@item foc_met\n\
Focusing method, options are: 'off', 'x', 'y', 'xy', and 'x+y'.\n\
@item  focal\n\
Focal distance [mm].\n\
@end table\n\
\n\
Error Handling: err_level;\n\
err_level is an optional text string parameter for controlling the error behavior, options are:\n\
\n\
@table @code\n\
@item 'ignore'\n\
An error is ignored (no error message is printed and the program is not stopped) but the err output \n\
argument is negative if an error occured.\n\
@item 'warn'\n\
An error message is printed but the program in not stopped (and err is negative).\n\
@item 'stop'\n\
An error message is printed and the program is stopped.\n\
@end table\n\
\n\
dreamcirc_f is an oct-function that is a part of the DREAM Toolbox available at\n\
@url{http://www.signal.uu.se/Toolbox/dream/}.\n\
\n\
Copyright @copyright{} 2006-2019 Fredrik Lingvall.\n\
@seealso {dreamcirc}\n\
@end deftypefn")
{
  double *ro,*geom_par,*s_par,*m_par;
  octave_idx_type nt, no;
  FocusMet foc_met=FocusMet::none;
  double R, dx, dy, dt;
  double *delay,v,cp,alpha,focal=0;
  double *h, *err_p;
  ErrorLevel err_level=ErrorLevel::stop;
  DATA *D;
  octave_idx_type start, stop;
  std::thread *threads;
  unsigned int thread_n, nthreads;
  sighandler_t   old_handler, old_handler_abrt, old_handler_keyint;
  octave_value_list oct_retval;

  int nrhs = args.length ();

  ArgParser ap;

  // Check for proper number of arguments

  if (!ap.check_arg_in("dreamcirc_f", nrhs, 7, 8)) {
    return oct_retval;
  }

  if (!ap.check_arg_out("dreamcirc_f", nlhs, 0, 2)) {
    return oct_retval;
  }

  //
  // Observation point.
  //

  if (!ap.check_obs_points("dreamcirc_f", args, 0)) {
    return oct_retval;
  }

  no = mxGetM(0); // Number of observation points.
  const Matrix tmp0 = args(0).matrix_value();
  ro = (double*) tmp0.fortran_vec();

  //
  // Transducer geometry
  //

  if (!ap.check_geometry("dreamcirc_f", args, 1, 1)) {
    return oct_retval;
  }

  const Matrix tmp1 = args(1).matrix_value();
  geom_par = (double*) tmp1.fortran_vec();
  R = geom_par[0];              // Radius of the transducer.

  //
  // Temporal and spatial sampling parameters.
  //

  if (!ap.check_sampling("dreamcirc_f", args, 2, 4)) {
    return oct_retval;
  }

  const Matrix tmp2 = args(2).matrix_value();
  s_par = (double*) tmp2.fortran_vec();
  dx = s_par[0];  // Spatial x-direction discretization size.
  dy = s_par[1];  // Spatial y-direction discretization size.
  dt = s_par[2];  // Temporal discretization size (= 1/sampling freq).
  nt = (octave_idx_type) s_par[3]; // Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //

  if (!ap.check_delay("dreamcirc_f", args, 3, no)) {
    return oct_retval;
  }

  const Matrix tmp3 = args(3).matrix_value();
  delay = (double*) tmp3.fortran_vec();

  //
  // Material parameters
  //

  if (!ap.check_material("dreamcirc_f", args, 4, 3)) {
    return oct_retval;
  }

  const Matrix tmp4 = args(4).matrix_value();
  m_par = (double*) tmp4.fortran_vec();
  v     = m_par[0]; // Normal velocity of transducer surface.
  cp    = m_par[1]; // Sound speed.
  alpha  = m_par[2]; // Attenuation coefficient [dB/(cm MHz)].

  //
  // Focusing parameters.
  //

  if (nrhs >= 6) {
    if (!ap.parse_focus_arg("dreamcirc_f", args, 5, foc_met, &focal)) {
      return oct_retval;
    }
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
  if (nthreads > (unsigned int) no) {
    nthreads = no;
  }

  //
  // Error reporting.
  //

  if (nrhs == 8) {
    if (!ap.parse_error_arg("dreamcirc_f", args, 7, err_level)) {
      return oct_retval;
    }
  } else {
    err_level = ErrorLevel::stop; // Default.
  }

  // Create an output matrix for the impulse response
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
  // Call the DREAM subroutine.
  //

  out_err = ErrorLevel::none;
  running = true;

  // Check if we have attenuation
  Attenuation att(nt, dt, alpha);
  Attenuation *att_ptr = nullptr;
  if (alpha > std::numeric_limits<double>::epsilon() ) {
    att_ptr = &att;
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

    if (mxGetM(3) * mxGetN(3) == 1)
      D[thread_n].delay_type = DelayType::single; // delay is a scalar.
    else
      D[thread_n].delay_type = DelayType::multiple; // delay is a vector.

    D[thread_n].delay = delay;
    D[thread_n].v = v;
    D[thread_n].cp = cp;
    D[thread_n].att = att_ptr;
    D[thread_n].h = h;
    D[thread_n].err_level = err_level;

    if (nthreads>1) {
      // Start the threads.
      threads[thread_n] = std::thread(smp_dream_circ_f, &D[thread_n]);
      set_dream_thread_affinity(thread_n, nthreads, threads);
    } else {
      smp_dream_circ_f(&D[0]);
    }
  }

  // Wait for all threads to finish.
  if (nthreads>1) {
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
    error("CTRL-C pressed!\n"); // Bail out.
    return oct_retval;
  }

  //
  // Check for Error.
  //

  if ( (err_level == ErrorLevel::stop) && (out_err != ErrorLevel::none)) {
    error(""); // Bail out if error.
    return oct_retval;
  }

  oct_retval.append(h_mat);

  //
  // Return error.
  //

  if (nlhs == 2) {
    Matrix err_mat(1, 1);
    err_p = err_mat.fortran_vec();
    err_p[0] = (double) out_err;
    oct_retval.append(err_mat);
  }

  return oct_retval;
}
