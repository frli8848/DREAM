/***
*
* Copyright (C) 2006,2007,2008,2009,2012,2014,2015,2016,2021,2021 Fredrik Lingvall
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

#include "dream_arr_cylind.h"
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
  dream_idx_type no;
  dream_idx_type start;
  dream_idx_type stop;
  double *ro;
  double a;
  double b;
  double Rcurv;
  double dx;
  double dy;
  double dt;
  dream_idx_type nt;
  DelayType delay_type;
  double *delay;
  double v;
  double cp;
  Attenuation *att;
  dream_idx_type num_elements;
  double *G;
  FocusMet foc_met;
  SteerMet steer_met;
  bool do_apod;
  ApodMet apod_met;
  double *focal;
  double *apod;
  double theta;
  double phi;
  double apod_par;
  double *h;
  ErrorLevel err_level;
} DATA;


typedef void (*sighandler_t)(int);

//
// Function prototypes.
//
void* smp_dream_arr_cylind(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

/***
 *
 * Thread function.
 *
 ***/

void* smp_dream_arr_cylind(void *arg)
{
  ErrorLevel tmp_err = ErrorLevel::none, err = ErrorLevel::none;
  DATA D = *(DATA *)arg;
  double xo, yo, zo;
  double *h = D.h;
  double a=D.a, b=D.b, Rcurv=D.Rcurv, dx=D.dx, dy=D.dy, dt=D.dt;
  dream_idx_type no=D.no, nt=D.nt;
  ErrorLevel tmp_lev, err_level=D.err_level;
  double *delay=D.delay, *ro=D.ro, v=D.v, cp=D.cp;
  Attenuation *att = D.att;
  dream_idx_type start=D.start, stop=D.stop;
  FocusMet foc_met=D.foc_met;
  SteerMet steer_met=D.steer_met;
  int do_apod = D.do_apod;
  ApodMet apod_met = D.apod_met;
  double *focal=D.focal, *apod=D.apod, theta=D.theta,phi=D.phi,apod_par=D.apod_par;
  int    num_elements = D.num_elements;

  double *gx = D.G;               // First column in the matrix.
  double *gy = gx + num_elements; // Second column in the matrix.
  double *gz = gy + num_elements; // Third column in the matrix.

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

  for (dream_idx_type n=start; n<stop; n++) {
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
      err = dream_arr_cylind(xo, yo, zo,
                             a, b, Rcurv,
                             dx, dy, dt, nt,
                             dlay, v, cp,
                             num_elements, gx, gy, gz,
                             foc_met, focal,
                             steer_met, theta, phi,
                             apod, do_apod, apod_met, apod_par,
                             &h[n*nt], tmp_lev);
    } else {
      err = dream_arr_cylind(*att, *xc_vec, *x_vec,
                             xo, yo, zo,
                             a, b, Rcurv,
                             dx, dy, dt, nt,
                             dlay, v, cp,
                             num_elements, gx, gy, gz,
                             foc_met, focal,
                             steer_met, theta, phi,
                             apod, do_apod, apod_met, apod_par,
                             &h[n*nt], tmp_lev);

    }

    if (err != ErrorLevel::none || out_err ==  ErrorLevel::parallel_stop) {
      tmp_err = err;
      if (err == ErrorLevel::parallel_stop || out_err ==  ErrorLevel::parallel_stop)
        break; // Jump out when a ErrorLevel::stop error occurs.
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
 *  Octave (OCT) gateway function for parallel dream_arr_cylind.
 *
 ***/

DEFUN_DLD (dream_arr_cylind, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} [H,err] = dream_arr_cylind(ro,geom_par,G,s_par,delay,m_par,foc_met,focal,...\n\
                     steer_met,steer_par,apod_met,apod,win_par,err_level);\n\
\n\
DREAM_ARR_CYLIND - Computes the spatial impulse response\n\
for an array with cylindrical concave (focused) elements.\n\
using parallel processing (using threads).\n\
\n\
Observation point(s) ([mm]):\n\
\n\
@table @code\n\
@item Ro\n\
An N x 3 matrix, Ro = [xo1 yo1 zo2; xo2 yo2 zo2; ... xoN yoN zoN]; where N is the number of observation points.\n\
@end table\n\
\n\
Geometrical parameters: geom_par = [a b Rcurv];\n\
\n\
@table @code\n\
@item a\n\
Element width in x-direction [mm].\n\
@item b\n\
Element width in y-direction [mm].\n\
@item Rcurv\n\
Radius of the curvature [mm]. If Rcurv > 0 then it is focsued/concave and if Rcurv < 0 it is convex/defocused.\n\
@end table\n\
\n\
Array grid parameter:\n\
\n\
@table @code\n\
@item G\n\
An L x 3 grid function matrix. Column 1 contain the x-postions \n\
of the elements, column 2 the y-positions, and column 3\n\
the z-positions, and where L is the number of elements.\n\
@end table\n\
\n\
\n\
@table @code\n\
@item dx\n\
Spatial x-direction discretization size [mm] .\n\
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
Attenuation coefficient [dB/(cm MHz)] .\n\
\n\
@end table\n\
\n\
Focusing parameters: foc_met and focal:\n\
\n\
@table @code\n\
@item foc_met\n\
Focusing method, options are: 'off', 'x', 'y', 'xy', 'x+y', and 'ud'.\n\
@item  focal\n\
Focal distance [mm]. If foc_met = 'ud' (user defined) then focal is a vector of focusing delays.\n\
@end table\n\
\n\
Beam steering parameters: steer_met and steer_par:\n\
\n\
@table @code\n\
@item steer_met\n\
Beam steering method, options are: 'off', 'x', 'y', and 'xy'.\n\
@item  steer_par =  [theta phi];\n\
theta [deg] is the x-direction steer angle and \n\
phi [deg] is the y-direction steer angle.\n\
@end table\n\
\n\
@table @code\n\
\n\
@item 'off'\n\
No apodization.\n\
@item 'ud'\n\
User defined apodization.\n\
@item  'triangle'\n\
Triangle window.\n\
@item 'gauss'\n\
Gaussian (bell-shaped) window.\n\
@item 'raised'\n\
Raised cosine.\n\
@item 'simply'\n\
Simply supported.\n\
@item 'clamped'\n\
Clamped.\n\
@end table\n\
\n\
and the apod and win_par parameters are:\n\
\n\
@table @code\n\
@item apod\n\
Vector of apodiztion weights (used for the 'ud' option).\n\
@item win_par\n\
A scalar parameter for raised cosine and Gaussian apodization functions.\n\
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
dream_arr_cylind is an oct-function that is a part of the DREAM Toolbox available at\n\
@url{https://github.com/frli8848/DREAM}.\n\
\n\
Copyright @copyright{} 2006-2021 Fredrik Lingvall.\n\
@seealso {dream_arr_cylind_d, dreamcylind}\n\
@end deftypefn")
{
  double *ro;
  double apod_par=0.0, *delay;
  double *G;
  FocusMet foc_met=FocusMet::none;
  SteerMet steer_met=SteerMet::none;
  double theta=0.0, phi=0.0;
  bool   do_apod=false;
  ApodMet apod_met=ApodMet::gauss;
  double *h;
  ErrorLevel err_level=ErrorLevel::stop;
  DATA *D;
  dream_idx_type start, stop;
  std::thread *threads;
  unsigned int thread_n, nthreads;
  sighandler_t old_handler, old_handler_abrt, old_handler_keyint;
  octave_value_list oct_retval;

  int nrhs = args.length ();

  ArgParser ap;

  // Check for proper number of arguments

  if (!ap.check_arg_in("dream_arr_cylind", nrhs, 13, 14)) {
    return oct_retval;
  }

  if (!ap.check_arg_out("dream_arr_cylind", nlhs, 0, 2)) {
    return oct_retval;
  }

  //
  // Observation points.
  //

  if (!ap.check_obs_points("dream_arr_cylind", args, 0)) {
    return oct_retval;
  }

  dream_idx_type no = mxGetM(0); // Number of observation points.
  const Matrix tmp0 = args(0).matrix_value();
  ro = (double*) tmp0.fortran_vec();

  //
  // Transducer geometry
  //

  double a=0.0, b=0.0, Rcurv=0.0;
  if (!ap.parse_geometry("dream_arr_cylind", args, 1, 3, a, b, Rcurv)) {
    return oct_retval;
  }

  //
  // Grid function (position vectors of the elements).
  //

  if (!ap.check_array("dream_arr_cylind", args, 2)) {
    return oct_retval;
  }

  dream_idx_type num_elements = (dream_idx_type) mxGetM(2); // Number of elementents in the array.

  const Matrix tmp2 = args(2).matrix_value();
  G = (double*) tmp2.fortran_vec(); // First column in the matrix.
  //gy    = gx + num_elements;		// Second column in the matrix.
  //gz    = gy + num_elements;		// Third column in the matrix.

  //
  // Temporal and spatial sampling parameters.
  //

  double dx=0.0, dy=0.0, dt=0.0;
  dream_idx_type nt=0;
  if (!ap.parse_sampling("dream_arr_cylind", args, 3, 4, dx, dy, dt, nt)) {
    return oct_retval;
  }

  //
  // Start point of impulse response vector ([us]).
  //

  if (!ap.check_delay("dream_arr_cylind", args, 4, no)) {
    return oct_retval;
  }

  const Matrix tmp4 = args(4).matrix_value();
  delay = (double*) tmp4.fortran_vec();

  //
  // Material parameters
  //

  double v=1.0, cp=1000.0, alpha=0.0;
  if (!ap.parse_material("dream_arr_cylind", args, 5, v, cp, alpha)) {
    return oct_retval;
  }

  //
  // Focusing parameters.
  //

  // Allocate space for the user defined focusing delays
  std::unique_ptr<double[]> focal = std::make_unique<double[]>(num_elements);

  if (nrhs >= 7) {
    if (!ap.parse_focus_args("dream_arr_cylind", args, 6, foc_met, focal.get())) {
      return oct_retval;
    }
  } else {
    foc_met = FocusMet::none;
  }

  //
  // Beam steering.
  //

  if (nrhs >= 9) {
    if (!ap.parse_steer_args("dream_arr_cylind", args, 8, steer_met, theta, phi)) {
      return oct_retval;
    }
  } else {
    steer_met = SteerMet::none;
  }

  //
  // Apodization.
  //

  // Allocate space for the user defined apodization weights
  std::unique_ptr<double[]> apod = std::make_unique<double[]>(num_elements);

  if (nrhs >= 11) {

    if (!ap.parse_apod_args("dream_arr_cylind", args, 10, num_elements,
                            do_apod, apod.get(), apod_met, apod_par)) {
      return oct_retval;
    }
  } else {
    do_apod = false;
  }

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hyper threading, C++11)
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

  if (nrhs == 14) {
    if (!ap.parse_error_arg("dream_arr_cylind", args, 13, err_level)) {
      return oct_retval;
    }
  } else {
    err_level = ErrorLevel::stop; // Default.
  }

  // Create an output matrix for the impulse response
  Matrix h_mat(nt, no);
  h = h_mat.fortran_vec();

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
    D[thread_n].a = a;
    D[thread_n].b = b;
    D[thread_n].Rcurv = Rcurv;
    D[thread_n].dx = dx;
    D[thread_n].dy = dy;
    D[thread_n].dt = dt;
    D[thread_n].nt = nt;

    if (mxGetM(4) * mxGetN(4) == 1)
      D[thread_n].delay_type = DelayType::single; // delay is a scalar.
    else
      D[thread_n].delay_type = DelayType::multiple; // delay is a vector.

    D[thread_n].delay = delay;
    D[thread_n].v = v;
    D[thread_n].cp = cp;
    D[thread_n].att = att_ptr;
    D[thread_n].num_elements = num_elements;
    D[thread_n].G = G;
    D[thread_n].foc_met = foc_met;
    D[thread_n].steer_met = steer_met;
    D[thread_n].do_apod = do_apod;
    D[thread_n].apod_met = apod_met;
    D[thread_n].focal = focal.get();
    D[thread_n].apod = apod.get();
    D[thread_n].theta = theta;
    D[thread_n].phi = phi;
    D[thread_n].apod_par = apod_par;
    D[thread_n].h = h;
    D[thread_n].err_level = err_level;

    if (nthreads > 1) {
      // Start the threads.
      threads[thread_n] = std::thread(smp_dream_arr_cylind, &D[thread_n]); // Start the threads.
      set_dream_thread_affinity(thread_n, nthreads, threads);
    } else {
      smp_dream_arr_cylind(&D[0]);
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

  if (std::signal(SIGTERM, old_handler) == SIG_ERR) {
    std::cerr << "Couldn't register old SIGTERM signal handler!" << std::endl;
  }

  if (std::signal(SIGABRT,  old_handler_abrt) == SIG_ERR) {
    std::cerr << "Couldn't register old SIGABRT signal handler!" << std::endl;
  }

  if (std::signal(SIGINT, old_handler_keyint) == SIG_ERR) {
    std::cerr << "Couldn't register old SIGINT signal handler!" << std::endl;
  }

  if (!running) {
    error("CTRL-C pressed!\n"); // Bail out.
    return oct_retval;
  }

  //
  // Check for Error.
  //

  if ( (err_level == ErrorLevel::stop) && (out_err != ErrorLevel::none)) {
    error("Error in dream_arr_cylind!"); // Bail out if error.
    return oct_retval;
  }

  oct_retval.append(h_mat);

  //
  // Return error.
  //

  if (nlhs == 2) {
    Matrix err_mat(1, 1);
    double *err_p = err_mat.fortran_vec();
    err_p[0] = (double) out_err;
    oct_retval.append(err_mat);
  }

  return oct_retval;
}
