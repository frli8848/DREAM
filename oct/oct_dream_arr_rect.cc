/***
*
* Copyright (C) 2006,2007,2008,2009,2012,2014,2015,2016,2019,2021 Fredrik Lingvall
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

#include <thread>
#include <mutex>

#include "dream_arr_rect.h"
#include "affinity.h"
#include "dream_error.h"
#include "arr_functions.h"

#define SINGLE 0
#define MULTIPLE 1

//
// Octave headers.
//

#include <octave/oct.h>

//
// Macros
//

#define mxGetM(N)   args(N).matrix_value().rows()
#define mxGetN(N)   args(N).matrix_value().cols()
#define mxIsChar(N) args(N).is_string()

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
  octave_idx_type no;
  octave_idx_type start;
  octave_idx_type stop;
  double *ro;
  double a;
  double b;
  double dx;
  double dy;
  double dt;
  octave_idx_type nt;
  int delay_method;
  double *delay;
  double v;
  double cp;
  Attenuation *att;
  int num_elements;
  double *G;
  int focus_met;
  int steer_met;
  bool do_apod;
  int apod_met;
  double *focal;
  double *apod;
  double theta;
  double phi;
  double param;
  double *h;
  int err_level;
} DATA;


typedef void (*sighandler_t)(int);

//
// Function prototypes.
//
void* smp_dream_arr_rect(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

/***
 *
 * Thread function.
 *
 ***/

void* smp_dream_arr_rect(void *arg)
{
  int tmp_err = NONE, err = NONE, n;
  DATA D = *(DATA *)arg;
  double xo, yo, zo;
  double *h = D.h;
  double a=D.a, b=D.b, dx=D.dx, dy=D.dy, dt=D.dt;
  octave_idx_type no=D.no, nt=D.nt;
  int    tmp_lev, err_level=D.err_level;
  double *delay=D.delay, *ro=D.ro, v=D.v, cp=D.cp;
  Attenuation *att = D.att;
  octave_idx_type start=D.start, stop=D.stop;
  int    focus_met=D.focus_met, steer_met=D.steer_met, do_apod = D.do_apod,apod_met=D.apod_met;
  double *focal=D.focal, *apod=D.apod, theta=D.theta,phi=D.phi,param=D.param;
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
  if (err_level == STOP)
    tmp_lev = PARALLEL_STOP;
  else
    tmp_lev = err_level;

  for (n=start; n<stop; n++) {
    xo = ro[n];
    yo = ro[n+1*no];
    zo = ro[n+2*no];

    double dlay = 0.0;
    if (D.delay_method == SINGLE) {
      dlay = delay[0];
    } else { // MULTIPLE delays.
      dlay = delay[n];
    }

    if (att == nullptr) {
      err = dream_arr_rect(xo, yo, zo,
                           a, b,
                           dx, dy, dt, nt,
                           dlay, v, cp,
                           num_elements, gx, gy, gz,
                           focus_met, focal, steer_met, theta, phi, apod, do_apod, apod_met, param,
                           &h[n*nt],tmp_lev);

    } else {
      err = dream_arr_rect(*att, *xc_vec, *x_vec,
                           xo, yo, zo,
                           a, b,
                           dx, dy, dt, nt,
                           dlay, v, cp,
                           num_elements, gx, gy, gz,
                           focus_met, focal, steer_met, theta, phi, apod, do_apod, apod_met, param,
                           &h[n*nt],tmp_lev);
    }

    if (err != NONE || out_err ==  PARALLEL_STOP) {
      tmp_err = err;
      if (err == PARALLEL_STOP || out_err ==  PARALLEL_STOP)
        break; // Jump out when a STOP error occurs.
    }

    if (!running) {
      octave_stdout << "Thread for observation points " << start+1 << " -> " << stop << " bailing out!\n";
      return(NULL);
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
 *  Octave (oct) gateway function for (parallel) dreamrect.
 *
 ***/
DEFUN_DLD (dream_arr_rect, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} [H,err] = [H,err] = dream_arr_rect(Ro,geom_par,G,s_par,delay,m_par,foc_met,...\n\
                focal,steer_met,steer_par,apod_met,apod,win_par,err_level);\n\
\n\
DREAM_ARR_RECT - Computes the spatial impulse response\n\
for a rectangular transducer using parallel processing\n\
(using threads).\n\
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
Element width in x-direction [mm].\n\
@item b\n\
Element width in y-direction [mm].\n\
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
odization parameters: apod_met, apod, and win_par. The apod_met (apodization method) options are:\n\
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
dream_arr_rect is an oct-function that is a part of the DREAM Toolbox available at\n\
@url{http://www.signal.uu.se/Toolbox/dream/}.\n\
\n\
Copyright @copyright{} 2006-2019 Fredrik Lingvall.\n\
@seealso {dreamrect}\n\
@end deftypefn")
{
  double *ro,*geom_par,*s_par,*m_par;
  double *steer_par;
  char   apod_str[50], foc_str[50], steer_str[50];
  int    buflen;
  double a,b,dx,dy,dt;
  octave_idx_type nt,no,n;
  double param=0,*delay,v,cp,alpha;
  int    num_elements;
  double *G;
  int    focus_met=0;
  double *focal=nullptr;
  int    steer_met=0;
  double theta=0.0, phi=0.0, *apod=nullptr;
  bool   do_apod=false;
  int    apod_met=0;
  double *h, *err_p;
  int    err_level=STOP, is_set = false;
  char   err_str[50];
  DATA   *D;
  octave_idx_type start, stop;
  std::thread *threads;
  unsigned int thread_n, nthreads;
  sighandler_t old_handler, old_handler_abrt, old_handler_keyint;
  octave_value_list oct_retval;

  int nrhs = args.length ();

  // Check for proper number of arguments

  if (!((nrhs == 13) || (nrhs == 14))) {
    error("dream_arr_rect requires 13 or 14 input arguments!");
    return oct_retval;
  }
  else
    if (nlhs > 2) {
      error("Too many output arguments for dream_arr_rect!");
      return oct_retval;
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

  // Check that arg 2 is a 2 element vector
  if (!((mxGetM(1)==2 && mxGetN(1)==1) || (mxGetM(1)==1 && mxGetN(1)==2))) {
    error("Argument 2 must be a vector of length 2!");
    return oct_retval;
  }
  const Matrix tmp1 = args(1).matrix_value();
  geom_par = (double*) tmp1.fortran_vec();
  a = geom_par[0];		// x-width of the transducer.
  b  = geom_par[1];		// y-width of the transducer.

  //
  // Grid function (position vectors of the elements).
  //

  num_elements = (int) mxGetM(2); // Number of elementents in the array.
  if (mxGetN(2) !=3 ) {
    error("Argument 3  must a (number of array elements) x 3 matrix!");
    return oct_retval;
  }
  const Matrix tmp2 = args(2).matrix_value();
  G = (double*) tmp2.fortran_vec(); // First column in the matrix.
  //gy    = gx + num_elements;		// Second column in the matrix.
  //gz    = gy + num_elements;		// Third column in the matrix.

  //
  // Temporal and spatial sampling parameters.
  //

  // Check that arg 4 is a 4 element vector
  if (!((mxGetM(3)==4 && mxGetN(3)==1) || (mxGetM(3)==1 && mxGetN(3)==4))) {
    error("Argument 4 must be a vector of length 4!");
    return oct_retval;
  }
  const Matrix tmp3 = args(3).matrix_value();
  s_par = (double*) tmp3.fortran_vec();
  dx    = s_par[0];		// Spatial x discretization size.
  dy    = s_par[1];		// Spatial dy iscretization size.
  dt    = s_par[2];		// Temporal discretization size (= 1/sampling freq).
  nt    = (octave_idx_type) s_par[3];	// Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //

  // Check that arg 5 is a scalar or a vector.
  if ( (mxGetM(4) * mxGetN(4) !=1) && ((mxGetM(4) * mxGetN(4)) != no)) {
    error("Argument 5 must be a scalar or a vector with a length equal to the number of observation points!");
    return oct_retval;
  }
  const Matrix tmp4 = args(4).matrix_value();
  delay = (double*) tmp4.fortran_vec();

  //
  // Material parameters
  //

 // Check that arg 6 is a 3 element vectora
  if (!((mxGetM(5)==3 && mxGetN(5)==1) || (mxGetM(5)==1 && mxGetN(5)==3))) {
    error("Argument 6 must be a vector of length 3!");
    return oct_retval;
  }
  const Matrix tmp5 = args(5).matrix_value();
  m_par = (double*) tmp5.fortran_vec();
  v     = m_par[0]; // Normal velocity of transducer surface.
  cp    = m_par[1]; // Sound speed.
  alpha  = m_par[2]; // Attenuation coefficient [dB/(cm MHz)].

  //
  // Focusing parameters.
  //

  if (nrhs >= 7) {

    if (!mxIsChar(6)) {
      error("Argument 7 must be a string");
      return oct_retval;
    }
    std::string strin = args(6).string_value();
    buflen = strin.length();
    for (n=0; n<=buflen; n++ ) {
      foc_str[n] = strin[n];
    }
    foc_str[buflen] = '\0';

    is_set = false;

    if (!strcmp(foc_str,"off")) {
      focus_met =  NO_FOCUS;
      is_set = true;
    }

    if (!strcmp(foc_str,"x")) {
      focus_met = FOCUS_X;
      is_set = true;
    }

    if (!strcmp(foc_str,"y")) {
      focus_met = FOCUS_Y;
      is_set = true;
    }

    if (!strcmp(foc_str,"xy")) {
      focus_met = FOCUS_XY;
      is_set = true;
    }

    if (!strcmp(foc_str,"x+y")) {
      focus_met = FOCUS_X_Y;
      is_set = true;
    }

    if (!strcmp(foc_str,"ud")) {
      focus_met = FOCUS_UD;
      is_set = true;

      if (mxGetM(7) * mxGetN(7) != num_elements ) {
        error("The time delay vector (argument 8) for user defined ('ud') focusing\n") ;
        error("delays must have the same length as the number of array elements.!");
        return oct_retval;
      }
      const Matrix tmp4 = args(7).matrix_value();
      focal = (double*) tmp4.fortran_vec();
    }
    else {

      // Check that arg 8 is a scalar.
      if (mxGetM(7) * mxGetN(7) !=1 ) {
        error("Argument 8  must be a scalar!");
        return oct_retval;
      }
      // Focal point (in mm).
      const Matrix tmp4 = args(7).matrix_value();
      focal = (double*) tmp4.fortran_vec();
    }

    if (is_set == false) {
      error("Unknown focusing method!");
      return oct_retval;
    }

  } else {
    focus_met = NO_FOCUS;
  }

  //
  // Beam steering.
  //

  if (nrhs >= 9) {

    if (!mxIsChar(8)) {
      error("Argument 9 must be a string");
      return oct_retval;
    }
    std::string strin = args(8).string_value();
    buflen = strin.length();
    for ( n=0; n<=buflen; n++ ) {
      steer_str[n] = strin[n];
    }
    steer_str[buflen] = '\0';

    steer_met = NO_STEER;      // Default no steering
    is_set = false;

    if (!strcmp(steer_str,"off")) {
      steer_met = NO_STEER;
      is_set = true;
    }

    if (!strcmp(steer_str,"x")) {
      steer_met = STEER_X;
      is_set = true;
    }

    if (!strcmp(steer_str,"y")) {
      steer_met = STEER_Y;
      is_set = true;
    }

    if (!strcmp(steer_str,"xy")) {
      steer_met = STEER_XY;
      is_set = true;
    }

    if (is_set == false) {
      error("Unknown beamsteering method!");
      return oct_retval;
    }

    // Check that arg 10 is a 2 element vector
    if (!((mxGetM(9)==2 && mxGetN(9)==1) || (mxGetM(9)==1 && mxGetN(9)==2))) {
      dream_err_msg("Argument 10 must be a vector of length 2!");
      return oct_retval;
    }
    const Matrix tmp5 = args(9).matrix_value();
    steer_par = (double*) tmp5.fortran_vec();
    theta  = steer_par[0];		// Angle in x-direction.
    phi    = steer_par[1];		// Angle in y-direction.

  } else {
    steer_met = NO_STEER;
  }

  //
  // Apodization.
  //

  if (nrhs >= 11) {

    if (!mxIsChar(10)) {
      error("Argument 11 must be a string");
      return oct_retval;
    }
    std::string strin = args(10).string_value();
    buflen = strin.length();
    for ( n=0; n<=buflen; n++ ) {
      apod_str[n] = strin[n];
    }
    apod_str[buflen] = '\0';

    do_apod = false;			// default off.
    is_set = false;

    if (!strcmp(apod_str,"off")) {
      do_apod = false;
      is_set = true;
    }

    if (!strcmp(apod_str,"ud")) {
      do_apod = true;
      apod_met = APOD_UD;
      is_set = true;

      // Vector of apodization weights.
      if (mxGetM(11) * mxGetN(11) != num_elements) {
        error("The length of argument 12 (apodization vector) must be the same as the number of array elements!");
        return oct_retval;
      }
      const Matrix tmp6 = args(11).matrix_value();
      apod = (double*) tmp6.fortran_vec();
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
      error("Unknown apodization!");
      return oct_retval;
    }

    // Parameter for raised cos and Gaussian apodization functions.
    if (mxGetM(12) * mxGetN(12) !=1 ) {
      error("Argument 13 must be a scalar");
      return oct_retval;
    }

    const Matrix tmp7 = args(12).matrix_value();
    param = (double) tmp7.fortran_vec()[0];

  } else {
    do_apod = false;
  }

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11).
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

    if (!mxIsChar(13)) {
      error("Argument 14 must be a string");
      return oct_retval;
    }

    std::string strin = args(13).string_value();
    buflen = strin.length();
    for (int n=0; n<=buflen; n++ ) {
      err_str[n] = strin[n];
    }
    err_str[buflen] = '\0';


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
      error("Unknown error level!");
      return oct_retval;
    }
  }
  else
    err_level = STOP; // Default.


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

  out_err = NONE;
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
    D[thread_n].dx = dx;
    D[thread_n].dy = dy;
    D[thread_n].dt = dt;
    D[thread_n].nt = nt;

    if (mxGetM(4) * mxGetN(4) == 1)
      D[thread_n].delay_method = SINGLE; // delay is a scalar.
    else
      D[thread_n].delay_method = MULTIPLE; // delay is a vector.

    D[thread_n].delay = delay;
    D[thread_n].v = v;
    D[thread_n].cp = cp;
    D[thread_n].att = att_ptr;
    D[thread_n].num_elements = num_elements;
    D[thread_n].G = G;
    D[thread_n].focus_met = focus_met;
    D[thread_n].steer_met = steer_met;
    D[thread_n].do_apod = do_apod;
    D[thread_n].apod_met = apod_met;
    D[thread_n].focal = focal;
    D[thread_n].apod = apod;
    D[thread_n].theta = theta;
    D[thread_n].phi = phi;
    D[thread_n].param = param;
    D[thread_n].h = h;
    D[thread_n].err_level = err_level;

    // Start the threads.
    threads[thread_n] = std::thread(smp_dream_arr_rect, &D[thread_n]); // Start the threads.
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

  if ( (err_level == STOP) && (out_err != NONE)) {
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
