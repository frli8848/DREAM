/***
*
* Copyright (C) 2006,2007,2008,2009,2012,2015,2016,2019 Fredrik Lingvall
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
#include <thread>
#include <mutex>
#include <signal.h>
#include "dream_arr_annu.h"
#include "affinity.h"
#include "dream_error.h"

#define SINGLE 0
#define MULTIPLE 1

#ifdef USE_FFTW
#include "att.h"
#endif

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
  double dx;
  double dy;
  double dt;
  octave_idx_type nt;
  int delay_method;
  double *delay;
  double v;
  double cp;
  double alpha;
  int isize;
  double *gr;
  int ifoc;
  int iweight;
  int iapo;
  double focal;
  double *apod;
  double param;
  double  *ud_focal;
  double *h;
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
  double *h = D.h;
  double dx=D.dx, dy=D.dy, dt=D.dt;
  octave_idx_type n, no=D.no, nt=D.nt;
  int    tmp_lev, err_level=D.err_level;
  double *delay=D.delay, *ro=D.ro, v=D.v, cp=D.cp, alpha=D.alpha;
  octave_idx_type start=D.start, stop=D.stop;
  int    ifoc=D.ifoc, iweight=D.iweight,iapo=D.iapo;
  double focal=D.focal, *apod=D.apod, param=D.param;
  int    isize = D.isize;
  double *gr=D.gr;
  double *ud_focal=D.ud_focal;

  // Let the thread finish and then catch the error.
  if (err_level == STOP)
    tmp_lev = PARALLEL_STOP;
  else
    tmp_lev = err_level;

  if (ifoc != 6) {

    if (D.delay_method == SINGLE) {
      for (n=start; n<stop; n++) {
        xo = ro[n];
        yo = ro[n+1*no];
        zo = ro[n+2*no];

        err = dream_arr_annu(xo,yo,zo,dx,dy,dt,nt,delay[0],v,cp,alpha,isize,gr,
                             ifoc,focal,apod,iweight,iapo,param,&h[n*nt],tmp_lev);

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
    } else { // MULTIPLE delays.
      for (n=start; n<stop; n++) {
        xo = ro[n];
        yo = ro[n+1*no];
        zo = ro[n+2*no];

        err = dream_arr_annu(xo,yo,zo,dx,dy,dt,nt,delay[n],v,cp,alpha,isize,
                             gr,ifoc,focal,apod,iweight,iapo,param,&h[n*nt],tmp_lev);

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
    }
  }
  else { // User defined focusing.

    if (D.delay_method == SINGLE) {
      for (n=start; n<stop; n++) {
        xo = ro[n];
        yo = ro[n+1*no];
        zo = ro[n+2*no];

        err = dream_arr_annu_ud(xo,yo,zo,dx,dy,dt,nt,delay[0],v,cp,alpha,isize,gr,
                             ifoc,ud_focal,apod,iweight,iapo,param,&h[n*nt],tmp_lev);

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
    } else { // MULTIPLE delays.
      for (n=start; n<stop; n++) {
        xo = ro[n];
        yo = ro[n+1*no];
        zo = ro[n+2*no];

        err = dream_arr_annu_ud(xo,yo,zo,dx,dy,dt,nt,delay[n],v,cp,alpha,isize,
                                gr,ifoc,ud_focal,apod,iweight,iapo,param,&h[n*nt],tmp_lev);

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
 *  oct_dream_arr_annug.cc - Octave (oct) gateway function for parallel dreamrect.
 *
 ***/
DEFUN_DLD (dream_arr_annu, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} [H,err] = dream_arr_annu(Ro,G,s_par,delay,m_par,foc_met,focal,...\n\
                                               apod_met,apod,win_par,err_level);\n\
\n\
 DREAM_ARR_ANNU - Computes the spatial impulse response\n\
for an annular array using parallel processing (using threads).\n\
\n\
Observation point(s) ([mm]):\n\
\n\
@table @code\n\
@item Ro\n\
An N x 3 matrix, Ro = [xo1 yo1 zo2; xo2 yo2 zo2; ... xoN yoN zoN]; where N is the number of observation points.\n\
@end table\n\
\n\
 Array grid parameter:\n\
@table @code\n\
@item G\n\
  Vector of annulus radi [mm].\n\
@end table\n\
\n\
Sampling parameters: s_par = [dx dy dt nt]; \n\
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
Focusing method, options are: 'on', 'off', and 'ud'.\n\
@item  focal\n\
Focal distance [mm]. If foc_met = 'ud' (user defined) then focal is a vector of focusing delays.\n\
@end table\n\
\n\
Apodization parameters: apod_met, apod, and win_par. The apod_met (apodization method) options are:\n\
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
dream_arr_annu is an oct-function that is a part of the DREAM Toolbox available at\n\
@url{http://www.signal.uu.se/Toolbox/dream/}.\n\
\n\
Copyright @copyright{} 2006-2016 Fredrik Lingvall.\n\
@end deftypefn")
{
  double *ro,*s_par,*m_par;
  char   apod_met[50],foc_met[50];
  int    buflen;
  double dx,dy,dt;
  octave_idx_type    nt,no,n;
  double param=0,*delay,v,cp,alpha;
  int    isize=0;
  double *gr;
  int    ifoc=0;
  double focal=0,*ud_focal=NULL;
  double *apod=NULL;
  int    iweight=0,iapo=0;
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

    if (!((nrhs == 10) || (nrhs == 11))) {
    error("dream_arr_annu requires 10 or 11 input arguments!");
    return oct_retval;
  }
  else
    if (nlhs > 2) {
      error("Too many output arguments for dream_arr_annu!");
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
  // Grid function (position vector of the elements).
  //

  if ((mxGetM(1) > 1) & (mxGetN(1) > 1)) {
    error("Argument 2 must a vector (number of array elements)");
    return oct_retval;
  }
  isize = (int) mxGetM(1)*mxGetN(1); // Number of elementents in the array.
  const Matrix tmp1 = args(1).matrix_value();
  gr = (double*) tmp1.fortran_vec(); // Vector of annular radi.

  //
  // Temporal and spatial sampling parameters.
  //

  // Check that arg 3 is a 4 element vector
  if (!((mxGetM(2)==4 && mxGetN(2)==1) || (mxGetM(2)==1 && mxGetN(2)==4))) {
    error("Argument 3 must be a vector of length 4!");
    return oct_retval;
  }
  const Matrix tmp2 = args(2).matrix_value();
  s_par = (double*) tmp2.fortran_vec();
  dx    = s_par[0];		// Spatial x discretization size.
  dy    = s_par[1];		// Spatial dy iscretization size.
  dt    = s_par[2];		// Temporal discretization size (= 1/sampling freq).
  nt    = (octave_idx_type) s_par[3];	// Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //
  // Check that arg 4 is a scalar or a vector.
  if ( (mxGetM(3) * mxGetN(3) !=1) && ((mxGetM(3) * mxGetN(3)) != no)) {
    error("Argument 4 must be a scalar or a vector with a length equal to the number of observation points!");
    return oct_retval;
  }
  const Matrix tmp3 = args(3).matrix_value();
  delay = (double*) tmp3.fortran_vec();

  //
  // Material parameters
  //

  // Check that arg 5 is a 3 element vectora
  if (!((mxGetM(4)==3 && mxGetN(4)==1) || (mxGetM(4)==1 && mxGetN(4)==3))) {
    error("Argument 5 must be a vector of length 3!");
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

  //  ifoc = 1 - no foc, 2 foc xy.

  if (nrhs >= 6) {

    if (!mxIsChar(5)) {
      error("Argument 6 must be a string");
      return oct_retval;
    }

    std::string strin = args(5).string_value();
    buflen = strin.length();
    for ( n=0; n<=buflen; n++ ) {
      foc_met[n] = strin[n];
    }
    foc_met[buflen] = '\0';

    is_set = false;

    if (!strcmp(foc_met,"off")) {
      ifoc = 1;
      is_set = true;
    }

    if (!strcmp(foc_met,"on")) {
      ifoc = 2;
      is_set = true;
    }

    if (!strcmp(foc_met,"ud")) {
      ifoc = 6;
      is_set = true;

      if (mxGetM(6) * mxGetN(6) != isize ) {
        error("The time delay vector (argument 7) for user defined ('ud') focusing\n") ;
        error("delays must have the same length as the number of array elements.!");
        return oct_retval;
      }
      const Matrix tmp5 = args(6).matrix_value();
      ud_focal = (double*) tmp5.fortran_vec();
    }
    else {

      // Check that arg 7 is a scalar.
      if (mxGetM(6) * mxGetN(6) !=1 ) {
        error("Argument 7  must be a scalar!");
        return oct_retval;
      }
      // Focal point (in mm).
      const Matrix tmp5 = args(6).matrix_value();
      focal = (double) tmp5.fortran_vec()[0];
    }

    if (is_set == false) {
      error("Unknown focusing method!");
      return oct_retval;
    }

  } else
    ifoc = 1;

  //
  // Apodization.
  //

  // iweight = 1 - no apodization, 2  apodization.
  // iapo = 0 - user defined, 1 traingle, 2 Gauss, 3 raised cosine, 4 simply supported, 5 clamped.

  if (nrhs >= 8) {

    if (!mxIsChar(7)) {
      error("Argument 8 must be a string");
      return oct_retval;
    }

    std::string strin = args(7).string_value();
    buflen = strin.length();
    for ( n=0; n<=buflen; n++ ) {
      apod_met[n] = strin[n];
    }
    apod_met[buflen] = '\0';

    iweight = 1;			// default off.
    is_set = false;

    if (!strcmp(apod_met,"off")) {
      iweight = 1;
      is_set = true;
    }

    if (!strcmp(apod_met,"ud")) {
      iweight = 2;
      iapo = 0;
      is_set = true;

      // Vector of apodization weights.
      if (mxGetM(8) * mxGetN(8) != isize) {
        error("The length of argument 9 (apodization vector) must be the same as the number of array elements!");
        return oct_retval;
      }

      const Matrix tmp6 = args(8).matrix_value();
      apod = (double*) tmp6.fortran_vec();
    }

    if (!strcmp(apod_met,"triangle")) {
      iweight = 2;
      iapo = 1;
      is_set = true;
    }

    if (!strcmp(apod_met,"gauss")) {
      iweight = 2;
      iapo = 2;
      is_set = true;
    }

    if (!strcmp(apod_met,"raised")) {
      iweight = 2;
      iapo = 3;
      is_set = true;
    }

    if (!strcmp(apod_met,"simply")) {
      iweight = 2;
      iapo = 4;
      is_set = true;
    }

    if (!strcmp(apod_met,"clamped")) {
      iweight = 2;
      iapo = 5;
      is_set = true;
    }

    if (is_set == false) {
      error("Unknown apodization!");
      return oct_retval;
    }

    // Parameter for raised cos and Gaussian apodization functions.
    if (mxGetM(9) * mxGetN(9) !=1 ) {
      error("Argument 13 must be a scalar");
      return oct_retval;
    }

    const Matrix tmp7 = args(9).matrix_value();
    param = (double) tmp7.fortran_vec()[0];
  }
  else
    iweight = 1;

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
  if (nthreads > (unsigned int) no) {
    nthreads = no;
  }

  //
  // Error reporting.
  //

  if (nrhs == 11) {

    if (!mxIsChar(10)) {
      error("Argument 11 must be a string");
      return oct_retval;
    }

    std::string strin = args(10).string_value();
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

#ifdef USE_FFTW
  if (alpha != (double) 0.0)
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
    D[thread_n].dx = dx;
    D[thread_n].dy = dy;
    D[thread_n].dt = dt;
    D[thread_n].nt = nt;

    if (mxGetM(3) * mxGetN(3) == 1)
      D[thread_n].delay_method = SINGLE; // delay is a scalar.
    else
      D[thread_n].delay_method = MULTIPLE; // delay is a vector.

    D[thread_n].delay = delay;
    D[thread_n].v = v;
    D[thread_n].cp = cp;
    D[thread_n].alpha = alpha;
    D[thread_n].isize = isize;
    D[thread_n].gr = gr;
    D[thread_n].ifoc = ifoc;
    D[thread_n].iweight = iweight;
    D[thread_n].iapo = iapo;
    D[thread_n].focal = focal;
    D[thread_n].apod = apod;
    D[thread_n].param = param;
    D[thread_n].ud_focal = ud_focal;
    D[thread_n].h = h;
    D[thread_n].err_level = err_level;

    // Start the threads.
    threads[thread_n] = std::thread(smp_process, &D[thread_n]); // Start the threads.
    set_dream_thread_affinity(thread_n, nthreads, threads);
  }

  // Wait for all threads to finish.
  for (thread_n = 0; thread_n<nthreads; thread_n++)
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

#ifdef USE_FFTW
  if (alpha != (double) 0.0)
    att_close();
#endif

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
    Matrix err_mat(nt, no);
    err_p = err_mat.fortran_vec();
    err_p[0] = (double) out_err;
    oct_retval.append(err_mat);
  }

  return oct_retval;
}