/***
*
* Copyright (C) 2008,2009,2013,2014,2015,2016,2021,2021 Fredrik Lingvall
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
#include <cmath>
#include <thread>
#include <mutex>

#include "dream.h"
#include "affinity.h"
#include "dream_error.h"

//
// Octave headers.
//

#include <octave/oct.h>

/***
 *
 * (Monostatic) Synthetic Aperture Focusing Technique - Parallel version.
 *
 ***/

//
// Globals
//

int running;

//
// typedef:s
//

typedef struct
{
  int    no;
  int    start;
  int    stop;
  double *ro;
  double dt;
  DelayType delay_type;
  double *delay;
  double cp;
  double a;
  int    K;
  int    L;
  double *r_trans;
  double *B;
  double *Bsaft;
} DATA;

typedef void (*sighandler_t)(int);

//
// Function prototypes.
//

void* smp_dream_saft(void *arg);
void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

/***
 *
 * Thread function.
 *
 ***/

void* smp_dream_saft(void *arg)
{
  octave_idx_type n, l;
  int  k_shift;			// Can be negative!
  DATA D = *(DATA *)arg;
  double *Bsaft = D.Bsaft;
  double *B = D.B;
  double dt=D.dt;
  int    no=D.no, K=D.K, L=D.L;
  double *delay=D.delay, *ro=D.ro, cp=D.cp;
  double xo, yo, zo, r_xy;
  octave_idx_type  start=D.start, stop=D.stop;
  double *r_trans = D.r_trans, a = D.a;
  double x_trans, y_trans, z_trans, d, z, z_prim, tmp;

  if (D.delay_type == DelayType::single) {
    for (n=start; n<stop; n++) {
      xo = ro[n];
      yo = ro[n+1*no];
      zo = ro[n+2*no];

      // The SAFT loop.
      for (l=0; l < (octave_idx_type) L; l++) {

        x_trans = r_trans[l];
        y_trans = r_trans[l+1*L];
        z_trans = r_trans[l+2*L];

        // Horizontal distance between transducer and observation point.
        r_xy = std::sqrt( (x_trans - xo)*(x_trans - xo) + (y_trans - yo)*(y_trans - yo) );

        // if r_xy is inside the synthetic aperture.
        if (r_xy <= a/2) {

          // Vertical distance between transducer and observation point.
          z = zo - z_trans;

          z_prim = std::sqrt(z*z + r_xy*r_xy);
          tmp = (2*z_prim/cp * 1e3 - 2*delay[0])/dt;

          // Better to round just one time!
          k_shift = (int) rint(tmp);

          // Rounding err.
          d =  tmp - ((double) k_shift);

          // Linear interpolation.
          if ((k_shift+1 < K) && (k_shift-1 > 0) ) {
            if (d >=0) {
              Bsaft[n] += (1.0 - d) * B[k_shift     + l*K];
              Bsaft[n] += d * B[(k_shift+1) + l*K];
            }
            else {
              Bsaft[n] -= d * B[k_shift     + l*K];
              Bsaft[n] += (1.0 + d) * B[(k_shift-1) + l*K];
            }
          } // if
        }
      }

      if (!running) {
        octave_stdout << "Thread for observation points " << start+1 << " -> " << stop  << " bailing out!\n";
        return(NULL);
      }

    }

  } else { // DelayType::multiple.

    for (n=start; n<stop; n++) {
      xo = ro[n];
      yo = ro[n+1*no];
      zo = ro[n+2*no];

      // The SAFT loop.
      for (l=0; l<L; l++) {

        x_trans = r_trans[l];
        y_trans = r_trans[l+1*L];
        z_trans = r_trans[l+2*L];

        // Horizontal distance between transducer and observation point.
        r_xy = std::sqrt( (x_trans - xo)*(x_trans - xo) + (y_trans - yo)*(y_trans - yo) );

        // if r_xy is inside the synthetic aperture.
        if (r_xy <= a/2) {

          // Vertical distance between transducer and observation point.
          z = zo - z_trans;

          z_prim = std::sqrt(z*z + r_xy*r_xy);
          tmp = (2*z_prim/cp * 1e3 - 2*delay[n])/dt;

          // Better to round just one time!
          k_shift = (int) rint(tmp);

          // Rounding err.
          d =  tmp - ((double) k_shift);

          // Linear interpolation.
          if ((k_shift+1 < K) && (k_shift-1 > 0) ) {
            if (d >=0) {
              Bsaft[n] += (1.0 - d) * B[k_shift     + l*K];
              Bsaft[n] += d * B[(k_shift+1) + l*K];
            }
            else {
              Bsaft[n] -= d * B[k_shift     + l*K];
              Bsaft[n] += (1.0 + d) * B[(k_shift-1) + l*K];
            }
          } // if
        }
      }

      if (!running) {
        octave_stdout << "Thread for observation points " << start+1 << " -> " << stop  << " bailing out!\n";
        return(NULL);
      }
    }
  }

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
 * Octave (oct) gateway function for SAFT_P.
 *
 ***/

DEFUN_DLD (saft, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {}  [Y] = saft(B,To,delay,s_par,m_par,Ro,a).\n\
\n\
SAFT_P - The sythetic aperture focusing techique.\n\
\n\
Observation point(s) ([mm]):\n\
\n\
@table @code\n\
@item B\n\
An K x L matrix B-scan data matrix.\n\
@end table\n\
\n\
@table @code\n\
@item To\n\
An L x 3 matrix, Ro = [xo1 yo1 zo1; xo2 yo2 zo2; ... xoL yoL zoL]; where L is the number of transducer positions.\n\
@end table\n\
\n\
Start point of SIR:\n\
\n\
@table @code\n\
@item  delay\n\
Scalar delay for all observation points or a vector with individual delays for each observation point [us].\n\
@end table\n\
\n\
Sampling parameters: s_par = [dt]; \n\
\n\
@table @code\n\
@item dt\n\
Temporal discretization period (= 1/sampling freq) [us].\n\
@end table\n\
\n\
Material parameters: m_par = [v cp];\n\
\n\
@table @code\n\
@item cp\n\
Sound velocity [m/s].\n\
@end table\n\
\n\
@table @code\n\
@item Ro\n\
An N x 3 matrix, Ro = [xo1 yo1 zo1; xo2 yo2 zo2; ... xoN yoN zoN]; where N is the number of observation points.\n\
@end table\n\
\n\
@table @code\n\
@item  a\n\
The synthetic aperture [mm].\n\
@end table\n\
\n\
saft is an oct-function that is a part of the DREAM Toolbox available at\n\
@url{http://www.signal.uu.se/Toolbox/dream/}.\n\
\n\
Copyright @copyright{} 2008-2021 Fredrik Lingvall.\n\
@seealso {saft,das,das_arr}\n\
@end deftypefn")
{
  double *s_par, *m_par;
  octave_idx_type K, L, no;
  double *B, *Bsaft, cp , a, dt, *ro;
  double *r_trans, *delay;
  sighandler_t old_handler, old_handler_abrt, old_handler_keyint;
  octave_idx_type start, stop;
  std::thread *threads;
  unsigned int thread_n, nthreads;
  DATA   *D;
  octave_value_list oct_retval;

  int nrhs = args.length ();

  // Check for proper number of arguments

  if (nrhs != 7) {
    error("saft requires 7 input arguments!");
    return oct_retval;
  }
  else
    if (nlhs > 1) {
      error("Too many output arguments for saft!");
      return oct_retval;
    }

  //
  // B-scan.
  //

  K = mxGetM(0); // Number of temporal samples.
  L = mxGetN(0); // Number of transducer positions.
  const Matrix tmp0 = args(0).matrix_value();
  B = (double*) tmp0.fortran_vec();

  //
  // Transducer position matrix.
  //

  // Check that arg 2 is a L x 3 is a matrix.
  if (!(mxGetM(1)==L && mxGetN(1)==3)) {
    error("Arg 2 must be a (number of transducer positions) x 3 matrix!");
    return oct_retval;
  }
  const Matrix tmp1 = args(1).matrix_value();
  r_trans = (double*) tmp1.fortran_vec();

  //
  // Start point of impulse response vector ([us]).
  //

  // Check that arg 3 is a scalar (or vector).
  if ( (mxGetM(2) * mxGetN(2) !=1) && ((mxGetM(2) * mxGetN(2)) != L)) {
    error("Argument 3 (delay(s)) must be a scalar or a vector with a length equal to the number of A-scans!");
    return oct_retval;
  }

  const Matrix tmp2 = args(2).matrix_value();
  delay = (double*) tmp2.fortran_vec();

  //
  // Sampling parameter.
  //

  // Check that arg 4 is a scalar.
  if (!((mxGetM(3)==1 && mxGetN(3)==1))) {
    error("Argument 4 (temporal sampling period) must be a scalar");
    return oct_retval;
  }

  const Matrix tmp3 = args(3).matrix_value();
  s_par = (double*) tmp3.fortran_vec();
  dt    = s_par[0]; // Temporal discretization size (= 1/sampling freq).

  //
  // Material parameter.
  //

  // Check that arg 5 is a scalar.
  if (!((mxGetM(4)==1 && mxGetN(4)==1))) {
    dream_err_msg("Argument 5 (sound speed) must be a scalar!");
    return oct_retval;
  }

  const Matrix tmp4 = args(4).matrix_value();
  m_par = (double*) tmp4.fortran_vec();
  cp    = m_par[0]; // Sound speed.

  //
  // Observation point matrix.
  //

  // Check that arg 6 (number of observation points) x 3 matrix.
  if (mxGetN(5) !=3 ) {
    error("Argument 6 must be a (number of observation points) x 3 matrix!");
    return oct_retval;
  }

  no = mxGetM(5); // Number of observation points.
  const Matrix tmp5 = args(5).matrix_value();
  ro = (double*) tmp5.fortran_vec();

  //
  // Synthetic Aperture
  //

  // Check that arg 7 is scalar.
  if ((mxGetM(6)!=1) && (mxGetN(6)!=1)) {
    error("Argument 7 must be a scalar!");
    return oct_retval;
  }

  const Matrix tmp6 = args(6).matrix_value();
  a = (double) tmp6.fortran_vec()[0];

  //
  // Number of threads.
  //

  // Get number of CPU cores (including hypethreading, C++11)
  nthreads = std::thread::hardware_concurrency();

  // nthreads can't be larger then the number of observation points.
  if (nthreads > (unsigned int) no) {
    nthreads = no;
  }

  //
  // Create an output matrix for the processed image.
  //

  Matrix h_mat(no,1);
  Bsaft = h_mat.fortran_vec();

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
  // Call the SAFT subroutine.
  //

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
    D[thread_n].cp = cp;
    D[thread_n].dt = dt;

    if (mxGetM(3) * mxGetN(3) == 1)
      D[thread_n].delay_type = DelayType::single; // delay is a scalar.
    else
      D[thread_n].delay_type = DelayType::multiple; // delay is a vector.

    D[thread_n].delay = delay;
    D[thread_n].B = B;
    D[thread_n].Bsaft = Bsaft;
    D[thread_n].a = a;
    D[thread_n].K = K;
    D[thread_n].L = L;
    D[thread_n].r_trans = r_trans;

    // Start the threads.
    threads[thread_n] = std::thread(smp_dream_saft, &D[thread_n]);
    set_dream_thread_affinity(thread_n, nthreads, threads);
  }

  // Wait for all threads to finish.
  for (thread_n = 0; thread_n < nthreads; thread_n++) {
    threads[thread_n].join();
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
    error("CTRL-C pressed!\n"); // Bail out.
    return oct_retval;
  }

  oct_retval.append(h_mat);

  return oct_retval;
}
