/***
 *
 * Copyright (C) 2008,2009,2016 Fredrik Lingvall
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
#include <signal.h>

#include "das_arr.h"
#include "dream_error.h"

//
// Octave headers.
//

#include <octave/oct.h>

//
// Globals
//

int running;

//
// Function prototypes.
//

void sighandler(int signum);
void sig_abrt_handler(int signum);
void sig_keyint_handler(int signum);

//
// typedef:s
//

typedef void (*sighandler_t)(int);

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
 * das_arr Octave (oct) gateway function for DAS_ARR (delay-and-sum for arrays).
 *
 ***/

DEFUN_DLD (das_arr, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {}  [Y] =  das_arr(Ro,G,s_par,delay,m_par,foc_met,\n \
    focal,steer_met,steer_par,apod_met,apod,win_par,err_level);\n\
\n\
DAS_ARR Computes the delay reponse for array transducers. That is,\n\
DAS_ARR only computes the delay to each observation point which is\n\
represented by a '1' at the corresponding data point.\n\
\n\
Observation point(s) ([mm]):\n\
\n\
@table @code\n\
@item Ro\n\
An N x 3 matrix, Ro = [xo1 yo1 zo2; xo2 yo2 zo2; ... xoN yoN zoN]; where N is the number of observation points.\n\
@end table\n\
\n\
Array grid parameter:\n\
\n\
@table @code\n\
@item G\n\
An L x 3 grid function matrix. Column 1 contain the x-postions\n\
of the elements, column 2 the y-positions, column 3\n\
the z-positions, and where L is the number of elements.\n\
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
Material parameters: m_par = [cp];\n\
\n\
@table @code\n\
Sound velocity [m/s].\n\
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
das_arr is an oct-function that is a part of the DREAM Toolbox available at\n\
@url{http://www.signal.uu.se/Toolbox/dream/}.\n\
\n\
Copyright @copyright{} 2008-2019 Fredrik Lingvall.\n\
@seealso {das,saft,saft_p}\n\
@end deftypefn")
{
  double *ro,*s_par,*m_par;
  double *steer_par;
  char   apod_met[50];
  int    buflen;
  double xo,yo,zo,dt;
  dream_idx_type    nt,no,n;
  double param=0,*delay,cp;
  int    isize=0;
  double *gx,*gy,*gz;
  FocusMet foc_met=FocusMet::none;
  double *focal=nullptr;
  SteerMet steer_met=SteerMet::none;
  double theta=0,phi=0,*apod=NULL;
  int    iweight=0,iapo=0;
  double *h, *err_p;
  int    err_level=STOP, err=NONE, out_err = NONE, is_set = false;
  char   err_str[50];
  sighandler_t   old_handler, old_handler_abrt, old_handler_keyint;
  octave_value_list oct_retval;

  int nrhs = args.length ();

  // Check for proper number of arguments

  if (!((nrhs == 12) || (nrhs == 13))) {
    error("das_arr requires 12 or 13 input arguments!");
    return oct_retval;
  }
  else
    if (nlhs > 2) {
      error("Too many output arguments for das_arr !");
      return oct_retval;
    }

  //
  // Observation point.
  //

  // Check that arg 1 is a (number of observation points) x 3 matrix
  if ( mxGetN(0) != 3) {
    error("Argument 1 must be a (number of observation points) x 3 matrix!");
    return oct_retval;
  }

  no = mxGetM(0); // Number of observation points.
  const Matrix tmp0 = args(0).matrix_value();
  ro = (double*) tmp0.fortran_vec();

  //
  // Grid function (position vectors of the elements).
  //

  isize = (int) mxGetM(1); // Number of elementents in the array.
  if (mxGetN(1) !=3 )
    dream_err_msg("Argument 3  must a (number of array elements) x 3 matrix!");

  const Matrix tmp1 = args(1).matrix_value();
  gx = (double*) tmp1.fortran_vec(); // First column in the matrix.
  gy    = gx + isize;           // Second column in the matrix.
  gz    = gy + isize;           // Third column in the matrix.

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
  nt    = (dream_idx_type) s_par[1]; // Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //

  // Check that arg 4 is a scalar (or vector).

  DelayType delay_type=DelayType::single;

  if ( (mxGetM(3) * mxGetN(3) !=1) && ((mxGetM(3) * mxGetN(3)) != no)) {
    error("Argument 4 must be a scalar or a vector with a length equal to the number of observation points!");
    return oct_retval;
  } else {
    if (mxGetM(3) * mxGetN(3) ==1) {
      delay_type=DelayType::single;
    } else {
      delay_type=DelayType::multiple;
    }
  }

  const Matrix tmp3 = args(3).matrix_value();
  delay = (double*) tmp3.fortran_vec();

  //
  // Material parameters
  //

  // Check that arg 5 is a scalar.
  if (!(mxGetM(4)==1 && mxGetN(4)==1)) {
    error("Argument 5 must be scalar!");
    return oct_retval;
  }

  const Matrix tmp4 = args(4).matrix_value();
  m_par = (double*) tmp4.fortran_vec();
  cp    = m_par[0]; // Sound speed.

  //
  // Focusing parameters.
  //

  if (nrhs >= 6) {

    if (!mxIsChar(5)) {
      error("Argument 6 must be a string");
      return oct_retval;
    }

    std::string foc_str = args(5).string_value();

    is_set = false;

    if (foc_str == "off") {
      foc_met = FocusMet::none;
      is_set = true;
    }

    if (foc_str == "x") {
      foc_met = FocusMet::x;
      is_set = true;
    }

    if (foc_str == "y") {
      foc_met = FocusMet::y;
      is_set = true;
    }

    if (foc_str == "xy") {
      foc_met = FocusMet::xy;
      is_set = true;
    }

    if (foc_str == "x+y") {
      foc_met = FocusMet::x_y;
      is_set = true;
    }

    if (foc_str == "ud") {
      foc_met = FocusMet::ud;
      is_set = true;

      if (mxGetM(6) * mxGetN(6) != isize ) {
        error("The time delay vector (argument 7) for user defined ('ud') focusing\n") ;
        error("delays must have the same length as the number of array elements.!");
        return oct_retval;
      }
    } else {

      // Check that arg 7 is a scalar.
      if (mxGetM(6) * mxGetN(6) !=1 ) {
        error("Argument 7 to must be a scalar!");
        return oct_retval;
      }
    }

    if (is_set == false) {
      error("Unknown focusing method!");
      return oct_retval;
    }

  } else {
    foc_met = FocusMet::none;
  }

  const Matrix tmp6 = args(6).matrix_value();
  focal = (double*) tmp6.fortran_vec();

  //
  // Beam steering.
  //

  if (nrhs >= 8) {

    if (!mxIsChar(7)) {
      error("Argument 8 must be a string");
      return oct_retval;
    }

    std::string steer_str = args(7).string_value();

    steer_met = SteerMet::none;                  // Default no steering
    is_set = false;

    if (steer_str == "off") {
      steer_met = SteerMet::none;
      is_set = true;
    }

    if (steer_str == "x") {
      steer_met = SteerMet::x;
      is_set = true;
    }

    if (steer_str == "y") {
      steer_met = SteerMet::y;
      is_set = true;
    }

    if (steer_str == "xy") {
      steer_met = SteerMet::xy;
      is_set = true;
    }

    if (is_set == false) {
      error("Unknown beamsteering method!");
      return oct_retval;
    }

    // Check that arg 9 is a 2 element vector
    if (!((mxGetM(8)==2 && mxGetN(8)==1) || (mxGetM(8)==1 && mxGetN(8)==2))) {
      error("Argument 9 must be a vector of length 2!");
      return oct_retval;
    }

    const Matrix tmp5 = args(8).matrix_value();
    steer_par = (double*) tmp5.fortran_vec();
    theta  = steer_par[0];              // Angle in x-direction.
    phi    = steer_par[1];              // Angle in y-direction.
  } else {
    steer_met = SteerMet::none;
  }

  //
  // Apodization.
  //

  // iweight = 1 - no apodization, 2  apodization.
  // iapo = 0 - user defined, 1 traingle, 2 Gauss, 3 raised cosine, 4 simply supported, 5 clamped.

  if (nrhs >= 10) {

    if (!mxIsChar(9)) {
      dream_err_msg("Argument 10 must be a string");
      return oct_retval;
    }

    std::string strin = args(9).string_value();
    buflen = strin.length();
    for ( n=0; n<=buflen; n++ ) {
      apod_met[n] = strin[n];
    }
    apod_met[buflen] = '\0';

    iweight = 1;                        // default off.
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
      if (mxGetM(10) * mxGetN(10) != isize) {
        error("The length of argument 11 (apodization vector) must be the same as the number of array elements!");
        return oct_retval;
      }

      const Matrix tmp6 = args(10).matrix_value();
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
      error("Unknown apodization method!");
      return oct_retval;
    }

    // Parameter for raised cos and Gaussian apodization functions.
    if (mxGetM(11) * mxGetN(11) !=1 ) {
      error("Argument 12 must be a scalar");
      return oct_retval;
    }

    const Matrix tmp7 = args(11).matrix_value();
    param = (double) tmp7.fortran_vec()[0];

  }
  else
    iweight = 1;

  //
  // Error reporting.
  //
  if (nrhs == 13) {

    if (!mxIsChar(12)) {
      dream_err_msg("Argument 13 must be a string");
      return oct_retval;
    }

    std::string strin = args(12).string_value();
    buflen = strin.length();
    for ( n=0; n<=buflen; n++ ) {
      err_str[n] = strin[n];
    }
    err_str[buflen] = '\0';

    is_set = false;

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
    printf("Couldn't register SIGTERM  signal handler.\n");
  }

  if (( old_handler_abrt=signal(SIGABRT, &sighandler)) == SIG_ERR) {
    printf("Couldn't register SIGABRT signal handler.\n");
  }

  if (( old_handler_keyint=signal(SIGINT, &sighandler)) == SIG_ERR) {
    printf("Couldn't register SIGINT signal handler.\n");
  }

  //
  // Call the DAS subroutine.
  //

  running = true;

  for (n=0; n<no; n++) {

    xo = ro[n];
    yo = ro[n+1*no];
    zo = ro[n+2*no];

    double dlay = 0.0;
    if (delay_type == DelayType::single) {
      dlay = delay[0];
    } else { // DelayType::multiple.
      dlay = delay[n];
    }

    err = das_arr(xo, yo, zo,
                  dt, nt,
                  dlay,
                  cp,isize,
                  gx, gy, gz,
                  foc_met, focal,
                  steer_met, theta, phi,
                  apod,iweight,iapo,param,
                  &h[n*nt],err_level);

    if (err != NONE) {
      out_err = err;
      if (err == STOP) {
        error("");
        return oct_retval;
      }

      if (!running) {
        break; // CTRL-C pressed.
      }

    }
  }

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

  oct_retval.append(h_mat);

  // Return error.
  if (nlhs == 2) {
    Matrix err_mat(nt, no);
    err_p = err_mat.fortran_vec();
    err_p[0] = (double) out_err;
    oct_retval.append(err_mat);
  }

  return oct_retval;
}
