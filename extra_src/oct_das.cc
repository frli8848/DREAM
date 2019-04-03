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
#include "das.h"
#include "dream_error.h"

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

/***
 *
 * das - Octave (oct) gateway function for DAS (delay-and-sum).
 *
 ***/

DEFUN_DLD (das, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {}  [Y] = das(Ro,s_par,delay,m_par).\n\
\n\
DAS Computes the delay reponse for single element transducer. That is,\n\
DAS only comutes the delay to each observation point which is\n\
represented by a '1' at the corresponding data point.\n\
\n\
Observation point(s) ([mm]):\n\
\n\
@table @code\n\
@item Ro\n\
An N x 3 matrix, Ro = [xo1 yo1 zo2; xo2 yo2 zo2; ... xoN yoN zoN]; where N is the number of observation points.\n\
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
@item cp\n\
Sound velocity [m/s].\n\
\n\
@end table\n\
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
das is an oct-function that is a part of the DREAM Toolbox available at\n\
@url{http://www.signal.uu.se/Toolbox/dream/}.\n\
\n\
Copyright @copyright{} 2008-2019 Fredrik Lingvall.\n\
@seealso {das_arr,saft,saft_p}\n\
@end deftypefn")
{
  double *ro,*s_par,*m_par;
  dream_idx_type nt,no,n;
  double xo,yo,zo,dt;
  double *delay,cp;
  double *h, *err_p;
  int    err_level=STOP, err=NONE, out_err = NONE, is_set = false;
  char   err_str[50];
  int    buflen;
  octave_value_list oct_retval;

  int nrhs = args.length ();

  // Check for proper number of arguments

  if (!((nrhs == 5) || (nrhs == 6))) {
    error("das requires 5 or 6 input arguments!");
    return oct_retval;
      }
  else
    if (nlhs > 2) {
      error("Too many output arguments for das!");
      return oct_retval;
    }

  //
  // Observation point.
  //

  // Check that arg (number of observation points) x 3 matrix.
  if (mxGetN(0) != 3) {
    error("Argument 1 must be a (number of observation points) x 3 matrix!");
    return oct_retval;
  }

  no = mxGetM(0); // Number of observation points.
  const Matrix tmp0 = args(0).matrix_value();
  ro = (double*) tmp0.fortran_vec();

  //
  // Temporal and spatial sampling parameters.
  //

  if (!((mxGetM(1)==2 && mxGetN(2)==1) || (mxGetM(1)==1 && mxGetN(1)==2))) {
    error("Argument 2 must be a vector of length 2!");
    return oct_retval;
  }

  const Matrix tmp1 = args(1).matrix_value();
  s_par = (double*) tmp1.fortran_vec();
  dt    = s_par[0]; // Temporal discretization size (= 1/sampling freq).
  nt    = (dream_idx_type) s_par[1];	// Length of SIR.

  //
  // Start point of impulse response vector ([us]).
  //

  // Check that arg 3 is a scalar (or vector).
    if ( (mxGetM(2) * mxGetN(2) !=1) && ((mxGetM(2) * mxGetN(2)) != no)) {
    error("Argument 3 must be a scalar or a vector with a length equal to the number of observation points!");
    return oct_retval;
  }

  const Matrix tmp2 = args(2).matrix_value();
  delay = (double*) tmp2.fortran_vec();

  //
  // Material parameters
  //

 // Check that arg 4 is a scalar.
  if (!(mxGetM(3)==1 && mxGetN(3)==1)) {
    error("Argument 4 must be a scalar!");
    return oct_retval;
  }

  const Matrix tmp3 = args(3).matrix_value();
  m_par = (double*) tmp3.fortran_vec();
  cp    = m_par[0]; // Sound speed.

  // Error reporting.
  if (nrhs == 5) {

    if (!mxIsChar(4)) {
      error("Argument 5 must be a string");
      return oct_retval;
    }

    std::string strin = args(4).string_value();
    buflen = strin.length();
    for ( n=0; n<=buflen; n++ ) {
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

  // Create an output matrix for the impulse response.
  Matrix h_mat(nt, no);
  h = h_mat.fortran_vec();

  // Call the DAS subroutine.
  if (mxGetM(2) * mxGetN(2) == 1) {
    for (n=0; n<no; n++) {
      xo = ro[n];
      yo = ro[n+1*no];
      zo = ro[n+2*no];

      err = das(xo,yo,zo,dt,nt,delay[0],cp,&h[n*nt],err_level);

      if (err != NONE) {
        out_err = err;
        if (err == STOP) {
          error("");
          return oct_retval;
        }
      }

    }
  } else {
    for (n=0; n<no; n++) {
      xo = ro[n];
      yo = ro[n+1*no];
      zo = ro[n+2*no];

      err = das(xo,yo,zo,dt,nt,delay[n],cp,&h[n*nt],err_level);

      if (err != NONE) {
        out_err = err;
        if (err == STOP) {
          error("");
          return oct_retval;
        }
      }

    }
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
