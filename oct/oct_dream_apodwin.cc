/***
*
* Copyright (C) 2006,2007,2008,2009,2015,2016,2019 Fredrik Lingvall
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
#include <stdio.h>
#include "arr_functions.h"
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
 * Octave (oct) gateway function for dream_apodwin (for the apodization functions in
 * arr_functions.c)
 *
 ***/

DEFUN_DLD (dream_apodwin, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} [w] = dream_apodwin(apod_met,apod,win_par)\n\
\n\
DREAM_APODWIN - Computes apodization using the same method\n\
as used in the array functions in the DREAM Toolbox.\n\
\n\
 Input parameters:\n\
\n\
@table @code\n\
@item apod_met\n\
The apodization method (a text string),\n\
@item apod\n\
Is the sample length of the apodization window (or a vector)\n\
of apodization weights for the 'ud' option),\n\
@item win_par\n\
Parameter for raised cosine and Gaussian apodization functions (a scalar).\n\
@end table\n\
\n\
The options for the @code{apod_met} parameter are:\n\
\n\
@table @code\n\
@item 'off'\n\
No apodizatio,\n\
@item 'ud'\n\
Used defined apodization.\n\
@item 'triangle'\n\
Triangle window,\n\
@item 'gauss'\n\
Gaussian (bell-shaped) window.\n\
@item 'raised'\n\
Raised cosine,\n\
@item 'simply'\n\
Simply supported.\n\
@item 'clamped'\n\
Clamped.\n\
@end table\n\
\n\
dream_apodwin is an oct-function that is a part of the DREAM Toolbox available at\n\
@url{http://www.signal.uu.se/Toolbox/dream/}.\n\
\n\
Copyright @copyright{} 2006-2019 Fredrik Lingvall.\n\
@seealso {dream_arr_rect, dream_arr_circ, dream_arr_cylind_f, dream_arr_cylind_d}\n\
@end deftypefn")
{
  int    iweight=0, apod_type=0, i, num_elements=0, is_set = false;
  double *apod=nullptr, weight, xs, ys, ramax, param;
  double *h;
  octave_value_list oct_retval;

  int nrhs = args.length ();

  // Check for proper number of arguments
  if (nrhs != 3) {
    error("dream_apodwin requires three input arguments!");
    return oct_retval;
  }
  else
    if (nlhs > 1) {
      error("dream_apodwin requires one output argument!");
      return oct_retval;
    }

  //
  // Apodization.
  //

  // apod_type = 0 - user defined, 1 traingle, 2 Gauss, 3 raised cosine, 4 simply supported, 5 clamped.

  if (!mxIsChar(0)) {
    error("Argument 1 must be a string");
    return oct_retval;
  }

  std::string apod_str = args(0).string_value();

  iweight = 1;			// default off.
  is_set = false;

  if (apod_str == "off") {
    iweight = 1;
    is_set = true;
  }

  if (apod_str == "ud") {
    iweight = 2;
    apod_type = APOD_UD;
    error(" 'ud'- (user defined) meaningless for this function!");
    return oct_retval;
  }

  if (apod_str == "triangle") {
    iweight = 2;
    apod_type = APOD_TRIANGLE;
    is_set = true;
  }

  if (apod_str == "gauss") {
    iweight = 2;
    apod_type = APOD_GAUSS;
    is_set = true;
  }

  if (apod_str == "raised") {
    iweight = 2;
    apod_type = APOD_RISED_COSINE;
    is_set = true;
  }

  if (apod_str == "simply") {
    iweight = 2;
    apod_type = APOD_SIMPLY_SUPPORTED;
      is_set = true;
  }

  if (apod_str == "clamped") {
    iweight = 2;
    apod_type = APOD_CLAMPED;
    is_set = true;
  }

  if (is_set == false) {
    error("Unknown apodization!");
    return oct_retval;
  }

  //
  // Number of elements.
  //

  if (mxGetM(1) * mxGetN(1) !=1) {
    dream_err_msg("Argument 2 must be a scalar!");
    return oct_retval;
  }
  const Matrix tmp1 = args(1).matrix_value();
  num_elements = (int) tmp1.fortran_vec()[0];

  if (num_elements < 0) {
    error("Argument 2 must be a positive integer!");
    return oct_retval;
  }

  //
  // param - Parameter used for raised cos and Gaussian apodization functions.
  //

  if (mxGetM(2) * mxGetN(2) !=1) {
    error("Argument 3 must be a scalar!");
    return oct_retval;
  }
  const Matrix tmp2 = args(2).matrix_value();
  param = (double) tmp2.fortran_vec()[0];

  // Create an output matrix for the impulse response
  Matrix h_mat(num_elements,1);
  h = h_mat.fortran_vec();

  ramax = 1;
  ys = 0;
  if (iweight != 1) {
    for (i=0; i<num_elements; i++) {
      xs = 2*ramax * (0.5 - ((double) i / (double) num_elements));
      apodization(apod_type, i, apod, &weight, xs, ys, ramax, param);
      h[i] = weight;
    }
  }

  oct_retval.append(h_mat);
  return oct_retval;
}
