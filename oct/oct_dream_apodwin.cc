/***
*
* Copyright (C) 2006,2007,2008,2009,2015,2016,2021,2021 Fredrik Lingvall
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

#include "arr_functions.h"
#include "arg_parser.h"

//
// Octave headers.
//

#include <octave/oct.h>

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
 Inputparameters:\n\
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
The options for the @code{apod_met}parameter are:\n\
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
@url{https://github.com/frli8848/DREAM}.\n\
\n\
Copyright @copyright{} 2006-2021 Fredrik Lingvall.\n\
@seealso {dream_arr_rect, dream_arr_circ, dream_arr_cylind_f, dream_arr_cylind_d}\n\
@end deftypefn")
{
  bool do_apod=false;
  ApodMet apod_met=ApodMet::gauss;
  double *h;
  octave_value_list oct_retval;

  dream_idx_type nrhs = args.length ();

  ArgParser ap;

  if (!ap.check_arg_in("dream_apodwin", nrhs, 3, 3)) {
    return oct_retval;
  }

  if (!ap.check_arg_out("dream_apodwin", nlhs, 0, 1)) {
    return oct_retval;
  }

  //
  // Apodization.
  //

  const Matrix tmp1 = args(1).matrix_value();
  dream_idx_type num_elements = (dream_idx_type) tmp1.fortran_vec()[0];

  // Allocate space for the user defined apodization weights
  std::unique_ptr<double[]> apod = std::make_unique<double[]>(num_elements);

  double apod_par;
  if (!ap.parse_apod_args("dream_apodwin", args, 0, num_elements,
                            do_apod, apod.get(), apod_met, apod_par)) {
    return oct_retval;
  }

  // Create an output matrix for the impulse response
  Matrix h_mat(num_elements,1);
  h = h_mat.fortran_vec();

  double weight=1.0;
  double ramax = 1.0;
  double ys = 0.0;
  if (do_apod) {
    for (dream_idx_type n=0; n<num_elements; n++) {
      double xs = 2*ramax * (0.5 - ((double) n / (double) num_elements));
      apodization(apod_met, n, apod.get(), &weight, xs, ys, ramax, apod_par);
      h[n] = weight;
    }
  }

  oct_retval.append(h_mat);
  return oct_retval;
}
