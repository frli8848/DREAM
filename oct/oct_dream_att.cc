/***
*
* Copyright (C) 2006,2007,2008,2009,2012,2015,2016 Fredrik Lingvall
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
#include <math.h>
#include "att.h"
#include "dream_error.h"

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


/***
 *
 * Octave (OCT) gateway function for the attenuation code on att.c
 *
 ***/

DEFUN_DLD (dream_att, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} [H] = dream_att(Ro,s_par,delay,m_par);\n\
\n\
DREAM_ATT - Computes attenuation impulse response(s) using the same method\n\
as the transducer functions in the DREAM Toolbox.\n\
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
Material parameters: m_par = [cp alpha];\n\
\n\
@table @code\n\
@item cp\n\
Sound velocity [m/s].\n\
@item alpha\n\
Attenuation coefficient [dB/(cm MHz)] .\n\
\n\
@end table\n\
\n\
dream_att is an oct-function that is a part of the DREAM Toolbox available at\n\
@url{http://www.signal.uu.se/Toolbox/dream/}.\n\
\n\
Copyright @copyright{} 2006-2019 Fredrik Lingvall.\n\
@seealso {all transducer functions: dreamline, dreamrect, dreamcirc, etc}\n\
@end deftypefn")
{
  double *ro,*s_par,*m_par;
  int    it,nt,no,n;
  double xo,yo,zo,dt;
  double *delay,cp,alpha,r;
  double *h;
  octave_value_list oct_retval;

  int nrhs = args.length ();

  // Check for proper number of arguments

  if (nrhs != 4) {
    error("dream_att requires 4 input arguments!");
    return oct_retval;
  }
  else
    if (nlhs > 1) {
      error("dream_att requires one output argument!");
      return oct_retval;
    }

  //
  // Observation point.
  //

  // Check that arg (number of observation points) x 3 matrix
  if ( mxGetN(0) != 3) {
    dream_err_msg("Argument 1 must be a (number of observation points) x 3 matrix!");
    return oct_retval;
  }
  no = mxGetM(0); // Number of observation points.
  const Matrix tmp0 = args(0).matrix_value();
  ro = (double*) tmp0.fortran_vec();

  //
  // Temporal and spatial sampling parameters.
  //


  if (!((mxGetM(1)==2 && mxGetN(1)==1) || (mxGetM(1)==1 && mxGetN(1)==2))) {
    error("Argument 2 must be a vector of length 2!");
    return oct_retval;
  }
  const Matrix tmp1 = args(1).matrix_value();
  s_par = (double*) tmp1.fortran_vec();
  dt    = s_par[0];		// Temporal discretization size (= 1/sampling freq).
  nt    = (int) s_par[1];	// Length of SIR.

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

  // Check that arg 4 is a 2 element vector.
 if (!((mxGetM(3)==2 && mxGetN(3)==1) || (mxGetM(3)==1 && mxGetN(3)==2))) {
    error("Argument 4 must be a vector of length 2!");
    return oct_retval;
  }
  const Matrix tmp3 = args(3).matrix_value();
  m_par = (double*) tmp3.fortran_vec();
  cp    = m_par[0]; // Sound speed.
  alpha  = m_par[1]; // Attenuation coefficient [dB/(cm MHz)].

  //
  // Create an output matrix for the impulse response.
  //

  // Create an output matrix for the impulse response(s)
  Matrix h_mat(nt, no);
  h = h_mat.fortran_vec();

  // Clear data (it looks like h_mat do not get out of
  // scope between calls!?).
  for (n=0; n<nt*no; n++)  {
    h[n] =  0.0;
  }

  {
    Attenuation att(nt, dt, cp, alpha);
    FFTCVec xc(nt);
    FFTVec x(nt);

    //
    // Call the attenuation subroutine.
    //

    if (mxGetM((2)) * mxGetN(2) == 1) {
      for (n=0; n<no; n++) {
        xo = ro[n];
        yo = ro[n+1*no];
        zo = ro[n+2*no];

        r = sqrt(xo*xo + yo*yo + zo*zo);
        it = (int) ( (r * 1.0e3/cp - delay[0])/dt + 1);
        att.att(xc, x, r, it, &h[n*nt], 1.0);
      }
    } else {
      for (n=0; n<no; n++) {
        xo = ro[n];
        yo = ro[n+1*no];
        zo = ro[n+2*no];

        r = sqrt(xo*xo + yo*yo + zo*zo);
        it = (int) ( (r*1.0e3/cp - delay[n])/dt + 1);
        att.att(xc, x, r, it, &h[n*nt], 1.0);
      }
    }
  }

  oct_retval.append(h_mat);
  return oct_retval;
}
