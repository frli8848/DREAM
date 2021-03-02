/***
*
* Copyright (C) 2003,2006,2007,2008,2009,2014,2021 Fredrik Lingvall
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

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "circ_sir.h"
#include "dream_error.h"

/***
 *
 *  Subroutine circ_sir - Continous spatial impulse response for a circular aperture.
 *
 *  See "Transient Radiation fom Pistons in an Infinite Planar Baffle", P.R. Stepanishen,
 *  Journal of the Acoustical Society of America}, volume={49}, pages={1629--38}, year={1971}.
 *
 ***/

void circ_sir(double xo_i,
              double yo_i,
              double zo_i,
              double r_i,
              double dt,
              dream_idx_type nt,
              double delay,
              double v,
              double cp,
              double *h)
{
  dream_idx_type it;
  double t, r_cyl, xo, yo, zo, r;
  double R, Rprim;

  // Convert to [m].
  xo = xo_i / 1000.0;
  yo = yo_i / 1000.0;
  zo = zo_i / 1000.0;
  r = r_i / 1000.0;

  for (it = 0; it < nt; it++) {
    h[it] = (double) 0.0;
  }

  r_cyl = sqrt(xo*xo + yo*yo); // Distance from z-axis.
  Rprim = sqrt(zo*zo + (r-r_cyl)*(r-r_cyl) );
  R     = sqrt(zo*zo + (r+r_cyl)*(r+r_cyl) );

  for (it=0; it<nt; it++) {

    t = (((double) it) * dt + delay)/1e6; // in [s].

    if (r > r_cyl) {

      if (cp*t < zo) {
        h[it] = 0.0;
      }
      else if ( (cp*t >= zo) && (cp*t < Rprim) ) {
        h[it] = v * cp;
      }
      else if ( (cp*t >= Rprim) && (cp*t < R) ) {
        h[it] = v * cp/M_PI * acos( (cp*t*cp*t - zo*zo + r_cyl*r_cyl - r*r) / (2*r_cyl*sqrt(cp*t*cp*t-zo*zo)));
      }
      else if (cp*t >= R) {
        h[it] = 0.0;
      }

    }
    else {  // r <=  r_cyl

      if (cp*t < Rprim) {
        h[it] = 0.0;
      }
      else if ( (cp*t >= Rprim) && (cp*t < R) ) {
        h[it] = v * cp/M_PI * acos( (cp*t*cp*t - zo*zo + r_cyl*r_cyl - r*r) / (2*r_cyl*sqrt(cp*t*cp*t-zo*zo)));
      }
      else if (cp*t >= R) {
        h[it] = 0.0;
      }

    }

  }
  return;
}
