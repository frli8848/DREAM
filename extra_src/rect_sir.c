/***
*
* Copyright (C) 2004,2006,2007,2008,2009,2014 Fredrik Lingvall
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

// $Revision: 771 $ $Date: 2014-05-16 09:32:48 +0200 (Fri, 16 May 2014) $ $LastChangedBy: frli8848 $

#include <math.h>
#include <string.h>
#include <stdio.h>
#include "rect_sir.h"
#include "dream_error.h"

/***
 *
 *  Subroutine rect_sir - Continous spatial impulse response for a rectangular aperture.
 *
 *  See "High-speed method for computing the exact solution for the pressure variations
 *  in the nearfield of a baffled piston", J.C. Lockwood and J.G. Willette, J. Acoust. Soc. Am.,
 *  year = {1973},  volume =  {53},  number = {3},  pages =  {735--741}, month =  {March}.
 *
 ***/

void rect_sir(double xo_i,
              double yo_i,
              double zo_i,
              double a_i,
              double b_i,
              double dt,
              dream_idx_type    nt,
              double delay,
              double v,
              double cp,
              double *RESTRICT h)
{
  dream_idx_type    it, k;
  double t, t_z, xo, yo, zo;
  double          pi;
  double          tau_1, tau_2, tau_3, tau_4, a, b;
  double          a_k=0, g_k=0, s_k=0, l_k=0;

  pi = atan( (double) 1.0) * 4.0;

  // Convert to [m].
  xo = fabs(xo_i) / 1000.0;	// Can take abs due to symmetry.
  yo = fabs(yo_i) / 1000.0;	// Can take abs due to symmetry.
  zo = zo_i / 1000.0;
  a = a_i / 1000;
  b = b_i / 1000;

  for (it = 0; it < nt; it++) {
    h[it] = (double) 0.0;
  }

  for  (k=1; k<=4; k++) { // loop over all 4 sub-rectangles.

    if ( (xo <= a/2) && (yo <= b/2) ) {// In shadow.

      g_k = 1; // Add all.
      switch (k) {

      case 1:
        s_k = fabs(a/2 - xo);
        l_k = fabs(b/2 - yo);
        break;

      case 2:
        s_k = fabs(xo + a/2);
        l_k = fabs(b/2 - yo);
        break;

      case 3:
        s_k = fabs(xo + a/2);
        l_k = fabs(yo + b/2);
        break;

      case 4:
        s_k = fabs(a/2 - xo);
        l_k = fabs(b/2 - yo);
        break;

      default:
        break;
      }
    } // if in shadow.

    if ( (xo <= a/2) && (yo > b/2) ) {// Inside a/2 but outside b/2.

      switch (k) {

      case 1:
        s_k = fabs(a/2 - xo);
        l_k = fabs(yo - b/2);
        g_k = -1; // Outside => subtract.
        break;

      case 2:
        s_k = fabs(a/2 + xo);
        l_k = fabs(yo - b/2);
        g_k = -1; // Outside => subtract.
        break;

      case 3:
        s_k = fabs(a/2 + xo);
        l_k = fabs(b/2 + yo);
        g_k = 1; // Inside => add.
        break;

      case 4:
        s_k = fabs(a/2 - xo);
        l_k = fabs(b/2 + yo);
        g_k = 1; // Inside => add.
        break;

      default:
        break;
      }
    } // if inside a/2 but outside b/2.

    if ( (xo > a/2) && (yo <= b/2) ) { // Inside b/2 but outside a/2.

      switch (k) {

      case 1:
        s_k = fabs(xo - a/2);
        l_k = fabs(b/2 - yo);
        g_k = -1; // Outside => subtract.
        break;

      case 2:
        s_k = fabs(a/2 + xo );
        l_k = fabs(b/2 - yo);
        g_k = 1; // Inside => add.
        break;

      case 3:
        s_k = fabs(a/2 + xo);
        l_k = fabs(b/2 + yo);
        g_k = 1; // Inside => add.
        break;

      case 4:
        s_k = fabs(xo - a/2);
        l_k = fabs(b/2 + yo);
        g_k = -1; // Outside => subtract.
        break;

      default:
        break;
      }
    } // if inside b/2 but outside a/2.

    if ( (xo > a/2) && (yo > b/2) ) {// Outside both  a/2 and b/2.
      switch (k) {

      case 1:
        s_k = fabs(xo - a/2);
        l_k = fabs(yo - b/2);
        g_k = 1; // Really outside but need to add since 1 is included both in 2 and 4.
        break;

      case 2:
        s_k = fabs(a/2 + xo);
        l_k = fabs(yo - b/2);
        g_k = -1; // Outside => subtract.
        break;

      case 3:
        s_k = fabs(a/2 + xo);
        l_k = fabs(b/2 + yo);
        g_k = 1; // Inside => add.
        break;

      case 4:
        s_k = fabs(xo - a/2);
        l_k = fabs(yo + b/2);
        g_k = -1; // Outside => subtract.
        break;

      default:
        break;
      }
    } // if outside both  a/2 and b/2.

    t_z = zo / cp;

    tau_1 = t_z;
    tau_2 = sqrt(zo*zo + s_k*s_k) / cp;
    tau_3 = sqrt(zo*zo + l_k*l_k) / cp;
    tau_4 = sqrt(zo*zo + l_k*l_k + s_k*s_k) / cp;

    for (it=0; it<nt; it++) {

      t = (((double) it) * dt + delay)/1e6; // in [s].

      a_k = 0;
      if ( (t >= tau_1) && (t <= tau_4) )
        a_k += cp/4;

      if ( (t >= tau_2) && (t <= tau_4) )
        a_k -= cp/(2*pi) * acos( s_k / (cp * sqrt(t*t - t_z*t_z)));

      if ( (t >= tau_3) && (t <= tau_4) )
        a_k -= cp/(2*pi) * acos( l_k / (cp * sqrt(t*t - t_z*t_z)));

      h[it] += g_k * a_k;

    } // for it
  } // for k

  return;
}
