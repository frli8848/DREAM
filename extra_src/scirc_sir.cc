/***
*
* Copyright (C) 2003,2006,2007,2008,2009,2014 Fredrik Lingvall
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
#include "scirc_sir.h"
#include "dream_error.h"

//
// Function prototypes.
//

double acos_acos(double c1,double c2);
double trapets_int(double t1, double t2, double zo, double r, double r_cyl, double cp,dream_idx_type int_len);

/***
 *
 *  Subroutine scirc_sir - Sampled spatial impulse response for a circular aperture.
 *
 * Sampling is performed by integrating, using the trapets formula, the continous spatial
 * impulse response for one sampling interval (zero-order hold sampling).
 *
 * The continous SIR can be found in "Transient Radiation fom Pistons in an Infinite Planar Baffle",
 * P.R. Stepanishen, Journal of the Acoustical Society of America}, volume={49}, pages={1629--38},
 * year={1971}.
 *
 ***/

void scirc_sir(double xo_i,
               double yo_i,
               double zo_i,
               double r_i,
               double dt,
               dream_idx_type    nt,
               double delay,
               double v,
               double cp,
               dream_idx_type int_len,
               double *RESTRICT h)
{
  dream_idx_type    it;
  double t1, t2, r_cyl, xo, yo, zo, r;
  double R, Rprim;
  double pi;
  pi = atan( (double) 1.0) * 4.0;

  // Convert to [m].
  xo = xo_i / 1000.0;
  yo = yo_i / 1000.0;
  zo = zo_i / 1000.0;
  r  = r_i  / 1000.0;

  for (it = 0; it < nt; it++) {
    h[it] = (double) 0.0 ;
  }

  r_cyl = sqrt(xo*xo + yo*yo); // Distance from z-axis.
  Rprim = sqrt(zo*zo + (r-r_cyl)*(r-r_cyl) );
  R     = sqrt(zo*zo + (r+r_cyl)*(r+r_cyl) );

  for (it=0; it<nt; it++) {

    // Sampling interval [t1,t2] = t +/- dt/2.
    t1 = ( (((double) it) * dt) + delay)/1.0e6; // in [s].
    t2 = t1 + dt/1.0e6;
    t1 -= (dt/1.0e6)/2.0;
    t2 -= (dt/1.0e6)/2.0;

    if (t1 < 0.0)
      t1 = 0.0;

    if (r > r_cyl) {

      if ( (cp*t1 < zo) && (cp*t2 < zo)) {
        //
        // t1,t2 < zo/cp.
        //
        h[it] = 0.0;
      }
      else if ( (cp*t1 < zo) && (cp*t2 >= zo) && (cp*t2 < Rprim) ) {
        //
        // t1 < zo/cp, t2 in [zo/cp,Rprim/cp).
        //
        h[it] = v * cp * (t2 - zo/cp);

      }
      else if ( (cp*t1 < zo) && (cp*t2 >= Rprim) && (cp*t2 < R) ) {
        //
        // t1 < zo/cp and t2 in [Rprim/cp,R/cp).
        //
        h[it] = 0.0;
        h[it] += v * cp * (Rprim - zo)/cp;

        if (R > Rprim)
          h[it] += v * cp/pi * trapets_int(Rprim/cp, t2, zo, r, r_cyl, cp,int_len);

      }
      else if ( (cp*t1 < zo) && (cp*t2 >= R) ) {
        //
        // t1 < zo/cp and t2 > R/cp.
        //
        h[it] = v * cp * (Rprim - zo)/cp;

        if (R > Rprim)
          h[it] += v * cp/pi * trapets_int(Rprim/cp, R/cp, zo, r, r_cyl, cp,int_len);

      } // **********************************
      else if ( (cp*t1 >= zo ) && (cp*t2 < Rprim) ) {
        //
        // t1,t2 in [zo/cp,Rprim/cp].
        //
        h[it] = v * cp * (t2-t1);

      }
      else if ( (cp*t1 >= zo) && (cp*t1 < Rprim) && (cp*t2 >= Rprim) && (cp*t2 < R) ) {
        //
        // t1 in [zo/cp,Rprim/cp) and t2 in [Rprim/cp,R/cp).
        //
        h[it] = v * cp * (Rprim/cp - t1);

        if (R > Rprim)
          h[it] += v * cp/pi * trapets_int(Rprim/cp, t2, zo, r, r_cyl, cp,int_len);

      }
      else if ( (cp*t1 >= zo) && (cp*t1 < Rprim) && (cp*t2 >= R) ) {
        //
        // t1 in [zo/cp,Rprim/cp) and t2 > R/cp.
        //
        h[it] = v * cp * (Rprim/cp - t1);

        if (R > Rprim)
          h[it] += v * cp/pi * trapets_int(Rprim/cp, R/cp, zo, r, r_cyl, cp,int_len);

      } // **********************************
      else if ( (cp*t1 >= Rprim) && (cp*t2 < R) ) {
        //
        // t1,t2 in [Rprim/cp,R/cp).
        //

        if (R > Rprim)
          h[it] += v * cp/pi * trapets_int(t1, t2, zo, r, r_cyl, cp,int_len);

      }
      else if ( (cp*t1 >= Rprim) &&  (cp*t1 < R) && (cp*t2 >= R) ) {
        //
        // t1 in [Rprim/cp,R/cp) and t2 > R/cp.
        //
        if (R > Rprim)
          h[it] += v * cp/pi * trapets_int(t1, R/cp, zo, r, r_cyl, cp,int_len);

      } // **********************************
      else if (cp*t1 > R) {
        //
        // t1,t2 > R/cp.
        //
        h[it] = 0.0;
      }
    }
    else {  // %%%%%%%%%%%%%%%%%% r <=  r_cyl %%%%%%%%%%%%%%%%%%%%%%

      //test = 0;

      //if (r_cyl <= 0.0)
      //	printf("Error!\n");

      if (cp*t2 < Rprim) {
        //
        // t1,t2 < Rprim/cp.
        //
        h[it] = 0.0;
      }
      else if ( (cp*t1 < Rprim) && (cp*t2 >= Rprim) && (cp*t2 < R) ) {
        //
        // t1 < Rprim/cp and t2 in [Rprim/cp,R/cp).
        //

        if (R > Rprim)
          h[it] += v * cp/pi * trapets_int(Rprim/cp, t2, zo, r, r_cyl, cp,int_len);

      }
      else if ( (cp*t1 < Rprim) && (cp*t2 >= R) ) {
        //
        // t1 < Rprim/cp and t2 > R/cp.
        //

        if (R > Rprim)
          h[it] += v * cp/pi * trapets_int(Rprim/cp, R/cp, zo, r, r_cyl, cp,int_len);

      } // **********************************
      else if ( (cp*t1 >= Rprim) && (cp*t2 < R) ) {
        //
        // t1,t2 in [Rprim/cp, Rprim/cp).
        //

        if (R > Rprim)
          h[it] += v * cp/pi * trapets_int(t1, t2, zo, r, r_cyl, cp,int_len);

      }
      else if ( (cp*t1 >= Rprim) && (cp*t1 < R) && (cp*t2 >= R) ) {
        //
        // t1 in [Rprim/cp, Rprim/cp) and t2 > R/cp.
        //

        if (R > Rprim)
          h[it] += v * cp/pi * trapets_int(t1, R/cp, zo, r, r_cyl, cp,int_len);

      }
      else if (cp*t1 >= R) {
        //
        // t1,t2 > R/cp.
        //
        h[it] = 0.0;
      }

    }

    h[it] /= dt/1.0e6; // Bogdan's "normalization".

  }

  return;
}


//
// Function to avoid numrical problem for acos near +/- 1.
//

double acos_acos(double c1,double c2)
{
  double ac1=0.0, ac2=0.0;
  double pi = atan( (double) 1.0) * 4.0;

  if (c1 >= -1.0 && c1 <= 1.0) {
    ac1 = acos(c1);
  }
  else if (c1 < -1.0) {
    ac1 = pi;
  }
  else if (c1 > 1.0) {
    ac1 = 0.0;
  }

  if (c2 >= -1.0 && c2 <= 1.0) {
    ac2 = acos(c2);
  }
  else if (c2 < -1.0) {
    ac2 = pi;
  }
  else if (c2 > 1.0) {
    ac2 = 0.0;
  }

  return ac1 + ac2;
}

//
// Use trapets rule to approximate the integral.
//

double trapets_int(double t1, double t2, double zo, double r, double r_cyl, double cp, dream_idx_type int_len)
{
  double dt = (t2 - t1)/((double) int_len), t, y;
  double c1 = 0.0, c2 = 0.0;
  dream_idx_type n;

  y = 0.0;
  t = t1;
  for (n = 0; n < int_len; n++) {

    c1 = (cp*t*cp*t - zo*zo + r_cyl*r_cyl - r*r) / (2*r_cyl*sqrt(cp*t*cp*t-zo*zo));
    c2 = (cp*(t+dt)*cp*(t+dt) - zo*zo + r_cyl*r_cyl - r*r)  / (2*r_cyl*sqrt(cp*(t+dt)*cp*(t+dt)-zo*zo));

    y += acos_acos(c1,c2)/2 * dt;

    t += dt;
  }

  return y;
}
