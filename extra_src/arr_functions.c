/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2014 Fredrik Lingvall
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
#include <stdio.h>
#include <stdlib.h>

#include "arr_functions.h"

/***
 *
 * Subroutines for common for array transducers.
 *
 *
 *
 ***/

/***
*
*  center_pos(xs,ys,zs,i,j) pour calculer le centrn element
*
***/

void center_pos(double *RESTRICT xs, double *RESTRICT ys, double *RESTRICT zs, int i,
                double *RESTRICT gx, double *RESTRICT gy, double *RESTRICT gz)
{
  *xs = gx[i];
  *ys = gy[i];
  *zs = gz[i];

  return;
} /* center_pos */


/***
 *
 *  max_dim_arr : Give les valeurs maximales d'aperture de l'array
 *
 *   Computes the maximum aperture of the array.
 *
 ***/

void max_dim_arr(double *RESTRICT xamax, double *RESTRICT yamax, double *RESTRICT ramax,
               double *RESTRICT gx, double *RESTRICT gy, double *RESTRICT gz, int isize)
{
  int i;
  double ret;

  *xamax = fabs(gx[0]);
  for (i=1; i<isize; i++) {
    ret =  fabs(gx[i]);
    if (ret > *xamax) {
      *xamax = ret;
    }
  }

  *yamax = fabs(gy[0]);
  for (i=1; i<isize; i++) {
    ret = fabs(gy[i]);
    if (ret > *yamax) {
      *yamax = ret;
    }
  }
  *ramax = sqrt(*xamax * *xamax + *yamax * *yamax);

  return;
} /* max_dim_arr */

/***
 *
 *   Focussing gives le retard retfoc du au focussing
 *
 *   focusing - computes the  focusing delay, retfoc, for point (xs,ys) on the
 *              transducer surface (?).
 *
 ***/

void focusing(int ifoc, double focal, double xs, double ys,
              double xamax, double yamax, double ramax, double cp, double *RESTRICT retfoc)
{
  double diff, rmax, retx, rety;

  //
  // ifoc = 1 No foc, 2 Foc x ,3 Foc y, 4 Foc xy 5 Foc x+y, 6 ud (user defined) */
  //

  switch(ifoc) {

  case 1:
    return;

  case 2:
    rmax = sqrt(xamax*xamax + focal*focal);
    diff = rmax - sqrt(xs*xs + focal*focal);
    *retfoc = diff*1000/cp;
    break;

  case 3:
    rmax = sqrt(yamax*yamax + focal*focal);
    diff  = rmax - sqrt(ys*ys + focal*focal);
    *retfoc = diff*1000/cp;
    break;

  case 4:
    rmax = sqrt(ramax*ramax + focal*focal);
    diff = rmax - sqrt(xs*xs + ys*ys + focal*focal);
    *retfoc = diff*1000/cp;
    break;

  case 5:
    rmax = sqrt(ramax*ramax + focal*focal);
    retx = sqrt(xs*xs + focal*focal);
    rety = sqrt(ys*ys + focal*focal);
    diff = rmax - (retx + rety);
    *retfoc = diff*1000/cp;
    break;

  case 6:
    *retfoc = focal; // Here focal is the user defined time delay in [us] (not the focal depth).
    break;

  default:
    break;

  }

  return;
} /* focusing */



/***
 *
 * Subroutine for beam steering
 *
 ***/

void beamsteering(int ister, double theta, double phi, double xs, double ys,
                  double xamax, double yamax, double ramax, double cp, double *RESTRICT retsteer)
{
  double diff, rmax, sinx, siny, retsteerx, retsteery, pi, pii;

  pi = atan((double) 1.0) * (double) 4.0;
  pii = pi / (double) 180.0;

  //
  // ister = 1 No steering, 2 Steer x ,3 Steer y, 4 Steer xy.
  //

  switch(ister) {

  case 1:
    break;

  case 2:
    sinx = sin(theta * pii);
    rmax = xamax * sinx;
    diff =  rmax + xs*sinx;
    *retsteer = diff*1000/cp;
    break;

  case 3:
    siny = sin(phi * pii);
    rmax = yamax * siny;
    diff = rmax + ys * siny;
    *retsteer = diff*1000/cp;
    break;

  case 4:
    sinx = sin(theta * pii);
    rmax = xamax * sinx;
    diff = rmax + xs*sinx;
    retsteerx = diff*1000/cp;
    siny = sin(phi * pii);
    rmax = yamax * siny;
    diff = rmax + ys * siny;
    retsteery = diff*1000/cp;
    *retsteer = retsteerx + retsteery;
    break;

  default:
    break;
  }

  return;
} /* beamsteering */

/***
 *
 *  Weighting
 *
 * apodization iapo = 0 apodization with imported apodization function apod(x,y)
 *
 *      iweight = 1 No apodization, 2  Weighting , param=input parameter
 *
 * iapo = 0 User defined
 * iapo = 1 Traingle
 * iapo = 2 Gauss
 * iapo = 3 Rised cosine
 * iapo = 4 Simply supported
 * iapo = 5 Clamped
 *
 ***/

void apodization(int iweight, int iapo, int i, double *RESTRICT apod, double *RESTRICT weight,
               double xs, double ys, double ramax, double param, int isize)
{
  static double r, pi;


  pi = atan((double) 1.0) * (double)4.;
  r = sqrt(xs*xs + ys*ys);

  if (iweight == 1) {
    return;
  }

  switch(iapo) {

  case 0:
    *weight = apod[i];
    break;

  case 1:
    *weight = (double) 1.0 - fabs(r) / ramax;
    break;

  case 2:
    *weight = exp(-(param * r*r) / (ramax*ramax));
    break;

  case 3:
    *weight = param + cos(r*pi/ramax);
    break;

  case 4:
    *weight = (double) 1.0 - r*r / (ramax*ramax);
    break;

  case 5:
    *weight = ((double) 1.0 - r*r / (ramax*ramax)) * ((double) 1.0  - r / (ramax*ramax));

  default:
    break;
  }

  return;
} /* apodization */

/***
 *
 * Subrutine distance(xi,xs,hs,ri,rx,rz) pour trouver le longeur du vecteur
 *
 ***/

void distance(double xo, double yo, double zo,double xs,double ys, double zs, double *RESTRICT ri)
{
  double rx, ry, rz;

  rx = xo - xs;
  ry = yo - ys;
  rz = zo - zs;
  *ri = sqrt(rx*rx + rz*rz + ry*ry);

  return;
} /* distance */

/***
 *
 * superpos(h,ha) Subroutine pour superposer les contributions des elements
 *
 * ha = output response
 * h  = input responce of actual element
 *
 ***/

void superpos(double *RESTRICT h, double *RESTRICT ha, dream_idx_type nt)
{
  dream_idx_type i;

  for (i=0; i< nt; i++)
    ha[i] += h[i];

  return;
} /* superpos */


/***
 *
 * suobroitine pour check le retard de la reponse
 *
 ***/

void check_delay(dream_idx_type it, double tt, int *icheck, dream_idx_type nt)
{
  if (tt < (double) 0.0)
    *icheck = 2;

  if (it > nt)
    *icheck = 2;

  return;
} /* check_delay */