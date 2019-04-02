/***
*
* Copyright (C) 2002,2003,2005,2006,2007,2008,2009,2014 Fredrik Lingvall
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

#include "att.h"
#include "dream_arr_annu.h"
#include "dream_error.h"

#if defined(_MSC_VER) || defined(__LCC__)
#include "msvc_rint.h"
#endif


//
// Function prototypes
//

void center_pos_annular(double *RESTRICT rs, double *RESTRICT gr, int isize, int nv, double *RESTRICT ramax);
void superpos_annular(double *RESTRICT hi, double *RESTRICT ha, dream_idx_type  nt, double weight, double retfoc, dream_idx_type j, double  dt);
void focusing_annular(int ifoc, double focal, double rs, double ramax, double cp, double *RESTRICT retfoc);
void apodization_annular(int iweight, int iapo, double *RESTRICT apod, double *RESTRICT weight, double rs,
                       double ramax, double param, dream_idx_type i);
int circ_annular(double xo, double  yo, double  zo, double  a, double dx, double dy, double dt,
                 dream_idx_type nt, double delay, double v, double cp, double alpha, double weight, double *RESTRICT h,
                 dream_idx_type k, int isize,int err_level);
void distance(double xo, double yo, double zo, double xs, double ys, double zs, double *RESTRICT ri);
void xlimit(double yi, double a, double *RESTRICT xsmin, double *RESTRICT xsmax);
void resp_annular(double *RESTRICT h, double *RESTRICT hi, dream_idx_type nt, dream_idx_type j);

/***
*
*     Routine for computation of spatial impulse response of an annular array.
*
***/

int dream_arr_annu(double xo, double yo, double zo, double dx, double dy, double dt,
                    dream_idx_type  nt, double delay, double v, double cp, double alpha,
                    int isize, double *RESTRICT gr, int ifoc, double focal, double *RESTRICT apod, int iweight,
                    int iapo, double param,double *RESTRICT ha,int err_level)
{
  double r, *h;
  double *RESTRICT hi;
  dream_idx_type i, j;
  double ramax;
  int    nv;
  double *RESTRICT rs, retfoc = 0, weight = 0;
  int err = NONE, out_err = NONE;

  h  = (double*) malloc( nt*isize * sizeof(double));
  hi = (double*) malloc( nt*isize * sizeof(double));
  rs = (double*) malloc(isize*sizeof(double));

  for (i=0; i< nt; i++) {
    ha[i] = (double) 0.0;
  }

  retfoc = (double) 0.0;
  weight = (double) 1.0;
  // nv - number of annulus.
  nv = (isize+1)/2;
  //nv = (isize)/2; // This is wrong!!

  for (i=0; i<nt; i++) {
    for (j=0; j<isize; j++) {
      h[i+j*nt] = (double) 0.0;
    }
  }

  /* ----------------------------------------- */

  for (i=0; i<isize; i++) {
    r = gr[i];

    if (i==0 && r == (double) 0.0 ) {
      for (j=0; j<nt; j++) {
        h[j] = (double) 0.0; // h(j,1) = 0.0;
      }
    } else {
      err = circ_annular(xo,yo,zo,r,dx,dy,dt,nt,delay,v,cp,alpha,weight,h,i,isize,err_level);
      if (err != NONE)
        out_err = err;
    }
  }

  /* --------------------------------------------- */

  center_pos_annular(rs, gr, isize, nv, &ramax);

  for (i=0; i<nv; i++) {
    focusing_annular(ifoc, focal, rs[i], ramax, cp, &retfoc);
    apodization_annular(iweight, iapo, apod, &weight, rs[i], ramax, param, i);
    resp_annular(h, hi, nt, i);
    superpos_annular(hi, ha, nt, weight, retfoc, i, dt);
  }

  free(h);
  free(hi);
  free(rs);

  return out_err;
} /* dream_arr_annu */


/***
 *
 * Routine for computation of spatial impulse response of an annular array -
 * user defined focusing.
 *
 ***/

int dream_arr_annu_ud(double xo, double yo, double zo, double dx, double dy, double dt,
                    dream_idx_type  nt, double delay, double v, double cp, double alpha,
                    int isize, double *RESTRICT gr, int ifoc, double *RESTRICT focal, double *RESTRICT apod, int iweight,
                    int iapo, double param, double *RESTRICT ha,int err_level)
{
  double r, *h;
  double *RESTRICT hi;
  dream_idx_type i, j;
  double ramax;
  int    nv;
  double *RESTRICT rs, retfoc = 0, weight = 0;
  int err = NONE, out_err = NONE;

  h  = (double*) malloc( nt*isize * sizeof(double));
  hi = (double*) malloc( nt*isize * sizeof(double));
  rs = (double*) malloc(isize*sizeof(double));

  for (i=0; i< nt; i++)
    ha[i] = (double) 0.0;

  retfoc = (double) 0.0;
  weight = (double) 1.0;

  // nv - number of annulus.
  nv = (isize+1)/2;
  //nv = (isize)/2; // This is wrong!!

  for (i=0; i<nt; i++) {
    for (j=0; j<isize; j++) {
      h[i+j*nt] = (double) 0.0;
    }
  }

  /* ----------------------------------------- */

  for (i=0; i<isize; i++) {
    r = gr[i];

    if (i==0 && r == (double) 0.0 ) {
      for (j=0; j<nt; j++) {
        h[j] = (double) 0.0; // h(j,1) = 0.0;
      }
    } else {
      err = circ_annular(xo,yo,zo,r,dx,dy,dt,nt,delay,v,cp,alpha,weight,h,i,isize,err_level);
      if (err != NONE)
        out_err = err;
    }
  }

  /* --------------------------------------------- */

  center_pos_annular(rs, gr, isize, nv, &ramax);

  for (i=0; i<nv; i++) {
    focusing_annular(ifoc, focal[i], rs[i], ramax, cp, &retfoc);  // Note ifoc must be 6 here!
    apodization_annular(iweight, iapo, apod, &weight, rs[i], ramax, param, i);
    resp_annular(h, hi, nt, i);
    superpos_annular(hi, ha, nt, weight, retfoc, i, dt);
  }

  free(h);
  free(hi);
  free(rs);

  return out_err;
} /* dream_arr_annu */


/***
 *
 * subroutine center_pos_annular(xs,ys,zs,i,j) pour calculer le centrn element
 *
 ***/

void center_pos_annular(double *RESTRICT rs, double *RESTRICT gr, int isize, int nv, double *RESTRICT ramax)
{
  int i, ns = 0;

  rs[0] = gr[0]/2;
  for (i=2; i<isize; i += 2) {
    ns = (i+1) / 2;
    rs[ns] = (gr[i]+gr[i-1]) / 2;
  }
  *ramax = rs[nv-1];

  return;
} /* center_pos_annular */


/***
*
* subroutine focussing gives le retard retfoc du au focussing
*
***/

void focusing_annular(int ifoc, double focal, double rs, double ramax, double cp, double *RESTRICT retfoc)
{
  double diff, rmax;

  // ifoc = 1 - no foc, 2 foc xy.
  if (ifoc == 1) {
    *retfoc = 0.0;
    return;
  }

  if (ifoc == 2) {
    rmax = sqrt(ramax*ramax + focal*focal);
    diff = rmax - sqrt(rs*rs + focal*focal);

    *retfoc = diff * 1000 / cp;
  }

  // User defined focusing.
  if (ifoc == 6)
    *retfoc = focal; // Here focal is the user defined time delay in [us] (not the focal depth)

  return;
} /* focusing_annular */

/***
 *
 *  apodization
 *
 * apodization iapo = 0 apodization with imported apodisation function apod(x,y)
 *
 * iweight = 1 - no apodization, 2  apodization , param=input parameter
 *
 * iapo = 0 - user defined.
 * iapo = 1 traingle.
 * iapo = 2 gauss.
 * iapo = 3 rised cosine
 * iapo = 4 simply supported.
 * iapo = 5 clamped.
 *
 ***/

void apodization_annular(int iweight, int iapo, double *RESTRICT apod, double *RESTRICT weight, double rs,
                       double ramax, double param, dream_idx_type i)
{
  double pi = atan((double) 1.0) * (double) 4.0;

  if (iweight == 1) {
    return;
  }

  switch(iapo) {

  case 0:
    *weight = apod[i];
    break;

  case 1:
    *weight = (double) 1.0 - fabs(rs) / ramax;
    break;

  case 2:
    *weight = exp(-(param * rs*rs) / (ramax*ramax));
    break;

  case 3:
    *weight = param + cos(rs * pi / ramax);
    break;

  case 4:
    *weight = (double) 1.0 - rs*rs / (ramax*ramax);
    break;

  case 5:
    *weight = ((double) 1.0 - rs*rs / (ramax*ramax)) *
      ((double) 1.0 - rs*rs / (ramax*ramax));
    break;

  default:
    break;
  }

  return;
} /* apodization_annular */

/***
 *
 * subroutine circ - pour calculer pulse respone of a circular aperture
 *
 * a = R1,R2,...,Rk,...,Risize = gr(k)
 *
 ***/

int circ_annular(double xo, double  yo, double  zo, double  r, double dx, double dy, double dt,
                  dream_idx_type nt, double delay, double v, double cp, double alpha, double weight, double *RESTRICT h,
                  dream_idx_type k, int isize, int err_level)
{
  double t;
  double xsmin, ysmin, xsmax, ysmax, ai, ds, pi, ri;
  dream_idx_type it;
  double zs, x, y;
  int err = NONE;

  pi = atan( (double) 1.0) * 4.0;
  ds = dx * dy;
  zs = (double) 0.0;

  ysmin = -r;
  ysmax =  r;

  //j = 0;
  //j++;
  //y = ysmin + (j-1) * dy + dy/2;
  y = ysmin + dy/2.0;
  while (y <= ysmax) {

    xlimit(y, r, &xsmin, &xsmax);

    //i = 0;
    //i++;
    //x = xsmin + (i-1) * dx + dx/2;
    x = xsmin + dx/2.0;
    while (x <= xsmax) {

      distance(xo, yo, zo, x, y, zs, &ri);
      ai = weight * v * ds / (2*pi*ri);
      ai /= dt;
      // Convert to SI units.
      ai *= 1000;
      // Propagation delay in micro seconds.
      t = ri * 1000 / cp;
      it = (dream_idx_type) rint((t - delay)/dt);

      if ((it < nt) && (it >= 0))

        if (alpha == (double) 0.0) {
          h[it + k*nt] += ai;
        }
        else {
          att_annu(alpha, ri, it, dt, cp, h, nt, ai, k, isize);
        }
      else  {
        if  (it >= 0)
          err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
        else
          err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

        if ( (err_level == PARALLEL_STOP) || (err_level == STOP) )
          return err; // Bail out.
      }

      //i++;
      //x = xsmin + (i-1)*dx + dx/2;
      x += dx;
    }
    //j++;
    //y = ysmin + (j-1)*dy + dy/2;
    y += dy;
  }

  return err;
} /* circ_annular */


/***
 *
 * subrutine distance(xi,xs,hs,ri,rx,rz) pour trouver le longeur du vecteur
 *
 ***/

void distance(double xo, double yo, double zo, double xs, double ys, double zs, double *RESTRICT ri)
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
 *  subroutine xlimit - pour definir les limits d integration en x
 *
 ***/

void xlimit(double yi, double a, double *RESTRICT xsmin, double *RESTRICT xsmax)
{
  double rs;

  rs = sqrt(a*a - yi*yi);
  *xsmin = -rs;
  *xsmax = rs;

  return;
} /* xlimit */


/***
 *
 * call superpos(h,ha) subroutine pour superposer les contributions des elements
 *
 * ha = output response
 * h  = input responce of actual element
 ***/

void superpos_annular(double *RESTRICT hi, double *RESTRICT ha, dream_idx_type  nt, double weight, double retfoc, dream_idx_type j, double  dt)
{
  double *RESTRICT buff;
  dream_idx_type    i,it1;

  buff = (double*) malloc(2*nt*sizeof(double));

  it1 = (dream_idx_type) (retfoc / dt) + 1;

  for (i=0; i<2*nt; i++)
    buff[i] = (double) 0.0;

  for (i=0; i<nt; i++)
    buff[i+it1] = hi[i+j*nt];

  for (i=0; i<nt; i++)
    ha[i] += weight * buff[i];

  free(buff);

  return;
} /* superpos_annular */


/***
 *
 * resp_annular
 *
 ***/

void resp_annular(double *RESTRICT h, double *RESTRICT hi, dream_idx_type nt, dream_idx_type j)
{
  dream_idx_type i, k;

  k = 2*j;

  if (k == 0) {
    for (i=0; i<nt; i++) { // Center element.
      // hi(i,j) = h(i,1)
      hi[i+j*nt] = h[i];
    }
  }
  else {
    for (i=0; i <nt; i++) { // Ring # j.
      // hi(i,j) = h(i,k)-h(i,k-1)
      hi[i+j*nt] = h[i+k*nt] - h[i+(k-1)*nt];
    }
  }

  return;
} /* resp_annular */
