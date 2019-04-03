/***
*
* Copyright (C) 2002,2003,2005,2006,2007,2008,2009,2014,2019 Fredrik Lingvall
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
#include "arr_functions.h"
#include "dream_arr_annu.h"
#include "dream_error.h"

#if defined(_MSC_VER) || defined(__LCC__)
#include "msvc_rint.h"
#endif

//
// Function prototypes
//

void center_pos_annular(double *RESTRICT rs, double *RESTRICT gr, int num_elements, int nv, double *RESTRICT ramax);
void superpos_annular(double *RESTRICT hi, double *RESTRICT ha, dream_idx_type  nt, double weight, double retfoc, dream_idx_type j, double  dt);
void focusing_annular(int foc_type, double focal, double rs, double ramax, double cp, double *RESTRICT retfoc);
void apodization_annular(int apod_type, int i, double *RESTRICT apod, double *RESTRICT weight, double rs,
                         double ramax, double param);
int circ_annular(double xo, double  yo, double  zo, double  a, double dx, double dy, double dt,
                 dream_idx_type nt, double delay, double v, double cp, double alpha, double weight, double *RESTRICT h,
                 dream_idx_type k, int num_elements,int err_level);
void xlimit_annular(double yi, double a, double *RESTRICT xsmin, double *RESTRICT xsmax);
void resp_annular(double *RESTRICT h, double *RESTRICT hi, dream_idx_type nt, dream_idx_type j);

/***
*
* dream_arr_annu
*
* Routine for computation of spatial impulse response of an annular array.
*
***/

int dream_arr_annu(double xo, double yo, double zo, double dx, double dy, double dt,
                    dream_idx_type  nt, double delay, double v, double cp, double alpha,
                    int num_elements, double *RESTRICT gr, int foc_type, double focal, double *RESTRICT apod, bool do_apod,
                    int apod_type, double param,double *RESTRICT ha,int err_level)
{
  double r, *h;
  double *RESTRICT hi;
  dream_idx_type i, j;
  double ramax;
  int    nv;
  double *RESTRICT rs, retfoc = 0, weight = 0;
  int err = NONE, out_err = NONE;

  h  = (double*) malloc( nt*num_elements * sizeof(double));
  hi = (double*) malloc( nt*num_elements * sizeof(double));
  rs = (double*) malloc(num_elements*sizeof(double));

  for (i=0; i< nt; i++) {
    ha[i] = 0.0;
  }

  retfoc = 0.0;
  weight = 1.0;
  // nv - number of annulus.
  nv = (num_elements+1)/2;
  //nv = (num_elements)/2; // This is wrong!!

  for (i=0; i<nt; i++) {
    for (j=0; j<num_elements; j++) {
      h[i+j*nt] = (double) 0.0;
    }
  }

  /* ----------------------------------------- */

  for (i=0; i<num_elements; i++) {
    r = gr[i];

    if (i==0 && r == (double) 0.0 ) {
      for (j=0; j<nt; j++) {
        h[j] = (double) 0.0; // h(j,1) = 0.0;
      }
    } else {
      err = circ_annular(xo,yo,zo,r,dx,dy,dt,nt,delay,v,cp,alpha,weight,h,i,num_elements,err_level);
      if (err != NONE)
        out_err = err;
    }
  }

  /* --------------------------------------------- */

  center_pos_annular(rs, gr, num_elements, nv, &ramax);

  for (i=0; i<nv; i++) {
    focusing_annular(foc_type, focal, rs[i], ramax, cp, &retfoc);
    if (do_apod){
      apodization_annular(apod_type, i, apod, &weight, rs[i], ramax, param);
    }
    resp_annular(h, hi, nt, i);
    superpos_annular(hi, ha, nt, weight, retfoc, i, dt);
  }

  free(h);
  free(hi);
  free(rs);

  return out_err;
}

/***
 *
 * dream_arr_annu_ud
 *
 * Routine for computation of spatial impulse response of an annular array -
 * user defined focusing.
 *
 ***/

int dream_arr_annu_ud(double xo, double yo, double zo, double dx, double dy, double dt,
                      dream_idx_type  nt, double delay, double v, double cp, double alpha,
                      int num_elements, double *RESTRICT gr, int foc_type, double *RESTRICT focal, double *RESTRICT apod, bool do_apod,
                      int apod_type, double param, double *RESTRICT ha,int err_level)
{
  double r, *h;
  double *RESTRICT hi;
  dream_idx_type i, j;
  double ramax;
  int    nv;
  double *RESTRICT rs, retfoc = 0, weight = 0;
  int err = NONE, out_err = NONE;

  h  = (double*) malloc( nt*num_elements * sizeof(double));
  hi = (double*) malloc( nt*num_elements * sizeof(double));
  rs = (double*) malloc(num_elements*sizeof(double));

  for (i=0; i< nt; i++) {
    ha[i] = 0.0;
  }

  retfoc = 0.0;
  weight = 1.0;

  // nv - number of annulus.
  nv = (num_elements+1)/2;
  //nv = (num_elements)/2; // This is wrong!!

  for (i=0; i<nt; i++) {
    for (j=0; j<num_elements; j++) {
      h[i+j*nt] = (double) 0.0;
    }
  }

  /* ----------------------------------------- */

  for (i=0; i<num_elements; i++) {
    r = gr[i];

    if (i==0 && r == (double) 0.0 ) {
      for (j=0; j<nt; j++) {
        h[j] = 0.0; // h(j,1) = 0.0;
      }
    } else {
      err = circ_annular(xo,yo,zo,r,dx,dy,dt,nt,delay,v,cp,alpha,weight,h,i,num_elements,err_level);
      if (err != NONE)
        out_err = err;
    }
  }

  /* --------------------------------------------- */

  center_pos_annular(rs, gr, num_elements, nv, &ramax);

  for (i=0; i<nv; i++) {
    focusing_annular(foc_type, focal[i], rs[i], ramax, cp, &retfoc);  // Note foc_type must be 6 here!
    if (do_apod){
      apodization_annular(apod_type, i, apod, &weight, rs[i], ramax, param);
    }
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
 * center_pos_annular
 *
 ***/

void center_pos_annular(double *RESTRICT rs, double *RESTRICT gr, int num_elements, int nv, double *RESTRICT ramax)
{
  int i, ns = 0;

  rs[0] = gr[0]/2;
  for (i=2; i<num_elements; i += 2) {
    ns = (i+1) / 2;
    rs[ns] = (gr[i]+gr[i-1]) / 2;
  }
  *ramax = rs[nv-1];

  return;
}

/***
*
* subroutine focussing gives le retard retfoc du au focussing
*
***/

void focusing_annular(int foc_type, double focal, double rs, double ramax, double cp, double *RESTRICT retfoc)
{
  double diff, rmax;

  // foc_type = 1 - no foc, 2 foc xy.
  if (foc_type == 1) {
    *retfoc = 0.0;
    return;
  }

  if (foc_type == 2) {
    rmax = sqrt(ramax*ramax + focal*focal);
    diff = rmax - sqrt(rs*rs + focal*focal);

    *retfoc = diff * 1000 / cp;
  }

  // User defined focusing.
  if (foc_type == 6)
    *retfoc = focal; // Here focal is the user defined time delay in [us] (not the focal depth)

  return;
}

/***
 *
 *  apodization
 *
 * apodization apod_type = 0 apodization with imported apodisation function apod(x,y)
 *
 * param=input parameter
 *
 * apod_type = 0 - user defined.
 * apod_type = 1 traingle.
 * apod_type = 2 gauss.
 * apod_type = 3 rised cosine
 * apod_type = 4 simply supported.
 * apod_type = 5 clamped.
 *
 ***/

void apodization_annular(int apod_type, int i, double *RESTRICT apod, double *RESTRICT weight, double rs,
                         double ramax, double param)
{
  double pi = atan((double) 1.0) * (double) 4.0;

  switch(apod_type) {

  case 0:
    *weight = apod[i];
    break;

  case 1:
    *weight = 1.0 - fabs(rs) / ramax;
    break;

  case 2:
    *weight = exp(-(param * rs*rs) / (ramax*ramax));
    break;

  case 3:
    *weight = param + cos(rs * pi / ramax);
    break;

  case 4:
    *weight = 1.0 - rs*rs / (ramax*ramax);
    break;

  case 5:
    *weight = (1.0 - rs*rs / (ramax*ramax)) * (1.0 - rs*rs / (ramax*ramax));
    break;

  default:
    break;
  }

  return;
}

/***
 *
 * subroutine circ - pour calculer pulse respone of a circular aperture
 *
 * a = R1,R2,...,Rk,...,Rnum_elements = gr(k)
 *
 ***/

int circ_annular(double xo, double  yo, double  zo, double  r, double dx, double dy, double dt,
                  dream_idx_type nt, double delay, double v, double cp, double alpha, double weight, double *RESTRICT h,
                  dream_idx_type k, int num_elements, int err_level)
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

  y = ysmin + dy/2.0;
  while (y <= ysmax) {

    xlimit_annular(y, r, &xsmin, &xsmax);

    x = xsmin + dx/2.0;
    while (x <= xsmax) {

      distance(xo, yo, zo, x, y, zs, &ri);
      ai = weight * v * ds / (2*pi*ri);
      ai /= dt;
      ai *= 1000;               // Convert to SI units.

      // Propagation delay in micro seconds.
      t = ri * 1000 / cp;
      it = (dream_idx_type) rint((t - delay)/dt);

      if ((it < nt) && (it >= 0))

        if (alpha == (double) 0.0) {
          h[it + k*nt] += ai;
        }
        else {
          att_annu(alpha, ri, it, dt, cp, h, nt, ai, k, num_elements);
        }
      else  {
        if  (it >= 0)
          err = dream_out_of_bounds_err("SIR out of bounds",it-nt+1,err_level);
        else
          err = dream_out_of_bounds_err("SIR out of bounds",it,err_level);

        if ( (err_level == PARALLEL_STOP) || (err_level == STOP) )
          return err; // Bail out.
      }

      x += dx;
    }
    y += dy;
  }

  return err;
}

/***
 *
 * xlimit_annular
 *
 * Computes the x-axis integration limits.
 *
 ***/

void xlimit_annular(double yi, double a, double *RESTRICT xsmin, double *RESTRICT xsmax)
{
  double rs;

  rs = sqrt(a*a - yi*yi);
  *xsmin = -rs;
  *xsmax = rs;

  return;
}

/***
 *
 * call superpos(h,ha) subroutine pour superposer les contributions des elements
 *
 * ha = output response
 * h  = input responce of actual element
 ***/

void superpos_annular(double *RESTRICT hi, double *RESTRICT ha, dream_idx_type  nt, double weight, double retfoc, dream_idx_type j, double  dt)
{
  double *RESTRICT buf;
  dream_idx_type    i,it1;

  buf = (double*) malloc(2*nt*sizeof(double));

  it1 = (dream_idx_type) (retfoc / dt) + 1;

  for (i=0; i<2*nt; i++)
    buf[i] = (double) 0.0;

  for (i=0; i<nt; i++)
    buf[i+it1] = hi[i+j*nt];

  for (i=0; i<nt; i++)
    ha[i] += weight * buf[i];

  free(buf);

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
