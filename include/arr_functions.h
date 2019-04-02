/***
*
* Copyright (C) 2002,2003,2006,2007,2008,2009,2019 Fredrik Lingvall
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
* Free
* Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
* 02110-1301, USA.
*
***/

#include "dream.h"

#define NO_FOCUS 1
#define FOCUS_X 2
#define FOCUS_Y 3
#define FOCUS_XY 4
#define FOCUS_X_Y 5             // X+Y
#define FOCUS_UD 6              // User defined

#define NO_STEER 1
#define STEER_X 2
#define STEER_Y 3
#define STEER_XY 4

#define IPOD_UD 0               // User defined
#define IPOD_TRIANGLE  1
#define IPOD_GAUSS  2
#define IPOD_RISED_COSINE 3
#define IPOD_SIMPLY_SUPPORTED 4
#define IPOD_CLAMPED 5

/***
 *
 * Header file for functions common to the array functions.
 *
 ***/

#ifdef __cplusplus
extern "C"
#endif
void center_pos(double *xs, double *ys, double *zs, int i, double *gx, double *gy, double *gz);

#ifdef __cplusplus
extern "C"
#endif
void max_dim_arr(double *xamax, double *yamax, double *ramax, double *gx, double *gy, double *gz, int isize);

#ifdef __cplusplus
extern "C"
#endif
void focusing(int ifoc, double focal, double xs, double ys,
              double xamax, double yamax, double ramax, double cp, double *retfoc);

#ifdef __cplusplus
extern "C"
#endif
void beamsteering(int ister, double theta, double phi, double xs, double ys,
                  double xamax, double yamax, double ramax, double cp, double *retsteer);

#ifdef __cplusplus
extern "C"
#endif
void apodization(int iweight, int iapo, int i, double  *apod, double *weight,
               double xs, double ys, double ramax, double param, int isize);

#ifdef __cplusplus
extern "C"
#endif
void distance(double xo, double yo, double zo,double xs,double ys, double zs,double *ri);

#ifdef __cplusplus
extern "C"
#endif
void superpos(double *h, double *ha, dream_idx_type nt);

#ifdef __cplusplus
extern "C"
#endif
void check_delay(dream_idx_type it, double tt, int *icheck, dream_idx_type nt);
