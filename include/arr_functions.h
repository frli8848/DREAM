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

#pragma once

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

#define APOD_UD 0               // User defined
#define APOD_TRIANGLE 1
#define APOD_GAUSS 2
#define APOD_RISED_COSINE 3
#define APOD_SIMPLY_SUPPORTED 4
#define APOD_CLAMPED 5

/***
 *
 * Header file for functions common to the array functions.
 *
 ***/

void center_pos(double *xs, double *ys, double *zs, int i, double *gx, double *gy, double *gz);
void max_dim_arr(double *xamax, double *y_max, double *r_max, double *gx, double *gy, double *gz, int num_elements);
void focusing(int foc_type, double focal, double xs, double ys,
              double x_max, double y_max, double r_max, double cp, double *retfoc);
void beamsteering(int steer_type, double theta, double phi, double xs, double ys,
                  double x_max, double y_max, double r_max, double cp, double *retsteer);
void apodization(int apod_type, int i, double  *apod_vec, double *weight,
                 double xs, double ys, double r_max, double param);
void distance(double xo, double yo, double zo,double xs,double ys, double zs, double *r);
void superpos(double *h, double *ha, dream_idx_type nt);
