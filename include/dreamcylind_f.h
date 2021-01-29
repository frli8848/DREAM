/***
*
* Copyright (C) 2002,2003,2005,2006,2007,2008,2009,2021 Fredrik Lingvall
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

/***
 *
 * dreamcylind_f
 *
 ***/

int dreamcylind_f(double xo,
                  double yo,
                  double zo,
                  double a,
                  double b,
                  double R,
                  double dx,
                  double dy,
                  double dt,
                  dream_idx_type nt,
                  double delay,
                  double v,
                  double cp,
                  double alpha,
                  double *h,
                  int err_level);
