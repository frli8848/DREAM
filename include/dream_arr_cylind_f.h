/***
*
* Copyright (C) 2002,2003,2005,2006,2007,2008,2009 Fredrik Lingvall
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

// $Revision: 565 $ $Date: 2009-09-17 22:24:06 +0200 (Thu, 17 Sep 2009) $ $LastChangedBy: dream $

#include "dream.h"

/***
 *
 * Header file for dream_arr_cylind_f.c
 *
 ***/

#ifdef __cplusplus 
extern "C" 
#endif
int dream_arr_cylind_f(double xo, double yo, double zo, double a, double b, double R, double dx, 
		       double dy, double dt,dream_idx_type nt, double delay, double v, double cp, double alfa, 
		       int isize,double *gx, double *gy, double *gz, int ifoc, double focal, int ister, 
		       double theta, double phi, double *apod, int iweight, int iapo, double param, 
		       double *ha, int err_level);

#ifdef __cplusplus 
extern "C" 
#endif
int dream_arr_cylind_udf(double xo, double yo, double zo, double a, double b, double R, 
			 double dx, double dy, double dt,
			 dream_idx_type nt, double delay, double v, double cp, double alfa, int isize, 
			 double *gx, double *gy, double *gz, int ifoc, double *focal,
			 int ister, double theta, double phi, 
			 double *apod, int iweight, int iapo, double param, double *ha, int err_level);
