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

// $Revision: 775 $ $Date: 2014-05-16 14:50:24 +0200 (Fri, 16 May 2014) $ $LastChangedBy: frli8848 $

#include <math.h>
#include <stdio.h>
#include "das.h"
#include "dream_error.h"

#if defined(_MSC_VER) || defined(__LCC__)
#include "msvc_rint.h"
#endif

/***
 *
 *  Delay-and-sum response.
 *
 ***/

/***
 *
 *  das - Delay-and-sum.
 *
 ***/

int das(double xo, double yo, double zo,
	double dt, dream_idx_type  nt, double delay,
	double cp, double *RESTRICT h, int err_level)
{
  dream_idx_type i,it;
  double t,tt;
  double ri;
  double qan;
  int err = NONE;

  for (i = 0; i < nt; i++) {
    h[i] = (double) 0.0 ;
  }
  
  ri = sqrt( xo*xo + yo*yo + zo*zo );

  t = ri * 1000.0/cp;
  tt = t - delay;
  qan = tt / dt;
  it = 0;
  it = (dream_idx_type) rint(qan);

  // Check if index is out of bounds.
  if ((it < nt) && (it >= 0)) {
    h[it] += 1.0;
  }
  else  {
    if  (it >= 0)
      err = dream_out_of_bounds_err("Delay out of bounds",it-nt+1,err_level);
    else
      err = dream_out_of_bounds_err("Delay out of bounds",it,err_level);

    if ( (err_level == PARALLEL_STOP) || (err_level == STOP) )    
      return err; // Bail out.
  }

  return err;
} /* das */

