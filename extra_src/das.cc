/***
*
* Copyright (C) 2003,2006,2007,2008,2009,2014,2021 Fredrik Lingvall
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

#include <cmath>

#include "das.h"
#include "dream_error.h"

/***
 *
 *  das - Delay-and-sum.
 *
 ***/

ErrorLevel das(double xo, double yo, double zo,
               double dt, dream_idx_type  nt, double delay,
               double cp, double *h, ErrorLevel err_level)
{
  ErrorLevel err = ErrorLevel::none;

  for (dream_idx_type it = 0; it<nt; it++) {
    h[it] = 0.0 ;
  }

  double ri = std::sqrt(xo*xo + yo*yo + zo*zo);
  double t = ri * 1000.0/cp;
  dream_idx_type it = (dream_idx_type) rint((t - delay)/dt);

  // Check if index is out of bounds.
  if ((it < nt) && (it >= 0)) {
    h[it] += 1.0;
  } else  {
    if  (it >= 0) {
      err = dream_out_of_bounds_err("Delay out of bounds",it-nt+1,err_level);
    } else {
      err = dream_out_of_bounds_err("Delay out of bounds",it,err_level);
    }

    if ( (err_level == ErrorLevel::parallel_stop) || (err_level == ErrorLevel::stop) ) {
      return err; // Bail out.
    }
  }

  return err;
}
