/***
*
* Copyright (C) 2003,2004,2006,2008,2009,2012,2016,2021,2022,2023,2024 Fredrik Lingvall
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

#include "dream_error.h"

//#if defined(DREAM_OCTAVE) || defined(DREAM_PYTHON)
#include <iostream>
//#endif

#ifdef DREAM_MATLAB
#include "mex.h"
#endif

/***
 *
 * Functions for error handling.
 *
 ***/

SIRError dream_out_of_bounds_err(const char *msg, dream_idx_type idx, ErrorLevel err_level, bool verbose)
{
  SIRError err = SIRError::none;

  switch(err_level) {

  case ErrorLevel::ignore:
    err = SIRError::ignore_out_of_bounds;
    break; // No printouts.

  case ErrorLevel::warn:
    {
      err = SIRError::warn_out_of_bounds;

      if (verbose) {
        std::cout << "Warning: " << msg << " (offset = " << idx << " samples)" << std::endl;
      }
    }
    break;

  case ErrorLevel::stop:
  case ErrorLevel::parallel_stop:
    {
      err = SIRError::out_of_bounds;

      if (verbose) {
        std::cout << "Error: " << msg << " (offset = " << idx << " samples)" << std::endl;
      }
    }
    break;

  default:
    break;
  }

  return err;
}

void dream_err_msg(const char *msg)
{
#if defined(DREAM_OCTAVE) || defined(DREAM_PYTHON) || defined(DREAM_JULIA)
  std::cout << msg << std::endl;
#endif

#ifdef DREAM_MATLAB
  mexErrMsgTxt(msg);
#endif

  return;
}
