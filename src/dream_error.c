/***
*
* Copyright (C) 2003,2004,2006,2008,2009,2012,2016 Fredrik Lingvall
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

// $Revision: 890 $ $Date: 2016-09-27 16:56:06 +0200 (Tue, 27 Sep 2016) $ $LastChangedBy: frli8848 $

#include "dream_error.h"

#ifdef OCTAVE
#include <stdio.h>
#else
#include <uchar.h>
#include <mex.h>
#endif

/***
 *
 * Functions for error handling.
 *
 ***/

int dream_out_of_bounds_err(const char *msg, int idx, int err_level)
{

  switch(err_level) {

  case STOP:

#ifdef OCTAVE
    printf("%s (offset = %d samples)\n",msg,idx);
#else
    mexPrintf("%s (offset = %d samples)\n",msg,idx);
    mexErrMsgTxt(""); // Bail out!
#endif
    break;

  case WARN:
#ifdef OCTAVE
    printf("%s (offset = %d samples)\n",msg,idx);
#else
    mexPrintf("%s (offset = %d samples)\n",msg,idx);
#endif
    break;

  case IGNORE:
    break; // Do nothing.

  case PARALLEL_STOP:
#ifdef OCTAVE
    printf("%s (offset = %d samples)\n",msg,idx);
#else
    mexPrintf("%s (offset = %d samples)\n",msg,idx);
#endif
    break;

  default:
    break;
  }

  return err_level;
}

void dream_err_msg(const char *msg)
{

#ifdef OCTAVE
  printf("%s\n",msg);
#else
  mexErrMsgTxt(msg);
#endif

  return;
}
