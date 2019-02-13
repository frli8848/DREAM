/***
*
* Copyright (C) 2008,2009 Fredrik Lingvall
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


#ifndef __RINT__
#define __RINT__

//
// MSVC don't have the rint function.
//

// From:
//   http://www.eecg.utoronto.ca/~aamodt/sourceware/MSVC.html

// Copyright (C) 2001 Tor M. Aamodt, University of Toronto
// Permisssion to use for all purposes commercial and otherwise granted.
// THIS MATERIAL IS PROVIDED "AS IS" WITHOUT WARRANTY, OR ANY CONDITION OR
// OTHER TERM OF ANY KIND INCLUDING, WITHOUT LIMITATION, ANY WARRANTY
// OF MERCHANTABILITY, SATISFACTORY QUALITY, OR FITNESS FOR A PARTICULAR
// PURPOSE.

double rint( double x)
{
  if( x > 0 ) {
    __int64 xint = (__int64) (x+0.5);
    if( xint % 2 ) {
      // then we might have an even number...
      double diff = x - (double)xint;
      if( diff == -0.5 )
        return (double) (xint-1);
    }
    return (double) (xint);
  } else {
    __int64 xint = (__int64) (x-0.5);
    if( xint % 2 ) {
      // then we might have an even number...
      double diff = x - (double) xint;
      if( diff == 0.5 )
        return (double) (xint+1);
    }
    return (double) (xint);
  }
};
#endif
