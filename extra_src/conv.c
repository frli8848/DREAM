/***
*
* Copyright (C) 2009,2015 Fredrik Lingvall
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

// $Revision: 879 $ $Date: 2015-05-08 15:05:55 +0200 (Fri, 08 May 2015) $ $LastChangedBy: frli8848 $

#include "dream.h"

#define TRUE 1
#define FALSE 0

#define EQU 0
#define SUM 1
#define NEG 2

/***
 *
 * Convolution of two vectors.
 *
 ***/

void conv(double *xr, dream_idx_type nx, double *yr, dream_idx_type ny, double *zr, int in_place, int mode)
{
  dream_idx_type i,j;

  // Don't clear if in-place '+=' and '-=' modes.
  if (in_place == FALSE || mode == EQU) {
    for (i=0; i < (nx+ny)-1; i++) {
      zr[i]=0.0;
    }
  }

  if (in_place == FALSE) { // Normal mode.

    for(i=0; i<nx; i++) {
      for(j=0; j<ny; j++) {
        zr[i+j] += xr[i] * yr[j];
      }
    }
  }
  else { // in-place

    switch (mode) {

    case EQU:
    case SUM:
      // in-place '=' and '+=' operation.
      for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
          zr[i+j] += xr[i] * yr[j];
        }
      }
      break;

    case NEG:
      // in-place '-=' operation.
      for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
          zr[i+j] -= xr[i] * yr[j];
        }
      }
      break;

    default:
      // in-place '=' operation.
      for(i=0; i<nx; i++) {
        for(j=0; j<ny; j++) {
          zr[i+j] += xr[i] * yr[j];
        }
      }
      break;
    } // switch
  } // in-place
}
