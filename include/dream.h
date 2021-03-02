/***
*
* Copyright (C) 2007-2009,2012,2014-x2016,2019,2021 Fredrik Lingvall
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


#ifndef __DREAM__
#define __DREAM__

#define SINGLE_DELAY 0
#define MULTIPLE_DELAYS 1

enum class DelayType {
  single,
    multiple};

#define NO_FOCUS 1
#define FOCUS_X 2
#define FOCUS_Y 3
#define FOCUS_XY 4
#define FOCUS_X_Y 5             // X+Y
#define FOCUS_UD 6              // User defined

enum class FocusMet {
  none,
    x,
    y,
    xy,
    x_y,
    ud};

#define NO_STEER 1
#define STEER_X 2
#define STEER_Y 3
#define STEER_XY 4

enum class SteerMet {
  none,
    x,
    y,
    xy};

#define APOD_UD 0               // User defined
#define APOD_TRIANGLE 1
#define APOD_GAUSS 2
#define APOD_RISED_COSINE 3
#define APOD_SIMPLY_SUPPORTED 4
#define APOD_CLAMPED 5

enum class ApodMet {
  triangle,
    gauss,
    raised_cosine,
    simply_supported,
    clamped,
    ud};

#ifdef DREAM_OCTAVE

// Use the index type that Octave was build with.
#include <octave-config.h>
typedef octave_idx_type dream_idx_type;

// Macros
#define mxGetM(N)   args(N).matrix_value().rows()
#define mxGetN(N)   args(N).matrix_value().cols()
#define mxIsChar(N) args(N).is_string()

#else

#include <stdlib.h>
// Matlab use this
typedef size_t dream_idx_type;

#endif

#endif // __DREAM__
