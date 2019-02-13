/***
*
* Copyright (C) 2007-2009,2012,2014-x2016,2019 Fredrik Lingvall
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

// $Revision: 910 $ $Date: 2019-02-06 07:51:39 +0100 (Wed, 06 Feb 2019) $ $LastChangedBy: frli8848 $

#ifndef __DREAM__
#define __DREAM__

// Check if we support C99.
#if __STDC_VERSION__ >= 199901L
#define RESTRICT restrict
//#define RESTRICT
#else
#define RESTRICT
#endif

#define DOUBLE_PRECISION 0
#define SINGLE_PRECISION 1

#ifdef OCTAVE

// Use the index type that Octave was build with.
#include <octave-config.h>
typedef octave_idx_type dream_idx_type;

#else

#include <stdlib.h>
// Matlab use this
typedef size_t dream_idx_type;

#endif

#endif // __DREAM__
