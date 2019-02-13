/***
*
* Copyright (C) 2002,2003,2004,2006,2007,2008,2009,2012,2014 Fredrik Lingvall
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

// $Revision: 774 $ $Date: 2014-05-16 09:40:39 +0200 (Fri, 16 May 2014) $ $LastChangedBy: frli8848 $

#ifdef USE_FFTW

#ifdef __cplusplus 
#include <complex> // Must be included before fftw3.h
#else
#include <complex.h> // Must be included before fftw3.h
#endif

#include <fftw3.h>
#endif

#include "dream.h"


#ifdef USE_FFTW

#ifdef __cplusplus
extern "C" 
#endif
void fft_init(dream_idx_type n, fftw_complex *xc, double *y);

#ifdef __cplusplus 
extern "C" 
#endif
void fft_close();

#endif

#ifdef __cplusplus 
extern "C" 
#endif
void dream_fft(double *x, double *y, dream_idx_type n);
#ifdef __cplusplus 
extern "C" 
#endif
void dream_ifft(double *x, double *y, dream_idx_type n);


#ifdef USE_FFTW

#ifdef __cplusplus 
extern "C" 
#endif
void cr_ifft(fftw_complex *xc, double *y, dream_idx_type n);

#ifdef __cplusplus 
extern "C" 
#endif
void cc_ifft(fftw_complex *xc, fftw_complex *yc, dream_idx_type n);

#else

#ifdef __cplusplus 
extern "C" 
#endif
void cr_ifft(double *xir, double *xii, double *y, dream_idx_type n);

#ifdef __cplusplus 
extern "C" 
#endif
void cc_ifft(double *xir, double *xii, double *yor, double *yoi, dream_idx_type n);

#endif
