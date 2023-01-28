/***
*
* Copyright (C) 2021,2022,2023 Fredrik Lingvall
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

#pragma once

#include "dream.h"
#include "fft.h"

void fftconv(FFT &fft,
             double *xr, dream_idx_type nx, double *yr, dream_idx_type ny, double *zr,
             double *a,  double *b, double *c,
             std::complex<double> *af, std::complex<double> *bf, std::complex<double> *cf,
             ConvMode conv_mode);

void add_fftconv(FFT &fft,
                 double **H, dream_idx_type L, dream_idx_type h_len, // 3D impulse response matrix.
                 dream_idx_type n, // observation point index.
                 double *U, dream_idx_type u_len, // Input signal matrix.
                 double *z,                       // Output vector.
                 double *a, double *b, double *c,
                 std::complex<double> *af, std::complex<double> *bf, std::complex<double> *cf,
                 ConvMode conv_mode);
