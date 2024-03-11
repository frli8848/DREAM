/***
*
* Copyright (C) 2024 Fredrik Lingvall
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

#include <thread>

#include "dream.h"
#include "dream_error.h"

class FFTConvP
{
 public:

 FFTConvP()
   : m_out_err(ErrorLevel::none)
    {;}

  ~FFTConvP()  = default;

  ErrorLevel run(double *Y,
                 const double *A, dream_idx_type A_M, dream_idx_type A_N,
                 const double *B, dream_idx_type B_M, dream_idx_type B_N,
                 ConvMode mode=ConvMode::equ);

 static void abort(int signum);
 bool is_running();

 private:

 void* smp_dream_fftconv(void *arg);
 std::thread conv_thread(void *arg) {
   return std::thread(&FFTConvP::smp_dream_fftconv, this, arg);
 }

 ErrorLevel m_out_err;
};
