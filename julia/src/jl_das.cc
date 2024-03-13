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

#include <csignal>
#include <memory>

#include "das.h"

// Persistent smart pointer to DAS object.
#ifdef USE_FLOAT
std::unique_ptr<DAS<float>> das=nullptr;
#else
std::unique_ptr<DAS<double>> das=nullptr;
#endif

/***
 *
 *  Julia gateway function for (parallel) DAS.
 *
 ***/

#include "jl_das.hpp"

#ifdef USE_FLOAT

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.method("das", jl_das<float>);
}

#else

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.method("das", jl_das<double>);
}

#endif
