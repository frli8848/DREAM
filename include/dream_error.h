/***
 *
 * Copyright (C) 2003,2004,2005,2006,2008,2009,2012,2021,2023,2024 Fredrik Lingvall
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

SIRError dream_out_of_bounds_err(const char *msg, dream_idx_type idx, ErrorLevel err_level, bool verbose=false);
void dream_err_msg(const char *msg);
