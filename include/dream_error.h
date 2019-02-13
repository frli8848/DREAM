/***
*
* Copyright (C) 2003,2004,2005,2006,2008,2009,2012 Fredrik Lingvall
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

// $Revision: 730 $ $Date: 2012-01-25 13:00:34 +0100 (Wed, 25 Jan 2012) $ $LastChangedBy: frli8848 $

// Error levels.
#define NONE   0
#define STOP   -1
#define WARN   -2
#define IGNORE -3
#define PARALLEL_STOP -4

#define TRUE 1
#define FALSE 0

#ifdef __cplusplus 
extern "C" 
#endif
int dream_out_of_bounds_err(const char *msg, int idx, int err_level);

#ifdef __cplusplus 
extern "C" 
#endif
void dream_err_msg(const char *msg);
