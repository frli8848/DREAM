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

#include "py_fftconv_p.hpp"

PYBIND11_MODULE(fftconv_p, m)
{
  m.doc() = "Y = fftconv_p(A,B)\n\
\n\
fftconv_p is a part of the DREAM Toolbox available at\n\
https://github.com/frli8848/DREAM.\n\
\n\
Copyright 2006-2024 Fredrik Lingvall.";

  m.def("fftconv_p", &py_fftconv_p<double>,
        "Y = fftconv_p(A,B)",
        py::arg("A"),
        py::arg("B"));
}
