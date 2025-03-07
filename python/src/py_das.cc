/***
*
* Copyright (C) 2024,2025 Fredrik Lingvall
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

#include "py_das.hpp"

#ifdef USE_FLOAT

PYBIND11_MODULE(das_f, m)
{
  m.doc() = "H = das_f(Y,Gt,Gr,Ro,dt,delay,cp,method,err_level,device,verbose)\n\
\n\
das is a part of the DREAM Toolbox available at\n\
https://github.com/frli8848/DREAM.\n\
\n\
Copyright 2006-2025 Fredrik Lingvall.";

  m.def("das", &py_das<float>,
        "Im = das(Y,Gt,Gr,Ro,dt,delay,cp,method,err_level,device,verbose)",
        py::arg("Y"),
        py::arg("Gt"),
        py::arg("Gr"),
        py::arg("Ro"),
        py::arg("dt"),
        py::arg("delay"),
        py::arg("cp"),
        py::arg("method"),
        py::arg("err_level"),
        py::arg("device"),
        py::arg("verbose"));
}

#else

PYBIND11_MODULE(das, m)
{
  m.doc() = "H = das(Y,Gt,Gr,Ro,dt,delay,cp,method,err_level,device,verbose)\n\
\n\
das is a part of the DREAM Toolbox available at\n\
https://github.com/frli8848/DREAM.\n\
\n\
Copyright 2006-2025 Fredrik Lingvall.";

  m.def("das", &py_das<double>,
        "Im = das(Y,Gt,Gr,Ro,dt,delay,cp,method,err_level,device,verbose)",
        py::arg("Y"),
        py::arg("Gt"),
        py::arg("Gr"),
        py::arg("Ro"),
        py::arg("dt"),
        py::arg("delay"),
        py::arg("cp"),
        py::arg("method"),
        py::arg("err_level"),
        py::arg("device"),
        py::arg("vebose"));
}
#endif
