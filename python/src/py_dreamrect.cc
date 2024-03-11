/***
*
* Copyright (C) 2023 Fredrik Lingvall
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

#include "dreamrect.h"
#include "arg_parser_python.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

py::array_t<double,py::array::f_style> py_dreamrect(py::array_t<double,py::array::f_style> *py_ro,
                                                    py::array_t<double,py::array::f_style> *py_geom_par,
                                                    py::array_t<double,py::array::f_style> *py_s_par,
                                                    py::array_t<double,py::array::f_style> *py_delay,
                                                    py::array_t<double,py::array::f_style> *py_m_par,
                                                    std::string err_level_str);

py::array_t<double,py::array::f_style> py_dreamrect(py::array_t<double,py::array::f_style> *py_ro,
                                                    py::array_t<double,py::array::f_style> *py_geom_par,
                                                    py::array_t<double,py::array::f_style> *py_s_par,
                                                    py::array_t<double,py::array::f_style> *py_delay,
                                                    py::array_t<double,py::array::f_style> *py_m_par,
                                                    std::string err_level_str)
{
  ArgParser<double> ap;

  //
  // Observation points.
  //

  if (!ap.check_obs_points("dreamrect", py_ro)) {
    throw std::runtime_error("Error in dreamrect!");
  }

  py::buffer_info ro_info = py_ro->request();
  auto ro_ndim = ro_info.ndim;
  auto ro_m = ro_info.shape[0];
  auto ro_n = ro_info.shape[1];
  double *ro = static_cast<double*>(ro_info.ptr);

  dream_idx_type no = (dream_idx_type) ro_m;

  //
  // Transducer geometry
  //

  double a=0.0, b=0.0, dummy=0.0;
  if (!ap.parse_geometry("dreamrect", py_geom_par, 2, a, b, dummy)) {
    throw std::runtime_error("Error in dreamrect!");
  }

  //
  // Temporal and spatial sampling parameters.
  //

  double dx=0.0, dy=0.0, dt=0.0;
  dream_idx_type nt=0;
  if (!ap.parse_sampling("dreamrect", py_s_par, 4, dx, dy, dt, nt)) {
    throw std::runtime_error("Error in dreamrect!");
  }

  //
  // Start point of impulse response vector ([us]).
  //

  if (!ap.check_delay("dreamrect", py_delay, no)) {
    throw std::runtime_error("Error in dreamrect!");
  }

  py::buffer_info delay_info = py_delay->request();
  auto delay_ndim = delay_info.ndim;
  auto delay_m = delay_info.shape[0];
  auto delay_n = delay_info.shape[1];
  double *delay = static_cast<double*>(delay_info.ptr);

  DelayType delay_type = DelayType::single; // delay is a scalar.
  if (delay_m * delay_n != 1) {
    delay_type = DelayType::multiple; // delay is a vector.
  }

  //
  // Material parameters
  //

  double v=1.0, cp=1000.0, alpha=0.0;
  if (!ap.parse_material("dreamrect", py_m_par, v, cp, alpha)) {
    throw std::runtime_error("Error in dreamrect!");
  }

  //
  // Error reporting.
  //

  ErrorLevel err=ErrorLevel::none, err_level=ErrorLevel::stop;

  if (!ap.parse_error_arg("dreamrect", err_level_str, err_level)) {
    throw std::runtime_error("Error in dreamrect!");
  }

  // Create an output matrix for the impulse response
  auto py_h_mat = py::array_t<double,py::array::f_style>({nt,no});
  py::buffer_info h_mat_info = py_h_mat.request();
  double *h = static_cast<double*>(h_mat_info.ptr);

  SIRData hsir(h, nt, no);
  hsir.clear();

  Rect rect;

  // Register signal handler.
  std::signal(SIGINT, Rect::abort);

  //
  // Call DREAM function
  //

  err = rect.dreamrect(alpha,
                       ro, no,
                       a, b,
                       dx, dy,  dt, nt,
                       delay_type, delay,
                       v, cp,
                       h, err_level);

  if (!rect.is_running()) {
    dream_err_msg("CTRL-C pressed!\n"); // Bail out.
    throw std::runtime_error("Error in dreamrect!");
  }

  if (err == ErrorLevel::stop) {
    throw std::runtime_error("Error in dreamrect!");
  }

  return py_h_mat;
}

/***
 *
 *  Python gateway function for (parallel) dreamrect.
 *
 ***/

PYBIND11_MODULE(dreamrect, m)
{
  m.doc() = "H = dreamrect(Ro,geom_par,s_par,delay,m_par,err_level)\n\
\n\
DREAMRECT - Computes the spatial impulse response\n\
for a rectangular transducer using parallel processing\n\
(using threads).\n\
\n\
Observation point(s) ([mm]):\n\
\n\
'Ro'\n\
   An N x 3 matrix, Ro = np.asmatrix([xo1 yo1 zo2; xo2 yo2 zo2; ... xoN yoN zoN]) where N is the number of observation points.\n\
\n\
Geometrical parameters: geom_par = np.asmatrix([a,b]);\n\
\n\
'a'\n\
   x-direction size [mm].\n\
'b'\n\
   y-direction size [mm].\n\
\n\
Sampling parameters: s_par = np.asmatrix([dx,dy,dt,nt]) \n\
\n\
'dx'\n\
   Spatial x-direction discretization size [mm].\n\
'dy'\n\
   Spatial y-direction discretization size [mm].\n\
'dt'\n\
   Temporal discretization period (= 1/sampling freq) [us].\n\
'nt'\n\
   Length of impulse response vector.\n\
\n\
Start point of SIR:\n\
\n\
'delay'\n\
   Scalar delay for all observation points or a vector with individual delays for each observation point [us].\n \
\n\
Material parameters: m_par = np.asmatrix([v,cp,alpha])\n\
\n\
'v'\n\
   Normal velocity [m/s].\n\
'cp'\n\
   Sound velocity [m/s].\n\
'alpha'\n\
   Attenuation coefficient [dB/(cm MHz)].\n\
\n\
Error Handling: err_level;\n\
'err_level' is an optional text string parameter for controlling the error behavior, options are:\n\
\n\
\"ignore\"\n\
   An error is ignored (no error message is printed and the program is not stopped) but the err output \n\
   argument is negative if an error occured.\n\
\"warn\"\n\
   A warning message is printed but the program in not stopped.\n\
\"stop\"\n\
  An error message is printed and the program is stopped.\n\
\n\
dreamrect is a part of the DREAM Toolbox available at\n\
https://github.com/frli8848/DREAM.\n\
\n\
Copyright 2006-2023 Fredrik Lingvall.";

  m.def("dreamrect", &py_dreamrect,
        "H = dreamrect(Ro,geom_par,G,s_par,delay,m_par,err_level)",
        py::arg("Ro"),
        py::arg("geom_par"),
        py::arg("s_par"),
        py::arg("delay"),
        py::arg("m_par"),
        py::arg("err_level"));
}
