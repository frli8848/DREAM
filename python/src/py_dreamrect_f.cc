/***
*
* Copyright (C) 2023,2024 Fredrik Lingvall
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

#include "dreamrect_f.h"
#include "arg_parser_python.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

py::array_t<double,py::array::f_style> py_dreamrect_f(py::array_t<double,py::array::f_style> *py_ro,
                                                      py::array_t<double,py::array::f_style> *py_geom_par,
                                                      py::array_t<double,py::array::f_style> *py_s_par,
                                                      py::array_t<double,py::array::f_style> *py_delay,
                                                      py::array_t<double,py::array::f_style> *py_m_par,
                                                      std::string foc_str, py::array_t<double,py::array::f_style> *py_focal,
                                                      std::string err_level_str);

py::array_t<double,py::array::f_style> py_dreamrect_f(py::array_t<double,py::array::f_style> *py_ro,
                                                      py::array_t<double,py::array::f_style> *py_geom_par,
                                                      py::array_t<double,py::array::f_style> *py_s_par,
                                                      py::array_t<double,py::array::f_style> *py_delay,
                                                      py::array_t<double,py::array::f_style> *py_m_par,
                                                      std::string foc_str, py::array_t<double,py::array::f_style> *py_focal,
                                                      std::string err_level_str)
{
  ArgParser<double> ap;

  //
  // Observation points.
  //

  if (!ap.check_obs_points("dreamrect_f", py_ro)) {
    throw std::runtime_error("Error in dreamrect_f!");
  }

  py::buffer_info ro_info = py_ro->request();
  auto ro_m = ro_info.shape[0];
  double *ro = static_cast<double*>(ro_info.ptr);

  dream_idx_type no = (dream_idx_type) ro_m;

  //
  // Transducer geometry
  //

  double a=0.0, b=0.0, dummy=0.0;
  if (!ap.parse_geometry("dreamrect_f", py_geom_par, 2, a, b, dummy)) {
    throw std::runtime_error("Error in dreamrect_f!");
  }

  //
  // Temporal and spatial sampling parameters.
  //

  double dx=0.0, dy=0.0, dt=0.0;
  dream_idx_type nt=0;
  if (!ap.parse_sampling("dreamrect_f", py_s_par, 4, dx, dy, dt, nt)) {
    throw std::runtime_error("Error in dreamrect_f!");
  }

  //
  // Start point of impulse response vector ([us]).
  //

  if (!ap.check_delay("dreamrect_f", py_delay, no)) {
    throw std::runtime_error("Error in dreamrect_f!");
  }

  py::buffer_info delay_info = py_delay->request();
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
  if (!ap.parse_material("dreamrect_f", py_m_par, v, cp, alpha)) {
    throw std::runtime_error("Error in dreamrect_f!");
  }

  //
  // Focusing parameters.
  //

  FocusMet foc_met = FocusMet::none;

  double focal=0.0;
  if (!ap.parse_focus_args("dreamrect_f", foc_str, py_focal, 1, foc_met, &focal)) {
    throw std::runtime_error("Error in dreamrect_f!");
  }

  //
  // Error reporting.
  //

  ErrorLevel err_level=ErrorLevel::stop;

  if (!ap.parse_error_arg("dreamrect_f", err_level_str, err_level)) {
    throw std::runtime_error("Error in dreamrect_f!");
  }

  // Create an output matrix for the impulse response
  auto py_h_mat = py::array_t<double,py::array::f_style>({nt,no});
  py::buffer_info h_mat_info = py_h_mat.request();
  double *h = static_cast<double*>(h_mat_info.ptr);

  SIRData hsir(h, nt, no);
  hsir.clear();

  Rect_f rect_f;

  // Register signal handler.
  std::signal(SIGINT, Rect_f::abort);

  //
  // Call DREAM function
  //

  SIRError err = rect_f.dreamrect_f(alpha,
                                    ro, no,
                                    a, b,
                                    foc_met, focal,
                                    dx, dy,  dt, nt,
                                    delay_type, delay,
                                    v, cp,
                                    h, err_level);

  if (!rect_f.is_running()) {
    if (err == SIRError::out_of_bounds) {
      dream_err_msg("SIR out-of-bounds!\n"); // Bail out.
    } else {
      dream_err_msg("CTRL-C pressed!\n"); // Bail out.
    }

    throw std::runtime_error("Error in dreamrect_f!");
  }

  if (err == SIRError::warn_out_of_bounds) {
    std::cout << "Warning: SIR out-of-bounds!" << std::endl;
  }

  return py_h_mat;
}

/***
 *
 *  Python gateway function for (parallel) dreamrect_f.
 *
 ***/

PYBIND11_MODULE(dreamrect_f, m)
{
  m.doc() = "H = dreamrect_f(Ro,geom_par,s_par,delay,m_par,foc_met,focal,err_level)\n\
\n\
DREAMRECT_F - Computes the spatial impulse response\n\
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
Focusing parameters: foc_met and focal:\n\
\n\
'foc_met'\n\
Focusing method, options are: \"off\", \"x\", \"y\", \"xy\", and \"x+y\".\n\
'focal'\n\
   Focal distance [mm].\n\
\n\
Error Handling: err_level;\n\
'err_level' is an optional text string parameter for controlling the error behavior, options are:\n\
\n\
\"ignore\"\n\
   An error is ignored (no error message is printed and the program is not stopped) but the err output \n\
   argument is non-zero if an error occured.\n\
\"warn\"\n\
   A warning message is printed but the program in not stopped (and err argument is non-zero).\n\
\"stop\"\n\
  An error message is printed and the program is stopped.\n\
\n\
dreamrect_f is a part of the DREAM Toolbox available at\n\
https://github.com/frli8848/DREAM.\n\
\n\
Copyright 2006-2023 Fredrik Lingvall.";

  m.def("dreamrect_f", &py_dreamrect_f,
        "H = dreamrect_f(Ro,geom_par,G,s_par,delay,m_par,foc_met,focal,err_level)",
        py::arg("Ro"),
        py::arg("geom_par"),
        py::arg("s_par"),
        py::arg("delay"),
        py::arg("m_par"),
        py::arg("foc_met"),
        py::arg("focal"),
        py::arg("err_level"));
}
