/***
*
* Copyright (C) 2021,2023,2024 Fredrik Lingvall
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

#include "dream_arr_annu.h"
#include "arg_parser_python.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

py::array_t<double,py::array::f_style> py_dream_arr_annu(py::array_t<double,py::array::f_style> *py_ro,
                                                         py::array_t<double,py::array::f_style> *py_G,
                                                         py::array_t<double,py::array::f_style> *py_s_par,
                                                         py::array_t<double,py::array::f_style> *py_delay,
                                                         py::array_t<double,py::array::f_style> *py_m_par,
                                                         std::string foc_str, py::array_t<double,py::array::f_style> *py_focal,
                                                         std::string apod_str, py::array_t<double,py::array::f_style> *py_apod, double win_par,
                                                         std::string err_level_str);

py::array_t<double,py::array::f_style> py_dream_arr_annu(py::array_t<double,py::array::f_style> *py_ro,
                                                         py::array_t<double,py::array::f_style> *py_G,
                                                         py::array_t<double,py::array::f_style> *py_s_par,
                                                         py::array_t<double,py::array::f_style> *py_delay,
                                                         py::array_t<double,py::array::f_style> *py_m_par,
                                                         std::string foc_str, py::array_t<double,py::array::f_style> *py_focal,
                                                         std::string apod_str, py::array_t<double,py::array::f_style> *py_apod, double win_par,
                                                         std::string err_level_str)
{
  ArgParser<double> ap;

  //
  // Observation points.
  //

  if (!ap.check_obs_points("dream_arr_annu", py_ro)) {
    throw std::runtime_error("Error in dream_arr_annu!");
  }

  py::buffer_info ro_info = py_ro->request();
  auto ro_m = ro_info.shape[0];
  double *ro = static_cast<double*>(ro_info.ptr);

  dream_idx_type no = (dream_idx_type) ro_m;

  //
  // Grid function (position vectors of the elements).
  //

  if (!ap.check_array_annu("dream_arr_annu", py_G)) {
    throw std::runtime_error("Error in dream_arr_annu!");
  }

  py::buffer_info G_info = py_G->request();
  auto G_m = G_info.shape[0];
  double *G = static_cast<double*>(G_info.ptr);

  dream_idx_type num_elements = (dream_idx_type) G_m;

  //
  // Temporal and spatial sampling parameters.
  //

  double dx=0.0, dy=0.0, dt=0.0;
  dream_idx_type nt=0;
  if (!ap.parse_sampling("dream_arr_annu", py_s_par, 4, dx, dy, dt, nt)) {
    throw std::runtime_error("Error in dream_arr_annu!");
  }

  //
  // Start point of impulse response vector ([us]).
  //

  if (!ap.check_delay("dream_arr_annu", py_delay, no)) {
    throw std::runtime_error("Error in dream_arr_annu!");
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
  if (!ap.parse_material("dream_arr_annu", py_m_par, v, cp, alpha)) {
    throw std::runtime_error("Error in dream_arr_annu!");
  }

  //
  // Focusing parameters.
  //

  FocusMet foc_met=FocusMet::none;

  // Allocate memory for the user defined focusing delays
  std::unique_ptr<double[]> focal = std::make_unique<double[]>(num_elements);

  if (!ap.parse_focus_args("dream_arr_annu", foc_str, py_focal, num_elements, foc_met, focal.get())) {
    throw std::runtime_error("Error in dream_arr_annu!");
  }

  //
  // Apodization.
  //

  bool    do_apod=false;
  ApodMet apod_met=ApodMet::gauss;

  // Allocate space for the user defined apodization weights
  std::unique_ptr<double[]> apod = std::make_unique<double[]>(num_elements);

  double apod_par=0.0;
  if (!ap.parse_apod_args("dream_arr_annu", apod_str, py_apod, num_elements,
                          do_apod, apod.get(), apod_met, apod_par)) {
    throw std::runtime_error("Error in dream_arr_annu!");
  }

  //
  // Error reporting.
  //

  ErrorLevel err_level=ErrorLevel::stop;

  if (!ap.parse_error_arg("dream_arr_annu", err_level_str, err_level)) {
    throw std::runtime_error("Error in dream_arr_annu!");
  }

  // Create an output matrix for the impulse response
  auto py_h_mat = py::array_t<double,py::array::f_style>({nt,no});
  py::buffer_info h_mat_info = py_h_mat.request();
  double *h = static_cast<double*>(h_mat_info.ptr);

  SIRData hsir(h, nt, no);
  hsir.clear();

  // Register signal handler.
  std::signal(SIGINT, ArrAnnu::abort);

  //
  // Call DREAM function
  //

  ArrAnnu arr_annu;

  SIRError err = arr_annu.dream_arr_annu(alpha,
                                        ro, no,
                                        dx, dy,  dt, nt,
                                        delay_type, delay,
                                        v, cp,
                                        num_elements, G,
                                        foc_met, focal.get(),
                                        apod.get(), do_apod, apod_met, apod_par,
                                        h, err_level);

  if (!arr_annu.is_running()) {
    if (err == SIRError::out_of_bounds) {
      dream_err_msg("SIR out-of-bounds!\n"); // Bail out.
    } else {
      dream_err_msg("CTRL-C pressed!\n"); // Bail out.
    }

    throw std::runtime_error("Error in dream_arr_annu!");
  }

  if (err == SIRError::warn_out_of_bounds) {
    std::cout << "Warning: SIR out-of-bounds!" << std::endl;
  }

  return py_h_mat;
}

/***
 *
 *  Python gateway function for (parallel) dream_arr_annu.
 *
 ***/

PYBIND11_MODULE(dream_arr_annu, m)
{
  m.doc() = "H = dream_arr_annu(Ro,G,s_par,delay,m_par,foc_met,\n\
                focal,apod_met,apod,win_par,err_level)\n\
\n\
DREAM_ARR_ANNU - Computes the spatial impulse response\n\
for an annular array using parallel processing (using threads).\n\
\n\
Observation point(s) ([mm]):\n\
\n\
'Ro'\n\
   An N x 3 matrix, Ro = np.asmatrix([xo1 yo1 zo2; xo2 yo2 zo2; ... xoN yoN zoN]) where N is the number of observation points.\n\
\n\
Array grid parameter:\n\
\n\
'G'\n\
   Vector of annulus radii [mm].\n\
\n\
Sampling parameters: s_par = np.asmatrix([dx,dy,dt,nt]) \n\
a\n\
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
   Focusing method, options are: \"off\", \"on\", and \"ud\".\n\
'focal'\n\
   Focal distance [mm]. If foc_met = \"ud\" (user defined) then focal is a vector of focusing delays.\n\
\n\
Apodization parameters: apod_met, apod, and win_par. The apod_met (apodization method) options are:\n\
\n\
\"off\"\n\
   No apodization.\n\
\"ud\"\n\
   User defined apodization.\n\
 \"triangle\"\n\
   Triangle window.\n\
\"gauss\"\n\
   Gaussian (bell-shaped) window.\n\
\"raised\"\n\
   Raised cosine.\n\
\"hann\"\n\
   Hann window.\n\
\"hamming\"\n\
   Hamming window.\n\
\"simply\"\n\
   Simply supported.\n\
\"clamped\"\n\
   Clamped.\n\
\n\
and the apod and win_par parameters are:\n\
\n\
'apod'\n\
   Vector of apodiztion weights (used for the \"ud\" option).\n\
'win_par'\n\
   A scalar parameter for raised cosine and Gaussian apodization functions.\n\
\n\
Error Handling: err_level\n\
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
dream_arr_annu is a part of the DREAM Toolbox available at\n\
https://github.com/frli8848/DREAM.\n\
\n\
Copyright 2006-2023 Fredrik Lingvall.";

  m.def("dream_arr_annu", &py_dream_arr_annu,
        "H = dream_arr_annu(Ro,G,s_par,delay,m_par,foc_met,\n\
                focal,apod_met,apod,win_par,err_level)",
        py::arg("Ro"),
        py::arg("G"),
        py::arg("s_par"),
        py::arg("delay"),
        py::arg("m_par"),
        py::arg("foc_met"),
        py::arg("focal"),
        py::arg("apod_met"),
        py::arg("apod"),
        py::arg("win_par"),
        py::arg("err_level"));

  //py::arg("H"),
  //py::arg("err"),

  //m.attr("H") = world;

}
