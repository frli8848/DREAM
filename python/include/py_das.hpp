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

#include "arg_parser_python.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

template <class T>
py::array_t<T,py::array::f_style> py_das(py::array_t<T,py::array::f_style> *py_Y,
                                         py::array_t<T,py::array::f_style> *py_Gt,
                                         py::array_t<T,py::array::f_style> *py_Gr,
                                         py::array_t<T,py::array::f_style>*py_Ro,
                                         T dt,
                                         py::array_t<T,py::array::f_style> *py_delay,
                                         T cp,
                                         std::string das_method_str,
                                         std::string err_level_str,
                                         std::string device_str,
                                         std::string verbose_str)
{
  ArgParser<T> ap;

  //
  // Data
  //

  dream_idx_type a_scan_len = (dream_idx_type) ap.get_m(py_Y); // A-scan length

  T *Y = static_cast<T*>(ap.get_data(py_Y));

  //
  // Transmit array
  //

  if (ap.get_n(py_Gt) != 3) {
    throw std::runtime_error("Argument 2 must be a (num transmit elements) x 3 matrix!");
  }

  dream_idx_type num_t_elements = (dream_idx_type) ap.get_m(py_Gt);
  const T *Gt =  static_cast<T*>(ap.get_data(py_Gt));

  //
  // DAS method
  //

  DASType das_type;
  if (!ap.parse_das_arg("das", das_method_str, das_type)) {
    throw std::runtime_error("Error in das!");
  }

  //
  // Recieve array
  //

  dream_idx_type num_r_elements = 0;
  const T *Gr = nullptr;

  if (das_type != DASType::saft) { // We do not use Gr in SAFT mode.

    if (ap.get_n(py_Gr) != 3) {
      throw std::runtime_error("Argument 3 must be a (num recieve elements) x 3 matrix!");
    }

    num_r_elements =  (dream_idx_type) ap.get_m(py_Gr);
    Gr = static_cast<T*>(ap.get_data(py_Gr));
  }

  //
  // Observation point.
  //

  if (!ap.check_obs_points("das", py_Ro)) {
    throw std::runtime_error("Obs point error in das!");
  }

  T *Ro = static_cast<T*>(ap.get_data(py_Ro));
  dream_idx_type No = (dream_idx_type) ap.get_m(py_Ro);

  //
  // Temporal sampling parameter.
  //

  // No arg parsing needed here. Or, perhaps check for reasonbale values?

  //
  // Start point of impulse response vector ([us]).
  //

  if (!ap.check_delay("das", py_delay, No)) {
    throw std::runtime_error("Error in das!");
  }

  T *delay = static_cast<T*>(ap.get_data(py_delay));

  DelayType delay_type = DelayType::single; // delay is a scalar.
  if (ap.get_len(py_delay) != 1) {
    delay_type = DelayType::multiple; // delay is a vector.
  }

  //
  // Material parameter
  //

  // No arg parsing needed here. Or, perhaps check for reasonable values?

  //
  // Error reporting.
  //

  ErrorLevel err_level=ErrorLevel::stop;

  if (!ap.parse_error_arg("das", err_level_str, err_level)) {
    throw std::runtime_error("Error in das!");
  }

  //
  // Compute device
  //

  if (device_str != "cpu" && device_str == "gpu") {
    throw std::runtime_error("Error in das - unkown compute device string!");
  }

  //
  // Verbose status printouts
  //

  bool verbose = false;

  if (verbose_str == "verbose") {
    verbose = true;
  }

  //
  // Init DAS and output arg.
  //

  bool init_das = true;

  // We use this to make the GPU init code run
  bool use_gpu = false;
#ifdef USE_OPENCL
  if (device_str == "gpu") {
    use_gpu = true;
  }
#endif

  //
  // Create an output matrix for the processed data, Im
  //

  auto py_Im_mat = py::array_t<T,py::array::f_style>(No);
  py::buffer_info Im_mat_info = py_Im_mat.request();
  T *Im = static_cast<T*>(Im_mat_info.ptr);

  SIRError err = SIRError::none;

  //
  // Init DAS and output arg.
  //

  if (das) { // das (T) object exist - check if we can reuse previous das init
    if (!das->das_setup_has_changed(das_type, a_scan_len, No, num_t_elements, num_r_elements, use_gpu)) {
      init_das = false;
    }
  }

  if (init_das) {

    try {
      das = std::make_unique<DAS<T>>(das_type, a_scan_len, No, num_t_elements, num_r_elements, use_gpu, verbose);
    }

    catch (std::runtime_error &err) {
      std::cout << err.what();
      throw std::runtime_error("Failed to create an DAS object in das!");
    }
  }

  das->set_running();

  // Register signal handler.
  std::signal(SIGABRT, DAS<T>::abort);

#ifdef USE_OPENCL

  // Check if we should use the GPU
  if (device_str == "gpu") {

    das->cl_das(Y, Ro, Gt, Gr, dt, delay[0], cp, Im);

  } else { // Otherwise use the cpu

#endif

    if (device_str == "gpu") {
      std::cout << "Warning: Compute device set to 'gpu' but DREAM is build without OpenCL support!" << std::endl;
      std::cout << "Using the CPU backend!" << std::endl;
    }

    err = das->das(Y, Ro, Gt, Gr, dt, delay_type, delay, cp, Im, err_level);
    if (!das->is_running()) {
      if (err == SIRError::out_of_bounds) {
        dream_err_msg("DAS out-of-bounds!\n"); // Bail out.
      } else {
        dream_err_msg("CTRL-C pressed!\n"); // Bail out.
      }

      throw std::runtime_error("CTRL-C pressed!\n"); // Bail out.
    }

#ifdef USE_OPENCL
  }
#endif

  if (err == SIRError::warn_out_of_bounds) {
    std::cout << "Warning: DAS out-of-bounds!" << std::endl;
  }

  return py_Im_mat;
}
