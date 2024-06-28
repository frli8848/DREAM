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

std::unique_ptr<DAS<double>> das_d=nullptr;
//std::unique_ptr<DAS<float>> das_f=nullptr; TODO: not implemented yet

// NB. We only have implemented double support here.

namespace jl = jlcxx;

jl::ArrayRef<double, 2> jl_das(jl::ArrayRef<double, 2> jl_Y,
                               jl::ArrayRef<double,2> jl_Gt,
                               jl::ArrayRef<double,2> jl_Gr,
                               jl::ArrayRef<double,2> jl_Ro,
                               double dt,
                               jl::ArrayRef<double> jl_delay,
                               double cp,
                               std::string das_method_str,
                               std::string err_level_str,
                               std::string device_str);

jl::ArrayRef<double, 2> jl_das(jl::ArrayRef<double, 2> jl_Y,
                               jl::ArrayRef<double,2> jl_Gt,
                               jl::ArrayRef<double,2> jl_Gr,
                               jl::ArrayRef<double,2> jl_Ro,
                               double dt,
                               jl::ArrayRef<double> jl_delay,
                               double cp,
                               std::string das_method_str,
                               std::string err_level_str,
                               std::string device_str)
{
  ArgParser ap;

  //
  // Data
  //

  dream_idx_type a_scan_len = (dream_idx_type) ap.get_m(jl_Y); // A-scan length
  dream_idx_type num_a_scans =  (dream_idx_type) ap.get_m(jl_Y);

  double *Y = static_cast<double*>(ap.get_data(jl_Y));

  //
  // Transmit array
  //

  if (ap.get_n(jl_Gt) != 3) {
    throw std::runtime_error("Argument 2 must be a (num transmit elements) x 3 matrix!");
  }

  dream_idx_type num_t_elements = (dream_idx_type) ap.get_m(jl_Gt);
  const double *Gt =  static_cast<double*>(ap.get_data(jl_Gt));

  //
  // Recieve array
  //

  dream_idx_type num_r_elements = 0;
  const double *Gr = nullptr;

  if (ap.get_m(jl_Gr) != 0) {         // SAFT is using an empty Gr.

    if (ap.get_n(jl_Gr) != 3) {
      throw std::runtime_error("Argument 3 must be a (num recieve elements) x 3 matrix!");
    }

    num_r_elements =  (dream_idx_type) ap.get_m(jl_Gr);
    Gr = static_cast<double*>(ap.get_data(jl_Gr));
  }

  //
  // Observation point.
  //

  if (!ap.check_obs_points("das", jl_Ro)) {
    throw std::runtime_error("Obs point error in das!");
  }

  double *Ro = static_cast<double*>(ap.get_data(jl_Ro));
  dream_idx_type No = (dream_idx_type) ap.get_m(jl_Ro);

  //
  // Temporal sampling parameter.
  //

  // No arg parsing needed here. Or, perhaps check for reasonbale values?

  //
  // Start point of impulse response vector ([us]).
  //

  if (!ap.check_delay("das", jl_delay, No)) {
    throw std::runtime_error("Error in das!");
  }

  double *delay = static_cast<double*>(ap.get_data(jl_delay));

  DelayType delay_type = DelayType::single; // delay is a scalar.
  if (ap.get_len(jl_delay) != 1) {
    delay_type = DelayType::multiple; // delay is a vector.
  }

  //
  // Material parameter
  //

  // No arg parsing needed here. Or, perhaps check for reasonbale values?

  //
  // DAS method
  //

  DASType das_type;
  if (!ap.parse_das_arg("das", das_method_str, das_type)) {
    throw std::runtime_error("Error in das!");
  }

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

  //FIXME: We have no check for allowed devices names yet!

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

  jl_value_t *matrix_type = jl_apply_array_type((jl_value_t*) jl_float64_type,2);
  jl_array_t *jl_Im_array = jl_alloc_array_2d(matrix_type, No, 1);
  double *Im = (double *) jl_array_data(jl_Im_array);
  auto jl_Im_mat = jl::ArrayRef<double, 2>(jl_Im_array);

  //
  // Init DAS and output arg.
  //

  SIRError err = SIRError::none;

  if (das_d) { // das (double) object exist - check if we can reuse previous das init
    if (!das_d->das_setup_has_changed(das_type, a_scan_len, No, num_t_elements, num_r_elements, use_gpu)) {
      init_das = false;
    }
  }

  if (init_das) {

    try {
      das_d = std::make_unique<DAS<double>>(das_type, a_scan_len, No, num_t_elements, num_r_elements, use_gpu);
    }

    catch (std::runtime_error &err) {
      std::cout << err.what();
      throw std::runtime_error("Failed to create an DAS object in das!");
    }
  }

  das_d->set_running();

  // Register signal handler.
  std::signal(SIGABRT, DAS<double>::abort);

#ifdef USE_OPENCL

  // Check if we should use the GPU
  if (device_str == "gpu") {

    das_d->cl_das(Y, Ro, Gt, Gr, dt, delay[0], cp, Im);

  } else { // Otherwise use the cpu

#endif

    if (device_str == "gpu") {
      std::cout << "Warning: Compute device set to 'gpu' but DREAM is build without OpenCL support!" << std::endl;
      std::cout << "Using the CPU backend!" << std::endl;
    }

    err = das_d->das(Y, Ro, Gt, Gr, dt, delay_type, delay, cp, Im, err_level);
    if (!das_d->is_running()) {
      if (err != SIRError::out_of_bounds) {
        dream_err_msg("CTRL-C pressed!\n"); // Bail out.
      } else {
        dream_err_msg("SIR out-of-bounds!\n"); // Bail out.
      }

      throw std::runtime_error("CTRL-C pressed!\n"); // Bail out.
    }

#ifdef USE_OPENCL
  }
#endif

  if (err == SIRError::out_of_bounds) {
    throw std::runtime_error("Error in DAS"); // Bail out if error.
  }

  return jl_Im_mat;
}

>>>>>>> 21307a2 (Fixed processing stop on SIR out-of-bounds errors for the Julia bindings)
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
