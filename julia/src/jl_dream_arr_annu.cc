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

#include "dream_arr_annu.h"
#include "arg_parser_julia.h"

namespace jl = jlcxx;

jl::ArrayRef<double, 2> jl_dream_arr_annu(jl::ArrayRef<double, 2> jl_ro,
                                          jl::ArrayRef<double> jl_G,
                                          jl::ArrayRef<double> jl_s_par,
                                          jl::ArrayRef<double> jl_delay,
                                          jl::ArrayRef<double> jl_m_par,
                                          std::string foc_str, jl::ArrayRef<double> jl_focal,
                                          std::string apod_str, jl::ArrayRef<double> jl_apod, double win_par,
                                          std::string err_level_str);

jl::ArrayRef<double, 2> jl_dream_arr_annu(jl::ArrayRef<double, 2> jl_ro,
                                          jl::ArrayRef<double> jl_G,
                                          jl::ArrayRef<double> jl_s_par,
                                          jl::ArrayRef<double> jl_delay,
                                          jl::ArrayRef<double> jl_m_par,
                                          std::string foc_str, jl::ArrayRef<double> jl_focal,
                                          std::string apod_str, jl::ArrayRef<double> jl_apod, double win_par,
                                          std::string err_level_str)
{
  ArgParser<double> ap;

  //
  // Observation points.
  //

  if (!ap.check_obs_points("dream_arr_annu", jl_ro)) {
    throw std::runtime_error("Obs point error in dream_arr_annu!");
  }

  double *ro = static_cast<double*>(ap.get_data(jl_ro));
  dream_idx_type no = (dream_idx_type) ap.get_m(jl_ro);

  //
  // Grid function (position vectors of the elements).
  //

  if (!ap.check_array_annu("dream_arr_annu", jl_G)) {
    throw std::runtime_error("Error in dream_arr_annu!");
  }

  auto G_m = ap.get_len(jl_G);
  double *G = ap.get_data(jl_G);

  dream_idx_type num_elements = (dream_idx_type) G_m;

  //
  // Temporal and spatial sampling parameters.
  //

  double dx=0.0, dy=0.0, dt=0.0;
  dream_idx_type nt=0;
  if (!ap.parse_sampling("dream_arr_annu", jl_s_par, 4, dx, dy, dt, nt)) {
    throw std::runtime_error("Sampling par rrror in dream_arr_annu!");
  }

  //
  // Start point of impulse response vector ([us]).
  //

  if (!ap.check_delay("dream_arr_annu", jl_delay, no)) {
    throw std::runtime_error("Error in dream_arr_annu!");
  }

  double *delay = static_cast<double*>(ap.get_data(jl_delay));

  DelayType delay_type = DelayType::single; // delay is a scalar.
  if (ap.get_len(jl_delay) != 1) {
    delay_type = DelayType::multiple; // delay is a vector.
  }

  //
  // Material parameters
  //

  double v=1.0, cp=1000.0, alpha=0.0;
  if (!ap.parse_material("dream_arr_annu", jl_m_par, v, cp, alpha)) {
    throw std::runtime_error("Error in dream_arr_annu!");
  }

  //
  // Focusing parameters.
  //

  FocusMet foc_met=FocusMet::none;

  // Allocate memory for the user defined focusing delays
  std::unique_ptr<double[]> focal = std::make_unique<double[]>(num_elements);

  if (!ap.parse_focus_args("dream_arr_annu", foc_str, jl_focal, num_elements, foc_met, focal.get())) {
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
  if (!ap.parse_apod_args("dream_arr_annu", apod_str, jl_apod, num_elements,
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

  // Create an output matrix for the impulse responses
  jl_value_t *matrix_type = jl_apply_array_type((jl_value_t*) jl_float64_type,2);
  jl_array_t *jl_h_array = jl_alloc_array_2d(matrix_type, nt, no);
  double *h = (double *) jl_array_data(jl_h_array);
  auto jl_h_mat = jl::ArrayRef<double, 2>(jl_h_array);

  SIRData hsir(h, nt, no);
  hsir.clear();

  ArrAnnu arr_annu;

  // Register signal handler.
  std::signal(SIGINT, ArrAnnu::abort);

  //
  // Call DREAM function
  //

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

  return jl_h_mat;
}

/***
 *
 *  Julia gateway function for (parallel) dream_arr_annu.
 *
 ***/

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.method("dream_arr_annu", jl_dream_arr_annu);
}
