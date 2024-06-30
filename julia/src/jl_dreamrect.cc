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

#include "dreamrect.h"
#include "arg_parser_julia.h"

namespace jl = jlcxx;

jl::ArrayRef<double, 2> jl_dreamrect(jl::ArrayRef<double, 2> jl_ro,
                                     jl::ArrayRef<double> jl_geom_par,
                                     jl::ArrayRef<double> jl_s_par,
                                     jl::ArrayRef<double> jl_delay,
                                     jl::ArrayRef<double> jl_m_par,
                                     std::string err_level_str);

jl::ArrayRef<double, 2> jl_dreamrect(jl::ArrayRef<double, 2> jl_ro,
                                     jl::ArrayRef<double> jl_geom_par,
                                     jl::ArrayRef<double> jl_s_par,
                                     jl::ArrayRef<double> jl_delay,
                                     jl::ArrayRef<double> jl_m_par,
                                     std::string err_level_str)
{
  ArgParser<double> ap;

  //
  // Observation points.
  //

  if (!ap.check_obs_points("dreamrect", jl_ro)) {
    throw std::runtime_error("Obs point error in dreamrect!");
  }

  double *ro = static_cast<double*>(ap.get_data(jl_ro));
  dream_idx_type no = (dream_idx_type) ap.get_m(jl_ro);

  //
  // Transducer geometry
  //

  double a=0.0, b=0.0, dummy=0.0;
  if (!ap.parse_geometry("dreamrect", jl_geom_par, 2, a, b, dummy)) {
    throw std::runtime_error("Geometry par error in dreamrect!");
  }

  //
  // Temporal and spatial sampling parameters.
  //

  double dx=0.0, dy=0.0, dt=0.0;
  dream_idx_type nt=0;
  if (!ap.parse_sampling("dreamrect", jl_s_par, 4, dx, dy, dt, nt)) {
    throw std::runtime_error("Sampling par rrror in dreamrect!");
  }

  //
  // Start point of impulse response vector ([us]).
  //

  if (!ap.check_delay("dreamrect", jl_delay, no)) {
    throw std::runtime_error("Error in dreamrect!");
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
  if (!ap.parse_material("dreamrect", jl_m_par, v, cp, alpha)) {
    throw std::runtime_error("Error in dreamrect!");
  }

  //
  // Error reporting.
  //

  ErrorLevel err_level=ErrorLevel::stop;

  if (!ap.parse_error_arg("dreamrect", err_level_str, err_level)) {
    throw std::runtime_error("Error in dreamrect!");
  }

  // Create an output matrix for the impulse responses
  jl_value_t *matrix_type = jl_apply_array_type((jl_value_t*) jl_float64_type,2);
  jl_array_t *jl_h_array = jl_alloc_array_2d(matrix_type, nt, no);
  double *h = (double *) jl_array_data(jl_h_array);
  auto jl_h_mat = jl::ArrayRef<double, 2>(jl_h_array);

  SIRData hsir(h, nt, no);
  hsir.clear();

  Rect rect;

  // Register signal handler.
  std::signal(SIGINT, Rect::abort);

  //
  // Call DREAM function
  //

  SIRError err = rect.dreamrect(alpha,
                                ro, no,
                                a, b,
                                dx, dy,  dt, nt,
                                delay_type, delay,
                                v, cp,
                                h, err_level);

  if (!rect.is_running()) {
    if (err == SIRError::out_of_bounds) {
      dream_err_msg("SIR out-of-bounds!\n"); // Bail out.
    } else {
      dream_err_msg("CTRL-C pressed!\n"); // Bail out.
    }

    throw std::runtime_error("Error in dreamrect!");
  }

  if (err == SIRError::warn_out_of_bounds) {
    std::cout << "Warning: SIR out-of-bounds!" << std::endl;
  }

  return jl_h_mat;
}

/***
 *
 *  Julia gateway function for (parallel) dreamrect.
 *
 ***/

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.method("dreamrect", jl_dreamrect);
}
