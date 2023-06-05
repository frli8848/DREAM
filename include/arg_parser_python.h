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
* Free
* Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
* 02110-1301, USA.
*
***/

#pragma once

#include <sstream>
#include <cstring>
#include <cstdarg>
#include <memory>

#include "dream.h"
#include "dream_error.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

class ArgParser
{
public:

  ArgParser() = default;
  ~ArgParser() = default;

  /*
  bool check_arg_in(const char *func_name, dream_idx_type nrhs, dream_idx_type min_in, dream_idx_type max_in) {
    bool retval=true;
    std::ostringstream s;
    if ( (nrhs <  min_in) || (nrhs > max_in) ) {
      if (min_in == max_in) {
        s << func_name <<  " requires " << min_in << " input arguments!";
      } else {
        s << func_name <<  " requires " << min_in << " to " << max_in << " input arguments!";
      }
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };
  */
  /*
  bool check_arg_out(const char *func_name, dream_idx_type nrhs, dream_idx_type min_out, dream_idx_type max_out) {
    bool retval=true;
    std::ostringstream s;
    if ( (nrhs <  min_out) || (nrhs > max_out) ) {
      if (min_out == max_out) {
        s << func_name <<  " requires " << min_out << " output arguments!";
      } else {
        s << func_name <<  " requires " << min_out << " to " << max_out << " output arguments!";
      }
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };
  */

  bool check_obs_points(const char *func_name, py::array_t<double,py::array::f_style> *py_ro) {
    bool retval=true;
    std::ostringstream s;

    py::buffer_info ro_info = py_ro->request();
    auto ro_n = ro_info.shape[1];

    if (ro_n != 3) {
      s << func_name <<  " requires that observation point arg must be a (number of observation points) x 3 matrix!";
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };

  bool check_array(const char *func_name, py::array_t<double,py::array::f_style> *py_G) {
    bool retval=true;
    std::ostringstream s;

    py::buffer_info G_info = py_G->request();
    auto G_ndim = G_info.ndim;
    auto G_m = G_info.shape[0];
    auto G_n = G_info.shape[1];

    if ( G_n != 3 ) {
      s << func_name <<  " requires that array grid arg must be a (number of array elements) x 3 matrix!";
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };

  bool check_array_annu(const char *func_name, py::array_t<double,py::array::f_style> *py_G) {
    bool retval=true;
    std::ostringstream s;

    if ( (get_m(py_G) != 1) && (get_n(py_G) != 1) &&
         (get_n(py_G)*get_n(py_G) % 2 == 0) ) {
      s << func_name <<  " requires that array grid arg must be a (number of radii) vector with an odd number of elements!";
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };

  /*
  bool check_geometry(const char *func_name, args_t args, dream_idx_type arg_num, dream_idx_type num_pars) {
    bool retval=true;
    std::ostringstream s;
    if (!((get_m(args, arg_num)==num_pars && get_n(args, arg_num)==1) || (get_m(args, arg_num)==1 && get_n(args, arg_num)==num_pars))) {
      if (num_pars == 1) {
        s << func_name <<  " requires that arg " << arg_num+1 << " (geometry) must be a scalar!";
      } else {
        s << func_name <<  " requires that arg " << arg_num+1 << " (geometry) must be a " << num_pars << " element vector!";
      }
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };
  */

  bool parse_geometry(const char *func_name, py::array_t<double,py::array::f_style> *py_geom_pars, dream_idx_type num_pars,
                      double &arg1, double &arg2, double &arg3) {
    bool retval=true;
    std::ostringstream s;

    if (!((get_m(py_geom_pars) == num_pars && get_n(py_geom_pars) == 1) ||
          (get_m(py_geom_pars) == 1 && get_n(py_geom_pars) == num_pars))) {
      if (num_pars == 1) {
        s << func_name <<  " requires that the geometry arg must be a scalar!";
      } else {
        s << func_name <<  " requires that the geometry arg must be a " << num_pars << " element vector!";
      }
      dream_err_msg(s.str().c_str());
      retval=false;
    } else {

      double *pars = get_data(py_geom_pars);

      if (num_pars==1) { // Scalars needs special treatment
        arg1 = get_scalar(py_geom_pars);
      }

      if (num_pars>1) {
        arg1 = pars[0];
        arg2 = pars[1];
      }

      if (num_pars>2) {
        arg3 = pars[2];
      }

    }

    return retval;
  };

  /*
  bool check_sampling(const char *func_name, py::array_t<double,py::array::f_style> *py_s_par, dream_idx_type num_pars) {
    bool retval=true;
    std::ostringstream s;
    if (!((get_m(py_s_par)==num_pars && get_n(py_s_pa)==1) ||
          (get_m(py_s_par)==1 && get_n(py_s_par)==num_pars))) {
      s << func_name <<  " requires that sampling par arg must be a " << num_pars << " element vector!";
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };
  */

  bool parse_sampling(const char *func_name, py::array_t<double,py::array::f_style> *py_s_par, dream_idx_type num_pars,
                      double &dx, double &dy, double &dt, dream_idx_type &nt) {
    bool retval=true;
    std::ostringstream s;

    if (!((get_m(py_s_par) == num_pars && get_n(py_s_par) == 1) ||
          (get_m(py_s_par) == 1 && get_n(py_s_par) == num_pars))) {
      s << func_name <<  " requires that sampling par arg must be a " << num_pars << " element vector!";
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    if (retval) {
      dx = get_data(py_s_par)[0]; // Spatial x-direction discretization size.
      dy = get_data(py_s_par)[1]; // Spatial y-direction discretization size.
      dt = get_data(py_s_par)[2]; // Temporal discretization size (= 1/sampling freq).
      nt = (dream_idx_type) get_data(py_s_par)[3]; // Length of SIR.
    }

    return retval;
  };


  /*
  bool check_material(const char *func_name, args_t args, dream_idx_type arg_num, dream_idx_type num_pars) {
    bool retval=true;
    std::ostringstream s;
    if (!((get_m(args, arg_num)==num_pars && get_n(args, arg_num)==1) || (get_m(args, arg_num)==1 && get_n(args, arg_num)==num_pars))) {
      s << func_name <<  " requires that arg " << arg_num+1 << " (material parameters) must be a " << num_pars << " element vector!";
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };
  */

  bool parse_material(const char *func_name, py::array_t<double,py::array::f_style> *py_m_par,
                      double &v, double &cp, double &alpha) {
    bool retval=true;
    std::ostringstream s;
    if (!((get_m(py_m_par) == 3 && get_n(py_m_par) == 1) ||
          (get_m(py_m_par) == 1 && get_n(py_m_par) == 3))) {
      s << func_name <<  " requires that material parameters par arg must be a 3 element vector!";
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    if (retval) {
      v     = get_data(py_m_par)[0];
      cp    = get_data(py_m_par)[1];
      alpha = get_data(py_m_par)[2];
    }

    return retval;
  };


  // Check that the delay is a scalar or vector of length no.
  bool check_delay(const char *func_name, py::array_t<double,py::array::f_style> *py_focal, dream_idx_type No) {
    bool retval=true;
    std::ostringstream s;

    dream_idx_type delay_len = get_m(py_focal)*get_n(py_focal);
    if ( (delay_len > 1) && (delay_len != No) ) { // delay is a vector

      if (!((get_m(py_focal) == delay_len && get_n(py_focal) == 1) ||
            (get_m(py_focal) == 1 && get_n(py_focal) == delay_len))) {
        s << func_name <<  " requires that delay arg must be a " << delay_len << " element vector";
        s << " (with a length equal to the number of observation points) or a scalar!";
      }
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };

  // A string arg
  bool parse_error_arg(const char *func_name, std::string err_level_str, ErrorLevel &err_level) {
    bool retval=true;
    std::ostringstream s;

    bool is_set = false;

    if (err_level_str == "ignore") {
      err_level = ErrorLevel::ignore;
      is_set = true;
    }

    if (err_level_str == "warn") {
      err_level = ErrorLevel::warn;
      is_set = true;
    }

    if (err_level_str == "stop") {
      err_level = ErrorLevel::stop;
      is_set = true;
    }

    if (is_set == false) {
      retval=false;
      s << func_name <<  " Unknown error level!";
      dream_err_msg(s.str().c_str());
    }

    return retval;
  };

  /***
   *
   * Read two input args:
   *
   * For arrays: A string arg and a vector or scalar
   * For rect_f and circ_f : A string arg and a scalar
   *
   ***/

  bool parse_focus_args(const char *func_name, std::string foc_str, py::array_t<double,py::array::f_style> *py_focal, dream_idx_type num_elements,
                        FocusMet &foc_met, double *focal)
  {
    bool retval=true;
    std::ostringstream s;

    bool is_set = false;

    if (foc_str == "off") {
      foc_met = FocusMet::none;
      is_set = true;
    }

    if (foc_str == "x") {
      foc_met = FocusMet::x;
      is_set = true;
    }

    if (foc_str == "y") {
      foc_met = FocusMet::y;
      is_set = true;
    }

    if ( (foc_str == "xy") || (foc_str == "on") ) { // Annular array uses "on".
      foc_met = FocusMet::xy;
      is_set = true;
    }

    if (foc_str == "x+y") {
      foc_met = FocusMet::x_y;
      is_set = true;
    }

    if (foc_str == "ud") {
      foc_met = FocusMet::ud;
      is_set = true;

      // Vector of focusing delays.
      if (get_m(py_focal) * get_n(py_focal) != num_elements) {
        s << func_name <<  " the time delay vector arg for user defined ('ud') focusing\n!";
        s << "delays must have the same length as the number of array elements!";
        retval=false;
      } else { // Copy the user defined apodization weights.
        double *tmp_focal = get_data(py_focal);
        if (focal) {
          std::memcpy(focal, tmp_focal, num_elements*sizeof(double));
        } else {
          s << func_name << " focal vector not initialzed!";
          retval=false;
        }
      }
    } else { // Here focus arg must be a scalar.

      if (get_m(py_focal)*get_n(py_focal) != 1) {
        s << func_name << " arg must be a scalar for non-user defined focusing!";
        dream_err_msg(s.str().c_str());
        retval=false;
      } else {
        if (focal) {
          focal[0] = get_scalar(py_focal);
        } else {
          s << func_name << " focal vector not initialzed!";
          dream_err_msg(s.str().c_str());
          retval=false;
        }
      }
    }

    if (is_set == false) {
      retval=false;
      s << func_name <<  " Unknown focusing method!";
      dream_err_msg(s.str().c_str());
    }

    return retval;
  };

  // A string arg and a 2 element vector
  bool parse_steer_args(const char *func_name, std::string steer_str, py::array_t<double,py::array::f_style> *py_steer_par,
                        SteerMet &steer_met, double &theta, double &phi)
  {
    bool retval=true;
    std::ostringstream s;
    bool is_set = false;

    if (steer_str == "off") {
      steer_met = SteerMet::none;
      is_set = true;
    }

    if (steer_str == "x") {
      steer_met = SteerMet::x;
      is_set = true;
    }

    if (steer_str == "y") {
      steer_met = SteerMet::y;
      is_set = true;
    }

    if (steer_str == "xy") {
      steer_met = SteerMet::xy;
      is_set = true;
    }

    if (is_set == false) {
      retval=false;
      s << func_name <<  " Unknown beem steering method!";
      dream_err_msg(s.str().c_str());
    }

    // Check that steer arg is a 2 element vector.
    if (!((get_m(py_steer_par) == 2 && get_n(py_steer_par) == 1) ||
          (get_m(py_steer_par) == 1 && get_n(py_steer_par) == 2)))
      {
        retval=false;
        s << func_name <<  " steer_par must be a 2 element vector!";
        dream_err_msg(s.str().c_str());
      }

    if (is_set && retval) {
      const double *steer_pars = get_data(py_steer_par);
      theta = steer_pars[0];
      phi= steer_pars[1];
    }

    return retval;
  };

  // A string arg, a vector, and a scalar.
  bool parse_apod_args(const char *func_name,  std::string apod_str, py::array_t<double,py::array::f_style> *py_apod, dream_idx_type num_elements,
                       bool &do_apod, double *apod, ApodMet &apod_met, double &apod_par) {
    bool retval=true;
    std::ostringstream s;
    bool is_set = false;

    do_apod = false;			// default off.

    if (apod_str == "off") {
      do_apod = false;
      is_set = true;
    }

    if (apod_str == "triangle") {
      do_apod = true;
      apod_met = ApodMet::triangle;
      is_set = true;
    }

    if (apod_str == "gauss") {
      do_apod = true;
      apod_met = ApodMet::gauss;
      is_set = true;
    }

    if (apod_str == "raised") {
      do_apod = true;
      apod_met = ApodMet::raised_cosine;
      is_set = true;
    }

    if (apod_str == "hann") {
      do_apod = true;
      apod_met = ApodMet::hann;
      is_set = true;
    }

    if (apod_str == "hamming") {
      do_apod = true;
      apod_met = ApodMet::hamming;
      is_set = true;
    }

    if (apod_str == "simply") {
      do_apod = true;
      apod_met = ApodMet::simply_supported;
      is_set = true;
    }

    if (apod_str == "clamped") {
      do_apod = true;
      apod_met = ApodMet::clamped;
      is_set = true;
    }

    if (apod_str == "ud") {
      do_apod = true;
      apod_met = ApodMet::ud;
      is_set = true;

      // Vector of apodization weights.
      if (get_m(py_apod) * get_n(py_apod) != num_elements) {
        s << func_name <<  " the length of apodization vector argument must be the same as the number of array elements!";
        retval=false;
      } else { // Copy the user defined apodization weights.
        double *tmp_apod = get_data(py_apod);
        if (apod) {
          std::memcpy(apod, tmp_apod, num_elements*sizeof(double));
        } else {
          s << func_name << " apodization vector not initialzed!";
          retval=false;
        }
      }
    }

    if (is_set == false) {
      retval=false;
      s << func_name <<  " Unknown apodization method!";
      dream_err_msg(s.str().c_str());
    }

    // Check that apod par is a scalar.
    /*
      if (!(get_m(args, arg_num+2)==1 && get_n(args, arg_num+2)==1)) {
      retval=false;
      s << func_name <<  " arg " << arg_num+3 << " must be a scalar!";
      dream_err_msg(s.str().c_str());
      }
    */

    /*
    if (is_set && retval) {
      apod_par = get_scalar(args, arg_num+2);
    }
    */

    return retval;
  };

private:

  dream_idx_type get_m(const py::array_t<double,py::array::f_style> *py_arg) {
    py::buffer_info py_arg_info = py_arg->request();
    auto m = py_arg_info.shape[0];

    return (dream_idx_type) m;
  };

  dream_idx_type get_n(const py::array_t<double,py::array::f_style> *py_arg) {
    py::buffer_info py_arg_info = py_arg->request();
    auto n = py_arg_info.shape[1];

    return (dream_idx_type) n;
  };

  double* get_data(const py::array_t<double,py::array::f_style> *py_arg) {
    py::buffer_info py_arg_info = py_arg->request();

    return static_cast<double*>(py_arg_info.ptr);
  };

  double get_scalar(const py::array_t<double,py::array::f_style> *py_arg) {
    py::buffer_info py_arg_info = py_arg->request();
    double *retval = static_cast<double*>(py_arg_info.ptr);

    return retval[0];
  };

  /*
  bool is_string (args_t args, dream_idx_type arg_num) {
    return args(arg_num).is_string();
  };
  */
  /*
  std::string get_string_arg(args_t args, dream_idx_type arg_num) {
    return args(arg_num).string_value();
  };
  */

};
