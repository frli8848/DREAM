/***
*
* Copyright (C) 2021,2022,2023 Fredrik Lingvall
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

#ifdef DREAM_OCTAVE
#include <octave/oct.h>
typedef octave_value_list args_t;
#endif

#ifdef DREAM_MATLAB
#include "mex.h"
typedef const mxArray** args_t;
#endif

//! Argument parser class for MATLAB and Octave
class ArgParser
{
public:

  ArgParser() = default;
  ~ArgParser() = default;

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

  bool check_obs_points(const char *func_name, args_t args, dream_idx_type arg_num) {
    bool retval=true;
    std::ostringstream s;
    if ( get_n(args, arg_num) != 3 ) {
      s << func_name <<  " requires that arg " << arg_num+1 << "  must be a (number of observation points) x 3 matrix!";
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };

  bool check_array(const char *func_name, args_t args, dream_idx_type arg_num) {
    bool retval=true;
    std::ostringstream s;
    if ( get_n(args, arg_num) != 3 ) {
      s << func_name <<  " requires that arg " << arg_num+1 << "  must be a (number of array elements) x 3 matrix!";
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };

  bool check_array_annu(const char *func_name, args_t args, dream_idx_type arg_num) {
    bool retval=true;
    std::ostringstream s;
    if ( (get_m(args, arg_num) != 1) && (get_n(args, arg_num) != 1) && (get_n(args, arg_num)*get_n(args, arg_num) % 2 == 0) ) {
      s << func_name <<  " requires that arg " << arg_num+1 << "  must be a (number of radii) vector with an odd number of elements!";
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };

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

  bool parse_geometry(const char *func_name, args_t args, dream_idx_type arg_num, dream_idx_type num_pars,
                      double &arg1, double &arg2, double &arg3) {
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
    } else {

      double *pars = get_data(args, arg_num);

      if (num_pars==1) { // Scalars needs special treatment
        arg1 = get_scalar(args, arg_num);
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

  bool check_sampling(const char *func_name, args_t args, dream_idx_type arg_num, dream_idx_type num_pars) {
    bool retval=true;
    std::ostringstream s;
    if (!((get_m(args, arg_num)==num_pars && get_n(args, arg_num)==1) || (get_m(args, arg_num)==1 && get_n(args, arg_num)==num_pars))) {
      s << func_name <<  " requires that arg " << arg_num+1 << " (sampling) must be a " << num_pars << " element vector!";
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };

  bool parse_sampling(const char *func_name, args_t args, dream_idx_type arg_num, dream_idx_type num_pars,
                      double &dx, double &dy, double &dt, dream_idx_type &nt) {
    bool retval=true;
    std::ostringstream s;
    if (!((get_m(args, arg_num)==num_pars && get_n(args, arg_num)==1) || (get_m(args, arg_num)==1 && get_n(args, arg_num)==num_pars))) {
      s << func_name <<  " requires that arg " << arg_num+1 << " (sampling) must be a " << num_pars << " element vector!";
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    if (retval) {
      dx = get_data(args, arg_num)[0]; // Spatial x-direction discretization size.
      dy = get_data(args, arg_num)[1]; // Spatial y-direction discretization size.
      dt = get_data(args, arg_num)[2]; // Temporal discretization size (= 1/sampling freq).
      nt = (dream_idx_type) get_data(args, arg_num)[3]; // Length of SIR.
    }

    return retval;
  };


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

  bool parse_material(const char *func_name, args_t args, dream_idx_type arg_num,
                      double &v, double &cp, double &alpha) {
    bool retval=true;
    std::ostringstream s;
    if (!((get_m(args, arg_num)==3 && get_n(args, arg_num)==1) || (get_m(args, arg_num)==1 && get_n(args, arg_num)==3))) {
      s << func_name <<  " requires that arg " << arg_num+1 << " (material parameters) must be a 3 element vector!";
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    if (retval) {
      v     = get_data(args, arg_num)[0];
      cp    = get_data(args, arg_num)[1];
      alpha = get_data(args, arg_num)[2];
    }

    return retval;
  };

  // Check that the delay is a scalar or vector of length No.
  bool check_delay(const char *func_name, args_t args, dream_idx_type arg_num, dream_idx_type No) {
    bool retval=true;
    std::ostringstream s;

    dream_idx_type delay_len = get_m(args, arg_num)*get_n(args, arg_num);
    if ( (delay_len > 1) && (delay_len != No) ) { // delay is a vector

      if (!((get_m(args, arg_num)==delay_len && get_n(args, arg_num)==1) || (get_m(args, arg_num)==1 && get_n(args, arg_num)==delay_len))) {
        s << func_name <<  " requires that arg " << arg_num+1 << " (delay) must be a " << delay_len << " element vector";
        s << " (with a length equal to the number of observation points) or a scalar!";
      }
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };

  // A string arg
  bool parse_error_arg(const char *func_name, args_t args, dream_idx_type arg_num, ErrorLevel &err_level) {
    bool retval=true;
    std::ostringstream s;
    if (!is_string(args, arg_num)) {
      s << func_name <<  " requires that arg " << arg_num+1 << " (error type) must be a string!";
      dream_err_msg(s.str().c_str());
      err_level = ErrorLevel::stop; // Default
      retval=false;
    } else {
      std::string err_str = get_string_arg(args, arg_num);
      bool is_set = false;

      if (err_str == "ignore") {
        err_level = ErrorLevel::ignore;
        is_set = true;
      }

      if (err_str == "warn") {
        err_level = ErrorLevel::warn;
        is_set = true;
      }

      if (err_str == "stop") {
        err_level = ErrorLevel::stop;
        is_set = true;
      }

      if (is_set == false) {
        retval=false;
        s << func_name <<  " Unknown error level in arg " << arg_num+1 << "!";
        dream_err_msg(s.str().c_str());
      }
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

  bool parse_focus_args(const char *func_name, args_t args, dream_idx_type arg_num,
                       FocusMet &foc_met, double *focal, dream_idx_type num_elements=1) {
    bool retval=true;
    std::ostringstream s;

    if (!is_string(args, arg_num)) {
      s << func_name <<  " requires that arg " << arg_num+1 << " (focus method) must be a string!";
      dream_err_msg(s.str().c_str());
      retval=false;
    } else {

      std::string foc_str = get_string_arg(args, arg_num);
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
        if (get_m(args, arg_num+1)*get_n(args, arg_num+1) != num_elements) {
          s << func_name <<  " the time delay vector arg " << arg_num+1 << " for user defined ('ud') focusing\n";
          s << "delays must have the same length as the number of array elements!";
          dream_err_msg(s.str().c_str());
          retval=false;
        } else { // Copy the user defined apodization weights.
          double *tmp_focal = get_data(args, arg_num+1);
          if (focal) {
            std::memcpy(focal, tmp_focal, num_elements*sizeof(double));
          } else {
            s << func_name << " focal vector not initialzed!";
            dream_err_msg(s.str().c_str());
            retval=false;
          }
        }
      } else { // Here focus arg must be a scalar.

        if (get_m(args, arg_num+1)*get_n(args, arg_num+1) != 1) {
          s << func_name << " arg " << arg_num+2 << " must be a scalar for non-user defined focusing!";
          dream_err_msg(s.str().c_str());
          retval=false;
        } else {
          if (focal) {
            focal[0] = get_scalar(args, arg_num+1);
          } else {
            s << func_name << " focal vector not initialzed!";
            dream_err_msg(s.str().c_str());
            retval=false;
          }
        }
      }

      if (is_set == false) {
        retval=false;
        s << func_name <<  " Unknown focusing method in arg " << arg_num+1 << "!";
        dream_err_msg(s.str().c_str());
      }

    }

    return retval;
  };

  // A string arg and a 2 element vector
  bool parse_steer_args(const char *func_name, args_t args, dream_idx_type arg_num,
                        SteerMet &steer_met, double &theta, double &phi) {
    bool retval=true;
    std::ostringstream s;
    if (!is_string(args, arg_num)) {
      s << func_name <<  " requires that arg " << arg_num+1 << " (focus method) must be a string!";
      dream_err_msg(s.str().c_str());
      retval=false;
    } else {
      std::string steer_str = get_string_arg(args, arg_num);
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
        s << func_name <<  " Unknown beem steering method in arg " << arg_num+1 << "!";
        dream_err_msg(s.str().c_str());
      }

      // Check that steer arg is a 2 element vector.
      if (!((get_m(args, arg_num+1)==2 && get_n(args, arg_num+1)==1) || (get_m(args, arg_num+1)==1 && get_n(args, arg_num+1)==2))) {
        retval=false;
        s << func_name <<  " arg " << arg_num+2 << " must be a 2 element vector!";
        dream_err_msg(s.str().c_str());
      }

      if (is_set && retval) {
        const double *steer_pars = get_data(args, arg_num+1);
        theta = steer_pars[0];
        phi= steer_pars[1];
      }

    }

    return retval;
  };

  // A string arg, a vector, and a scalar.
  bool parse_apod_args(const char *func_name, args_t args, dream_idx_type arg_num, dream_idx_type num_elements,
                       bool &do_apod, double *apod, ApodMet &apod_met, double &apod_par) {
    bool retval=true;
    std::ostringstream s;
    if (!is_string(args, arg_num)) {
      s << func_name <<  " requires that arg " << arg_num+1 << " (focus method) must be a string!";
      dream_err_msg(s.str().c_str());
      retval=false;
    } else {
      std::string apod_str = get_string_arg(args, arg_num);
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
        if (get_m(args, arg_num+1) * get_n(args, arg_num+1) != num_elements) {
          s << func_name <<  " the length of argument arg " << arg_num+1 << " (apodization vector) must be the same as the number of array elements!";
          retval=false;
        } else { // Copy the user defined apodization weights.
          double *tmp_apod = get_data(args, arg_num+1);
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
        s << func_name <<  " Unknown apodization method in arg " << arg_num+1 << "!";
        dream_err_msg(s.str().c_str());
      }

      // Check that apod par is a scalar.
      if (!(get_m(args, arg_num+2)==1 && get_n(args, arg_num+2)==1)) {
        retval=false;
        s << func_name <<  " arg " << arg_num+3 << " must be a scalar!";
        dream_err_msg(s.str().c_str());
      }

      if (is_set && retval) {
        apod_par = get_scalar(args, arg_num+2);
      }

    }

    return retval;
  };

# ifdef DELAY_AND_SUM

  // A string arg
  bool parse_das_arg(const char *func_name, args_t args, dream_idx_type arg_num, DASType &das_type) {
    bool retval=true;
    std::ostringstream s;
    if (!is_string(args, arg_num)) {
      s << func_name <<  " requires that arg " << arg_num+1 << " (DAS type) must be a string!";
      dream_err_msg(s.str().c_str());
      das_type = DASType::saft; // Default
      retval=false;
    } else {
      std::string das_str = get_string_arg(args, arg_num);
      bool is_set = false;

      if (das_str == "saft") {
        das_type = DASType::saft;
        is_set = true;
      }

      if (das_str == "tfm") {
        das_type = DASType::tfm;
        is_set = true;
      }

      if ( (das_str == "rca") || (das_str == "rca_coltx") ) { // Default to column transmit.
        das_type = DASType::rca_coltx;
        is_set = true;
      }

      if (das_str == "rca_rowtx") {
        das_type = DASType::rca_rowtx;
        is_set = true;
      }

      if (is_set == false) {
        retval=false;
        s << func_name <<  " Unknown error level in arg " << arg_num+1 << "!";
        dream_err_msg(s.str().c_str());
      }
    }

    return retval;
  };

#endif

#ifdef DREAM_OCTAVE

  dream_idx_type get_m(args_t args, dream_idx_type arg_num) {
    return args(arg_num).matrix_value().rows();
  };

  dream_idx_type get_n(args_t args, dream_idx_type arg_num) {
    return args(arg_num).matrix_value().cols();
  };

  double* get_data(const args_t args, dream_idx_type arg_num) {
    const Matrix tmp_mat = args(arg_num).matrix_value();
    return (double*) tmp_mat.data();
  };

  double get_scalar(const args_t args, dream_idx_type arg_num) {
    return args(arg_num).double_value();
  };

  bool is_string (args_t args, dream_idx_type arg_num) {
    return args(arg_num).is_string();
  };

  std::string get_string_arg(args_t args, dream_idx_type arg_num) {
    return args(arg_num).string_value();
  };

#endif

#ifdef DREAM_MATLAB

  dream_idx_type get_m(args_t args, dream_idx_type arg_num) {
    return mxGetM(args[arg_num]);
  };

  dream_idx_type get_n(args_t args, dream_idx_type arg_num) {
    return mxGetN(args[arg_num]);
  };

  double* get_data(args_t args, dream_idx_type arg_num) {
    return (double*) mxGetPr(args[arg_num]);
  };

  double get_scalar(const args_t args, dream_idx_type arg_num) {
    return mxGetScalar(args[arg_num]);
  };

  bool is_string (args_t args, dream_idx_type arg_num) {
    return mxIsChar(args[arg_num]);
  };

  std::string get_string_arg(args_t args, dream_idx_type arg_num) {
    std::ostringstream s;
    char tmp_str[100];
    dream_idx_type buflen = (get_m(args, arg_num) * get_n(args, arg_num) * sizeof(mxChar)) + 1;
    if (buflen < 100) {
      mxGetString(args[arg_num], tmp_str, buflen);
      s << tmp_str;
    } else {
      s << "";
    }
    return s.str();
  };
#endif
};
