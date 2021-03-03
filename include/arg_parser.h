/***
*
* Copyright (C) 2021 Fredrik Lingvall
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

#include "dream.h"
#include "dream_error.h"

#ifdef DREAM_OCTAVE
#include <octave/oct.h>
typedef octave_value_list args_t;
#define GET_M(N)   args(N).matrix_value().rows()
#define GET_N(N)   args(N).matrix_value().cols()
#define IS_STRING(N) args(N).is_string()

std::string get_string_arg(args_t args, dream_idx_type arg_num) {
  return args(arg_num).string_value();
};

#endif

#ifdef DREAM_MATLAB
#include "mex.h"

typedef const mxArray** args_t;
#define GET_M(N)   mxGetM(args[N])
#define GET_N(N)   mxGetN(args[N])
#define GET_PR(N)  mxGetPr(args[N])
#define IS_STRING(N) mxIsChar(args[N])

std::string get_string_arg(args_t args, dream_idx_type arg_num) {
  std::ostringstream s;
  char tmp_str[100];
  dream_idx_type buflen = (GET_M(arg_num) * GET_N(arg_num) * sizeof(mxChar)) + 1;
  if (buflen < 100) {
    mxGetString(args[arg_num], tmp_str, buflen);
    s << tmp_str;
  } else {
    s << "";
  }
  return s.str();
}
#endif

class ArgParser
{
 public:

  ArgParser() = default;
  ~ArgParser() = default;

  bool check_arg_in(const char *func_name, dream_idx_type nrhs, dream_idx_type min_in, dream_idx_type max_in) {
    bool retval=true;
    std::ostringstream s;
    if ( (nrhs <  min_in) || (nrhs > max_in) ) {
      s << func_name <<  " requires " << min_in << " to " << max_in << " input arguments!";
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };

  bool check_arg_out(const char *func_name, dream_idx_type nrhs, dream_idx_type min_in, dream_idx_type max_in) {
    bool retval=true;
    std::ostringstream s;
    if ( (nrhs <  min_in) || (nrhs > max_in) ) {
      s << func_name <<  " requires " << min_in << " to " << max_in << " output arguments!";
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };


  bool check_obs_points(const char *func_name, args_t args, dream_idx_type arg_num) {
    bool retval=true;
    std::ostringstream s;
    if ( GET_N(arg_num) != 3 ) {
      s << func_name <<  " requires that arg " << arg_num+1 << "  must be a (number of observation points) x 3 matrix!";
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };

  bool check_geometry(const char *func_name, args_t args, dream_idx_type arg_num, dream_idx_type num_pars) {
    bool retval=true;
    std::ostringstream s;
    if (!((GET_M(arg_num)==num_pars && GET_N(arg_num)==1) || (GET_M(arg_num)==1 && GET_N(arg_num)==num_pars))) {
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

  bool check_sampling(const char *func_name, args_t args, dream_idx_type arg_num, dream_idx_type num_pars) {
    bool retval=true;
    std::ostringstream s;
    if (!((GET_M(arg_num)==num_pars && GET_N(arg_num)==1) || (GET_M(arg_num)==1 && GET_N(arg_num)==num_pars))) {
      s << func_name <<  " requires that arg " << arg_num+1 << " (sampling) must be a " << num_pars << " element vector!";
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };

  bool check_material(const char *func_name, args_t args, dream_idx_type arg_num, dream_idx_type num_pars) {
    bool retval=true;
    std::ostringstream s;
    if (!((GET_M(arg_num)==num_pars && GET_N(arg_num)==1) || (GET_M(arg_num)==1 && GET_N(arg_num)==num_pars))) {
      s << func_name <<  " requires that arg " << arg_num+1 << " (material parameters) must be a " << num_pars << " element vector!";
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };

  // Check that the delay is a scalar or vector of length no.
  bool check_delay(const char *func_name, args_t args, dream_idx_type arg_num, dream_idx_type no) {
    bool retval=true;
    std::ostringstream s;

    dream_idx_type delay_len = GET_M(arg_num)*GET_N(arg_num);
    if ( (delay_len > 1) && (delay_len != no) ) { // delay is a vector

      if (!((GET_M(arg_num)==delay_len && GET_N(arg_num)==1) || (GET_M(arg_num)==1 && GET_N(arg_num)==delay_len))) {
        s << func_name <<  " requires that arg " << arg_num+1 << " (delay) must be a " << delay_len << " element vector";
        s << " (with a length equal to the number of observation points) or a scalar!";
      }
      dream_err_msg(s.str().c_str());
      retval=false;
    }

    return retval;
  };

  bool parse_error_arg(const char *func_name, args_t args, dream_idx_type arg_num, ErrorLevel &err_level) {
    bool retval=true;
    std::ostringstream s;
    if (!IS_STRING(arg_num)) {
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

private:

};
