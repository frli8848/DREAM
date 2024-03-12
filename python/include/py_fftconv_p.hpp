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
#include <iostream>

#include "arg_parser_python.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "fftconv_p.h"

namespace py = pybind11;

template <class T>
py::array_t<T,py::array::f_style> py_fftconv_p(py::array_t<T,py::array::f_style> *py_A,
                                               py::array_t<T,py::array::f_style> *py_B)
{
  ArgParser<T> ap;

  // A
  dream_idx_type A_ndim = (dream_idx_type) ap.get_ndim(py_A);
  dream_idx_type A_M = (dream_idx_type) ap.get_m(py_A);
  dream_idx_type A_N = 1;
  if (A_ndim == 2) {
    A_N =  (dream_idx_type) ap.get_n(py_A);
  }
  T *A = static_cast<T*>(ap.get_data(py_A));

  // B
  dream_idx_type B_ndim = (dream_idx_type) ap.get_ndim(py_B);
  dream_idx_type B_M = (dream_idx_type) ap.get_m(py_B);
  dream_idx_type B_N = 1;
  if (B_ndim == 2) {
    B_N =  (dream_idx_type) ap.get_n(py_B);
  }
  T *B = static_cast<T*>(ap.get_data(py_B));

  // Check that arg 2 has the correct dim (matrix or vector allowed).
  if ( B_M != 1 && B_N !=1 && B_N != A_N) {
    throw std::runtime_error("Argument 2 must be a vector or a matrix with the same number of columns as arg 1!");
  }

  if (  B_M == 1 || B_N == 1 ) { // B is a vector.
    B_M = B_M*B_N;
    B_N = 1;
  }

  FFTConvP fftconvp;

  //
  // Create an output matrix for the processed data, Im
  //

  auto py_Y = py::array_t<T,py::array::f_style>({(A_M+B_M-1),A_N});
  py::buffer_info Y_mat_info = py_Y.request();

  T *Y = static_cast<T*>(Y_mat_info.ptr);

  fftconvp.run(Y,
               A, A_M, A_N,
               B, B_M, B_N);

  return py_Y;
}
