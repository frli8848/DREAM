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

#include "das_uni.h"

// xxd:ed kernel sources
#include "das_uni_float.cl.h"
#include "das_uni_double.cl.h"

template <class T>
DAS_UNI<T>::DAS_UNI(DASType das_type,
                    dream_idx_type a_scan_len,   // Data size
                    T min_t, T pitch_t, T max_t, // Transmit
                    T min_r, T pitch_r, T max_r, // Receive
                    T min_Rx, T dx, T max_Rx,    // Observation points
                    T min_Ry, T dy, T max_Ry,
                    T min_Rz, T dz, T max_Rz,
                    unsigned char *cl_kernel_str,
                    unsigned int cl_kernel_str_len)
  : m_out_err(ErrorLevel::none)
  , m_das_type(das_type)
  , m_a_scan_len(a_scan_len)
{
  if (sizeof(T) == sizeof(float) ) {
    init_opencl(min_t, pitch_t, max_t,
                min_r, pitch_r, max_r,
                min_Rx, dx, max_Rx,
                min_Ry, dy, max_Ry,
                min_Rz, dz, max_Rz,
                 __das_uni_float_cl, __das_uni_float_cl_len);
  } else {
    init_opencl(min_t, pitch_t, max_t,
                min_r, pitch_r, max_r,
                min_Rx, dx, max_Rx,
                min_Ry, dy, max_Ry,
                min_Rz, dz, max_Rz,
                __das_uni_double_cl, __das_uni_double_cl_len);
  }
}

/***
 *
 * OpenCL das_uni
 *
 ***/

// SAFT (recieve not used), TFM and RCA
template <class T>
int DAS_UNI<T>::cl_das_uni(const T *Y, // Data
                           T min_t, T pitch_t, T max_t, // Transmit
                           T min_r, T pitch_r, T max_r, // Receive
                           T min_Rx, T dx, T max_Rx,    // Observation points
                           T min_Ry, T dy, T max_Ry,
                           T min_Rz, T dz, T max_Rz,
                           T dt,
                           T delay,
                           T cp,
                           T *Im)
{

  //
  // Transfer input data to the memory buffer
  //

  m_queue.enqueueWriteBuffer(m_cl_buf_Y, CL_TRUE, 0, m_buflen_Y, Y);

  //
  // Set the arguments of the kernel
  //

  size_t arg_idx = 0;

  // Y, a_scan_len
  m_kernel.setArg(arg_idx, m_cl_buf_Y);
  arg_idx++;
  m_kernel.setArg(arg_idx, (int) m_a_scan_len);

  // min_t, pitch_t, max_t
  arg_idx++;
  m_kernel.setArg(arg_idx, min_t);
  arg_idx++;
  m_kernel.setArg(arg_idx, pitch_t);
  arg_idx++;
  m_kernel.setArg(arg_idx, max_t);

  // min_r, pitch_r, max_r
  if (m_num_r_elements > 0) { // SAFT do not use this one.
    arg_idx++;
    m_kernel.setArg(arg_idx, min_r);
    arg_idx++;
    m_kernel.setArg(arg_idx, pitch_r);
    arg_idx++;
    m_kernel.setArg(arg_idx, max_r);
  }

  // min_Rx, dx, max_Rx
  arg_idx++;
  m_kernel.setArg(arg_idx, min_Rx);
  arg_idx++;
  m_kernel.setArg(arg_idx, dx);
  arg_idx++;
  m_kernel.setArg(arg_idx, max_Rx);

  // min_Ry, dy, max_Ry
  arg_idx++;
  m_kernel.setArg(arg_idx, min_Ry);
  arg_idx++;
  m_kernel.setArg(arg_idx, dy);
  arg_idx++;
  m_kernel.setArg(arg_idx, max_Ry);

  // min_Rz, dz, max_Rz
  arg_idx++;
  m_kernel.setArg(arg_idx, min_Rz);
  arg_idx++;
  m_kernel.setArg(arg_idx, dz);
  arg_idx++;
  m_kernel.setArg(arg_idx, max_Rz);

  // dt
  arg_idx++;
  m_kernel.setArg(arg_idx, dt);

  // delay
  arg_idx++;
  m_kernel.setArg(arg_idx, delay);

  // cp
  arg_idx++;
  m_kernel.setArg(arg_idx, cp);

  // Im (output arg)
  arg_idx++;
  m_kernel.setArg(arg_idx, m_cl_buf_Im);

  //
  // Execute the OpenCL kernel
  //

  cl::NDRange global(m_No);
  cl::NDRange local(64); // Works on both Nvidia and AMD GPUs.
  try {
    m_queue.enqueueNDRangeKernel(m_kernel, cl::NullRange, global, local);
  }

  catch (const cl::Error &err) {
    std::cerr << "OpenCL kernel run failed: " << err.what() << std::endl;
  }

  //
  // Transfer the memory buffer on the device to the host buffer
  //

  try {
    m_queue.enqueueReadBuffer(m_cl_buf_Im, CL_TRUE, 0,  m_buflen_Im, Im);
  }

  catch (const cl::Error &err) {
    std::cerr << "Failed to get results buffer from GPU: " << err.what() << std::endl;
  }

  int retval = 0;

  return retval;
}

// Support doubles and floats
template class DAS_UNI<double>;
template class DAS_UNI<float>;
