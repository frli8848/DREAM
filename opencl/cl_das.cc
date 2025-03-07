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

#include "das.h"

// xxd:ed kernel sources
#include "das_float.cl.h"
#include "das_double.cl.h"

template <class T>
DAS<T>::DAS(DASType das_type,
            dream_idx_type a_scan_len,
            dream_idx_type No,
            dream_idx_type num_t_elements,
            dream_idx_type num_r_elements, // SAFT if num_r_elements = 0;
            bool use_gpu,
            bool verbose)
  : m_das_type(das_type)
  , m_a_scan_len(a_scan_len)
  , m_No(No)
  , m_num_t_elements(num_t_elements)
  , m_num_r_elements(num_r_elements)
  , m_use_gpu(use_gpu)
{
  if (m_use_gpu) { // We do not want to run the OpenCL init func when we use the CPU.
    if (sizeof(T) == sizeof(float) ) {
      init_opencl(__das_float_cl, __das_float_cl_len);
    } else {
      init_opencl(__das_double_cl, __das_double_cl_len);
    }
  }
}

/***
 *
 * OpenCL version of das
 *
 ***/

// SAFT (gr not used), TFM and RCA
template <class T>
int DAS<T>::cl_das(const T *Y, // Data
                   const T *Ro,
                   const T *gt,
                   const T *gr,
                   T dt,
                   T delay,
                   T cp,
                   T *Im)
{

  //
  // Transfer input vector/matrix args to their respective memory buffers
  //

  m_queue.enqueueWriteBuffer(m_cl_buf_Y, CL_TRUE, 0, m_buflen_Y, Y);
  m_queue.enqueueWriteBuffer(m_cl_buf_gt, CL_TRUE, 0, m_buflen_gt, gt);
  if (m_num_r_elements > 0) { // SAFT do not use this one.
    m_queue.enqueueWriteBuffer(m_cl_buf_gr, CL_TRUE, 0, m_buflen_gr, gr);
  }
  m_queue.enqueueWriteBuffer(m_cl_buf_Ro, CL_TRUE, 0, m_buflen_Ro, Ro);

  //
  // Set the arguments of the kernel
  //

  size_t arg_idx = 0;

  // Y
  m_kernel.setArg(arg_idx, m_cl_buf_Y);
  arg_idx++;
  m_kernel.setArg(arg_idx, (int) m_a_scan_len);

  // gt
  arg_idx++;
  m_kernel.setArg(arg_idx, m_cl_buf_gt);
  arg_idx++;
  m_kernel.setArg(arg_idx, (int) m_num_t_elements);

  // gr
  if (m_num_r_elements > 0) { // SAFT do not use this one.
    arg_idx++;
    m_kernel.setArg(arg_idx, m_cl_buf_gr);
    arg_idx++;
    m_kernel.setArg(arg_idx, (int) m_num_r_elements);
  }

  // Ro
  arg_idx++;
  m_kernel.setArg(arg_idx, m_cl_buf_Ro);
  arg_idx++;
  m_kernel.setArg(arg_idx, (int) m_No);

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

// https://bytefreaks.net/programming-2/c/c-undefined-reference-to-templated-class-function

// Support doubles and floats
template class DAS<double>;
template class DAS<float>;
