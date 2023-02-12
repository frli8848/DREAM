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
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 *
 ***/

#include "das.h"

// Support doubles and floats
template class DAS<double>;
template class DAS<float>;

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
