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

#include <thread>

#include "dream.h"
#include "attenuation.h"
#include "dream_error.h"

#include <iostream>
#include <fstream>
#include <streambuf>
#include <list>
#include <vector>

#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_CL_1_2_DEFAULT_BUILD
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY
#include <CL/opencl.hpp>

// This switches on das arg parsing
// in the arg_parser_xx:ers
#define DELAY_AND_SUM

template <class T>
class DAS_UNI
{
 public:

 DAS_UNI(DASType das_type,
         dream_idx_type a_scan_len,   // Data size
         T min_t, T pitch_t, T max_t, // Transmit
         T min_r, T pitch_r, T max_r, // Receive
         T min_Rx, T dx, T max_Rx,    // Observation points
         T min_Ry, T dy, T max_Ry,
         T min_Rz, T dz, T max_Rz,
         unsigned char *cl_kernel_str=nullptr,
         unsigned int cl_kernel_str_len=0);

  void init_opencl(T min_t, T pitch_t, T max_t, // Transmit
                   T min_r, T pitch_r, T max_r, // Receive
                   T min_Rx, T dx, T max_Rx,    // Observation points
                   T min_Ry, T dy, T max_Ry,
                   T min_Rz, T dz, T max_Rz,
                   unsigned char *cl_kernel_str, unsigned int cl_kernel_str_len)
  {
    // Array dims
    m_num_t_elements = (dream_idx_type) ((max_t - min_t)/pitch_t+1.0);
    m_num_r_elements = (dream_idx_type) ((max_r - min_r)/pitch_r+1.0);

    // Image dims
    dream_idx_type Nx = (dream_idx_type) ((max_Rx - min_Rx)/dx+1.0);
    dream_idx_type Ny = (dream_idx_type) ((max_Ry - min_Ry)/dy+1.0);
    dream_idx_type Nz = (dream_idx_type) ((max_Rz - min_Rz)/dz+1.0);
    m_No = Nx*Ny*Nz;

    std::string kernel_str;
    std::string kernel_name="das_uni_tfm"; // Default function name of the OpenCL kernel

    if (m_das_type == DASType::saft) {
      kernel_name="das_uni_saft";
      m_num_r_elements = 0;
    }

    if (m_das_type == DASType::tfm) {
      kernel_name="das_uni_tfm";
    }

    if (m_das_type == DASType::rca_coltx) {
      kernel_name="das_uni_rca_coltx";
    }

    if (m_das_type == DASType::rca_rowtx) {
      kernel_name="das_uni_rca_rowtx";
    }

    //
    // Read OpenCL kernels
    //

    if (cl_kernel_str_len == 0)  { // Fallback to read kernels from file when no xxd:ed kernel header files exist.
      std::string dream_cl_kernel = "";
      if(const char* env_p = std::getenv("DREAM_CL_KERNELS")) {
        dream_cl_kernel += env_p;
      } else {
        throw std::runtime_error("Error in das_uni - DREAM_CL_KERNELS env variable not set!");
      }

      std::string das_kernel = dream_cl_kernel;
      if (sizeof(T) == sizeof(float) ) {
        das_kernel += "/das_uni_float.cl";
      } else {
        das_kernel += "/das_uni_double.cl";
      }
      std::ifstream f_kernel(das_kernel);
      f_kernel.seekg(0, std::ios::end);
      kernel_str.reserve(f_kernel.tellg());
      f_kernel.seekg(0, std::ios::beg);
      kernel_str.assign((std::istreambuf_iterator<char>(f_kernel)),
                        std::istreambuf_iterator<char>());
    } else {
      kernel_str = std::string((char *) cl_kernel_str, cl_kernel_str_len);
    }


    //
    // Get platform
    //

    std::vector<cl::Platform> platforms;
    try {
      cl::Platform::get(&platforms);
    }

    catch (const cl::Error &err) {
      std::string the_err = "Error in das_uni - OpenCL error: " ;
      the_err += err.what();
      the_err += "\n";
      throw std::runtime_error(the_err);
    }

    if (platforms.empty()) {
      throw std::runtime_error("Error in das_uni - OpenCL error: No platforms found!");
    } else {
      std::cout << "Found " << platforms.size() << " platform(s).\n" << std::endl;
      size_t platformIndex = 0;
      for (auto &platform : platforms) {
        std::cout << "Platform[" << platformIndex++ << "]:" << std::endl;
        std::cout << "  Name: " << platform.getInfo<CL_PLATFORM_NAME>() << std::endl;
        std::cout << "  Vendor: " << platform.getInfo<CL_PLATFORM_VENDOR>() << std::endl;
      }
    }

    //
    // Select the default platform and create a context using this platform and the GPU.
    //

    cl::Context context;

    cl_context_properties properties[] = {
      CL_CONTEXT_PLATFORM,
      (cl_context_properties)(platforms[0])(),
      0};

    try {
      context = cl::Context(CL_DEVICE_TYPE_GPU, properties); // Use the GPU.
    }

    catch (const cl::Error &err) {
      std::string the_err = "Error in das_uni - OpenCL error: ";
      the_err += err.what();
      the_err += "\n";
      throw std::runtime_error(the_err);
    }

    //
    // Get a list of devices on this platform
    //

    std::vector<cl::Device> devices;

    try {
      devices = context.getInfo<CL_CONTEXT_DEVICES>();
    }

  catch (const cl::Error &err) {
    std::string the_err = "Error in das_uni - Failed to get OpenCL context info: ";
    the_err += err.what();
    the_err += "\n";
    throw std::runtime_error(the_err);
  }

  size_t dev_idx = 0;
  /*
  if (dev_idx >= devices.size()) {
   std::cerr << "Selected OpenCL device not found! Using the first device!";
    dev_idx = 0;
  }
  */

  if (devices.empty()) {
    throw std::runtime_error("Error in das_uni - OpenCL error: Found no devices(s)!");
  } else {
    std::cout << "Found " <<  devices.size() << " devices(s)." << std::endl;
    size_t deviceIndex = 0;
    for(auto &device : devices) {
      std::cout << "Device[" << deviceIndex++ << "]:" << std::endl;
      std::cout << "  Name: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
      if (device.getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_GPU) {
        std::cout << "  Type: GPU" << std::endl;
      }
      if (device.getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU) {
        std::cout << "  Type: CPU" << std::endl;
      }
      std::cout << "  Vendor: " << device.getInfo<CL_DEVICE_VENDOR>() << std::endl;
      std::cout << "  Max Compute Units: " << device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << std::endl;
      std::cout << "  Global Memory: " << device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>()<< std::endl;
      std::cout << "  Max Clock Frequency: " << device.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << std::endl;
      std::cout << "  Max Allocatable Memory: " << device.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>() << std::endl;
      std::cout << "  Local Memory: " << device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() << std::endl;
      std::cout << "  Available: " << device.getInfo<CL_DEVICE_AVAILABLE>() << std::endl;
    }
  }

  //
  // Create a command queue and use the first device.
  //

  try {
    m_queue = cl::CommandQueue(context, devices[dev_idx]);
  }

  catch (const cl::Error &err) {
    std::string the_err = "Error in das_uni - OpenCL error: ";
    the_err += err.what();
    the_err += "\n";
    throw std::runtime_error(the_err);
  }

  // Add the sources
  cl::Program::Sources source;
  source = cl::Program::Sources(1, std::make_pair(kernel_str.c_str(), kernel_str.length()+1));

  // Make program of the source code in the context
  cl::Program program;
  program = cl::Program(context, source);

  // Build program for these specific devices
  try {
    program.build(devices);
  }

  catch (const cl::Error &err) {
    std::string the_err = "Error in das_uni - OpenCL error - kernel build failed: ";
    the_err += err.what();
    the_err += "\n";
    std::cout << "Build log:" << std::endl << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[dev_idx]) << std::endl;
    throw std::runtime_error(the_err);
  }

  // Make kernel

  try {
    m_kernel = cl::Kernel(program, kernel_name.c_str());
  }

  catch (const cl::Error &err) {
    std::string the_err = "";
    std::cerr << "Error in das_uni - OpenCL error - Kernel: " << kernel_name << " creation failed: " << err.what();
    the_err += err.what();
    the_err += "\n";
    throw std::runtime_error(the_err);
  }

  //
  // Create memory buffers on the device for vector/matrix arguments to the kernel
  //

  // Read buffers

  // Data

  dream_idx_type n_r_el = m_num_r_elements;
  if (m_das_type == DASType::saft) {
    n_r_el = 1; // So we do not multiply by 0.
  }

  m_buflen_Y = m_a_scan_len*m_num_t_elements*n_r_el*sizeof(T);
  try {
    m_cl_buf_Y = cl::Buffer(context, CL_MEM_READ_ONLY, m_buflen_Y);
  }

  catch (const cl::Error &err) {
    std::string the_err = "Error in das_uni - OpenCL error - Buffer Y allocation failed: " ;
    the_err += err.what();
    the_err += "\n";
    throw std::runtime_error(the_err);
  }

  // Write buffer

  // Image
  m_buflen_Im  = m_No*sizeof(T);
  try {
    m_cl_buf_Im = cl::Buffer(context, CL_MEM_WRITE_ONLY, m_buflen_Im);
  }

  catch (const cl::Error &err) {
    std::string the_err = "Error in das_uni - OpenCL error - Buffer Im allocation failed:" ;
    the_err += err.what();
    the_err += "\n";
    throw std::runtime_error(the_err);
  }

  }

  ~DAS_UNI()  = default;

  // Function to check if we need to resize buffers and re-init.
  bool das_setup_has_changed(DASType das_type,
                             dream_idx_type a_scan_len,
                             T min_t, T pitch_t, T max_t, // Transmit
                             T min_r, T pitch_r, T max_r, // Receive
                             T min_Rx, T dx, T max_Rx,    // Observation points
                             T min_Ry, T dy, T max_Ry,
                             T min_Rz, T dz, T max_Rz) {
    bool retval = false;

    // Array dims
    dream_idx_type num_t_elements = (dream_idx_type) ((max_t - min_t)/pitch_t+1.0);
    dream_idx_type num_r_elements = (dream_idx_type) ((max_r - min_r)/pitch_r+1.0);

    // Image dims
    dream_idx_type Nx = (dream_idx_type) ((max_Rx - min_Rx)/dx+1.0);
    dream_idx_type Ny = (dream_idx_type) ((max_Ry - min_Ry)/dy+1.0);
    dream_idx_type Nz = (dream_idx_type) ((max_Rz - min_Rz)/dz+1.0);
    dream_idx_type No = Nx*Ny*Nz;

    if ( das_type != m_das_type ||
         a_scan_len != m_a_scan_len ||
         No != m_No ||
         num_t_elements != m_num_t_elements ||
         num_r_elements != m_num_r_elements) {
      retval = true;
    }

    return retval;
  }

  int cl_das_uni(const T *Y, // Data
                 T min_t, T pitch_t, T max_t, // Transmit
                 T min_r, T pitch_r, T max_r, // Receive
                 T min_Rx, T dx, T max_Rx,    // Observation points
                 T min_Ry, T dy, T max_Ry,
                 T min_Rz, T dz, T max_Rz,
                 T dt,
                 T delay,
                 T cp,
                 T *Im);

 private:

  DASType m_das_type;
  dream_idx_type m_a_scan_len;
  dream_idx_type m_No;
  dream_idx_type m_num_t_elements;
  dream_idx_type m_num_r_elements;
  ErrorLevel m_out_err;

  cl::CommandQueue m_queue;
  cl::Kernel m_kernel;
  cl::Buffer m_cl_buf_Y;
  cl::Buffer m_cl_buf_Im;
  cl_int m_buflen_Y;
  cl_int m_buflen_Im;
};
