/***
 *
 * Copyright (C) 2003,2006,2007,2008,2009,2021,2023,2024 Fredrik Lingvall
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

#include <iostream>
#include <thread>

#include "dream.h"
#include "dream_error.h"

#ifdef USE_OPENCL

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
#endif

// This switches on das arg parsing
// in the arg_parser_xx:ers
#define DELAY_AND_SUM

template <class T>
class DAS
{
public:

  DAS(DASType das_type,
      dream_idx_type a_scan_len,
      dream_idx_type No,
      dream_idx_type num_t_elements,
      dream_idx_type num_r_elements, // SAFT if num_r_elements = 0;
      bool use_gpu=false);

  ~DAS()  = default;

  SIRError das(const T *Y, const T *Ro, const T *gt, const T *gr,
               T dt,
               DelayType delay_type, const T *delay,
               T cp,
               T *Im,
               ErrorLevel err_level);

  static void abort(int signum);
  bool is_running();
  void set_running();
  bool das_setup_has_changed(DASType das_type, // Function to check if we need to resize buffers and re-init.
                             dream_idx_type a_scan_len,
                             dream_idx_type No,
                             dream_idx_type num_t_elements,
                             dream_idx_type num_r_elements,
                             bool use_gpu) {
    bool retval = false;
    if ( das_type != m_das_type ||
         a_scan_len != m_a_scan_len ||
         No != m_No ||
         num_t_elements != m_num_t_elements ||
         num_r_elements != m_num_r_elements ||
         m_use_gpu != use_gpu) {
      retval = true;
    }

    return retval;
  }

#ifdef USE_OPENCL

  void init_opencl(unsigned char *cl_kernel_str, unsigned int cl_kernel_str_len)
  {
    std::string kernel_str;
    std::string kernel_name="das_tfm"; // Default function name of the OpenCL kernel

    if (m_das_type == DASType::saft) {
      kernel_name="das_saft";
      m_num_r_elements = 0;
    }

    if (m_das_type == DASType::tfm) {
      kernel_name="das_tfm";
    }

    if (m_das_type == DASType::rca_coltx) {
      kernel_name="das_rca_coltx";
    }

    if (m_das_type == DASType::rca_rowtx) {
      kernel_name="das_rca_rowtx";
    }

    //
    // Read OpenCL kernels
    //

    if (cl_kernel_str_len == 0)  { // Fallback to read kernels from file when no xxd:ed kernel header files exist.
      std::string dream_cl_kernel = "";
      if(const char* env_p = std::getenv("DREAM_CL_KERNELS")) {
        dream_cl_kernel += env_p;
      } else {
        throw std::runtime_error("Error in das - DREAM_CL_KERNELS env variable not set!");
      }

      std::string das_kernel = dream_cl_kernel;
      if (sizeof(T) == sizeof(float) ) {
        das_kernel += "/das_float.cl";
      } else {
        das_kernel += "/das_double.cl";
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
      std::string the_err = "Error in das - OpenCL error: " ;
      the_err += err.what();
      the_err += "\n";
      throw std::runtime_error(the_err);
    }

    if (platforms.empty()) {
      throw std::runtime_error("Error in das - OpenCL error: No platforms found!");
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
      std::string the_err = "Error in das - OpenCL error: ";
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
      std::string the_err = "Error in das - Failed to get OpenCL context info: ";
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
      throw std::runtime_error("Error in das - OpenCL error: Found no devices(s)!");
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
      std::string the_err = "Error in das - OpenCL error: ";
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
      std::string the_err = "Error in das - OpenCL error - kernel build failed: ";
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
      std::cerr << "Error in das - OpenCL error - Kernel: " << kernel_name << " creation failed: " << err.what();
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
      std::string the_err = "Error in das - OpenCL error - Buffer Y allocation failed: " ;
      the_err += err.what();
      the_err += "\n";
      throw std::runtime_error(the_err);
    }

    // Transmit elements
    m_buflen_gt = m_num_t_elements*3*sizeof(T);
    try {
      m_cl_buf_gt = cl::Buffer(context, CL_MEM_READ_ONLY, m_buflen_gt);
    }

    catch (const cl::Error &err) {
      std::string the_err = "Error in das - OpenCL error - Buffer gt allocation failed: " ;
      the_err += err.what();
      the_err += "\n";
      throw std::runtime_error(the_err);
    }

    // Receive elements
    if (m_num_r_elements > 0) { // SAFT do not use this one.

      m_buflen_gr = m_num_r_elements*3*sizeof(T);
      try {
        m_cl_buf_gr = cl::Buffer(context, CL_MEM_READ_ONLY, m_buflen_gr);
      }

      catch (const cl::Error &err) {
        std::string the_err = "Error in das - OpenCL error - Buffer gr allocation failed: " ;
        the_err += err.what();
        the_err += "\n";
        throw std::runtime_error(the_err);
      }

    }

    // Observation points
    m_buflen_Ro  =  m_No*3*sizeof(T);
    try {
      m_cl_buf_Ro = cl::Buffer(context, CL_MEM_READ_ONLY, m_buflen_Ro);
    }

    catch (const cl::Error &err) {
      std::string the_err = "Error in das - OpenCL error - Buffer Ro allocation failed: " ;
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
      std::string the_err = "Error in das - OpenCL error - Buffer Im allocation failed: " ;
      the_err += err.what();
      the_err += "\n";
      throw std::runtime_error(the_err);
    }
  }

  int cl_das(const T *Y, // Data
             const T *Ro,
             const T *gt,
             const T *gr,
             T dt,
             T delay,
             T cp,
             T *Im);
#endif

private:

  void* smp_das(void *arg);
  std::thread das_thread(void *arg) {
    return std::thread(&DAS::smp_das, this, arg);
  }

  SIRError das_saft_serial(const T *Y, // Size: nt x num_elements
                           dream_idx_type a_scan_len,
                           const T *g, dream_idx_type num_elements,
                           T xo, T yo, T zo,
                           T dt, T delay,
                           T cp, T &im, ErrorLevel err_level);

  SIRError das_tfm_serial(const T *Y, // Size: nt x num_t_elements*num_r_elements (=FMC)
                          dream_idx_type a_scan_len,
                          const T *gt, dream_idx_type num_t_elements,
                          const T *gr, dream_idx_type num_r_elements,
                          T xo, T yo, T zo,
                          T dt, T delay,
                          T cp, T &im, ErrorLevel err_level);

  SIRError das_rca_serial_coltx(const T *Y, // Size: a_scan_len x num_cols*num_rows (=FMC)
                                dream_idx_type a_scan_len,
                                const T *G_col, dream_idx_type num_cols,
                                const T *G_row, dream_idx_type num_rows,
                                T xo, T yo, T zo,
                                T dt, T delay,
                                T cp, T &im, ErrorLevel err_level);

  SIRError das_rca_serial_rowtx(const T *Y, // Size: a_scan_len x num_cols*num_rows (=FMC)
                                dream_idx_type a_scan_len,
                                const T *G_col, dream_idx_type num_cols,
                                const T *G_row, dream_idx_type num_rows,
                                T xo, T yo, T zo,
                                T dt, T delay,
                                T cp, T &im, ErrorLevel err_level);

  DASType m_das_type;
  dream_idx_type m_a_scan_len;
  dream_idx_type m_No;
  dream_idx_type m_num_t_elements;
  dream_idx_type m_num_r_elements;
  ErrorLevel m_out_err;
  bool m_use_gpu;

#ifdef USE_OPENCL
  cl::CommandQueue m_queue;
  cl::Kernel m_kernel;
  cl::Buffer m_cl_buf_Y;
  cl::Buffer m_cl_buf_gt;
  cl::Buffer m_cl_buf_gr;
  cl::Buffer m_cl_buf_Ro;
  cl::Buffer m_cl_buf_Im;
  cl_int m_buflen_Y;
  cl_int m_buflen_gt;
  cl_int m_buflen_gr;
  cl_int m_buflen_Ro;
  cl_int m_buflen_Im;
#endif

};
