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

#include <iostream>
#include <fstream>
#include <streambuf>
// #include <cassert>

#include <list>
#include <vector>

#include "das.h"

#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_CL_1_2_DEFAULT_BUILD
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY
#include <CL/opencl.hpp>

/***
 *
 * OpenCL version of das
 *
 ***/

int DAS::cl_das(const double *Y, // Data
                dream_idx_type a_scan_len,
                const double *Ro, dream_idx_type No,
                const double *gt, dream_idx_type num_t_elements,
                const double *gr, dream_idx_type num_r_elements,
                double dt,
                double delay,
                double cp,
                DASType das_type,
                double *Im)
{
  std::string kernel_str;
  std::string kernel_name="das_tfm"; // Default function name of the OpenCL kernel

  if (das_type == DASType::tfm) {
    kernel_name="das_tfm";
  }

  if (das_type == DASType::rca) {
    kernel_name="das_rca";
  }

  //
  // Read OpenCL kernel file
  //

  std::string dream_cl_kernel = "";
  if(const char* env_p = std::getenv("DREAM_CL_KERNELS")) {
    dream_cl_kernel += env_p;
  } else {
    throw std::runtime_error("Error in das - DREAM_CL_KERNELS env variable not set!");
  }

  std::string das_kernel = dream_cl_kernel;
  das_kernel += "/das.cl";
  std::ifstream f_kernel(das_kernel);
  f_kernel.seekg(0, std::ios::end);
  kernel_str.reserve(f_kernel.tellg());
  f_kernel.seekg(0, std::ios::beg);
  kernel_str.assign((std::istreambuf_iterator<char>(f_kernel)),
                    std::istreambuf_iterator<char>());

  //std::cout << "CL Kernel: " <<  kernel_str.c_str() << std::endl;

  //
  // Get platform
  //

  std::vector<cl::Platform> platforms;
  cl::Platform::get(&platforms);

  if (platforms.empty()) {
    std::cerr << "No platforms found!" << std::endl;;
    return -1;
    //trow std::runtime_error(Platform size 0);
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
    std::cerr << "OpenCL error: " << err.what();
  }

  //
  // Get a list of devices on this platform
  //

  std::vector<cl::Device> devices;

  try {
    devices = context.getInfo<CL_CONTEXT_DEVICES>();
  }

  catch (const cl::Error &err) {
    std::cerr << "Failed to get OpenCL context info: " << err.what();
  }


  size_t dev_idx = 0;
  /*
  if (dev_idx >= devices.size()) {
   std::cerr << "Selected OpenCL device not found! Using the first device!";
    dev_idx = 0;
  }
  */

  if (devices.empty()) {
   std::cerr << "Found no devices(s)!";
    return -13;
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

  cl::CommandQueue queue;

  try {
    queue = cl::CommandQueue(context, devices[dev_idx]);
  }

  catch (const cl::Error &err) {
    std::cerr << "OpenCL error: " << err.what();
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
    std::cerr << "Kernel build failed: " << err.what() << std::endl;
    std::cout << "Build log:" << std::endl << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[dev_idx]) << std::endl;
  }

  // Make kernel
  cl::Kernel kernel;

  try {
    kernel = cl::Kernel(program, kernel_name.c_str());
  }

  catch (const cl::Error &err) {
    std::cerr << "Kernel: " << kernel_name << " creation failed: " << err.what();
  }

  //
  // Create memory buffers on the device for vector/matrix arguments to the kernel
  //

  // Read buffers

  // Data
  cl_int buflen_Y  =  a_scan_len*num_t_elements*num_r_elements*sizeof(double);
  cl::Buffer cl_buf_Y;
  try {
    cl_buf_Y = cl::Buffer(context, CL_MEM_READ_ONLY, buflen_Y);
  }

  catch (const cl::Error &err) {
    std::cerr << "Buffer Y allocation failed: " << err.what();
  }

  // Transmit elements
  cl_int buflen_gt = num_t_elements*3*sizeof(double);
  cl::Buffer cl_buf_gt = cl::Buffer(context, CL_MEM_READ_ONLY, buflen_gt);

  // Receive elements
  cl_int buflen_gr = num_r_elements*3*sizeof(double);
  cl::Buffer cl_buf_gr = cl::Buffer(context, CL_MEM_READ_ONLY, buflen_gr);

  // Observation points
  cl_int buflen_Ro  =  No*3*sizeof(double);
  cl::Buffer cl_buf_Ro = cl::Buffer(context, CL_MEM_READ_ONLY, buflen_Ro);

  // Write buffer

  // Image
  cl_int buflen_Im  = No*sizeof(double);
  cl::Buffer cl_buf_Im = cl::Buffer(context, CL_MEM_WRITE_ONLY, buflen_Im);

  //
  // Transfer input vector/matrix args to their respective memory buffers
  //

  queue.enqueueWriteBuffer(cl_buf_Y, CL_TRUE, 0, buflen_Y, Y);
  queue.enqueueWriteBuffer(cl_buf_gt, CL_TRUE, 0, buflen_gt, gt);
  queue.enqueueWriteBuffer(cl_buf_gr, CL_TRUE, 0, buflen_gr, gr);
  queue.enqueueWriteBuffer(cl_buf_Ro, CL_TRUE, 0, buflen_Ro, Ro);

  //
  // Set the arguments of the kernel
  //

  size_t arg_idx = 0;

  // Y
  kernel.setArg(arg_idx, cl_buf_Y);
  arg_idx++;
  kernel.setArg(arg_idx, (int) a_scan_len);

  // gt
  arg_idx++;
  kernel.setArg(arg_idx, cl_buf_gt);
  arg_idx++;
  kernel.setArg(arg_idx, (int) num_t_elements);

  // gr
  arg_idx++;
  kernel.setArg(arg_idx, cl_buf_gr);
  arg_idx++;
  kernel.setArg(arg_idx, (int) num_r_elements);

  // Ro
  arg_idx++;
  kernel.setArg(arg_idx, cl_buf_Ro);
  arg_idx++;
  kernel.setArg(arg_idx, (int) No);

  // dt
  arg_idx++;
  kernel.setArg(arg_idx, dt);

  // delay
  arg_idx++;
  kernel.setArg(arg_idx, delay);

  // cp
  arg_idx++;
  kernel.setArg(arg_idx, cp);

  // Im (output arg)
  arg_idx++;
  kernel.setArg(arg_idx, cl_buf_Im);

  //
  // Execute the OpenCL kernel
  //

  cl::NDRange global(No);
  cl::NDRange local(64); // Works on both Nvidia and AMD GPUs.
  try {
    queue.enqueueNDRangeKernel(kernel, cl::NullRange, global, local);
  }

  catch (const cl::Error &err) {
    std::cerr << "OpenCL kernel run failed: " << err.what() << std::endl;
  }

  //
  // Transfer the memory buffer on the device to the host buffer
  //

  try {
    queue.enqueueReadBuffer(cl_buf_Im, CL_TRUE, 0,  buflen_Im, Im);
  }

  catch (const cl::Error &err) {
    std::cerr << "Failed to get results buffer from GPU: " << err.what() << std::endl;
  }

  int retval = 0;

  return retval;
}
