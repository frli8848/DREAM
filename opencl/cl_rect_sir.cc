/***
*
* Copyright (C) 2020 Fredrik Lingvall
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


#include <math.h>
#include <iostream>
#include <fstream>
#include <streambuf>

#include "rect_sir.h"

#define CL_TARGET_OPENCL_VERSION 120
#include <CL/cl.h>

/***
 *
 * OpenCL version of rect_sir
 *
 ***/

int RectSir::cl_rect_sir(const double *Ro,
                         int No,
                         double a,
                         double b,
                         double dt,
                         int nt,
                         double delay,
                         double v,
                         double cp,
                         double *h)
{
  std::string kernel_str;

    std::string dream_cl_kernel = "";
  if(const char* env_p = std::getenv("DREAM_CL_KERNELS")) {
    dream_cl_kernel += env_p;
  } else {
    throw std::runtime_error("Error in rect_sir - DREAM_CL_KERNELS env variable not set!");
  }

  std::string rect_sir_kernel = dream_cl_kernel;
  rect_sir_kernel += "/rect_sir.cl";
  std::ifstream f_kernel(rect_sir_kernel);
  f_kernel.seekg(0, std::ios::end);
  kernel_str.reserve(f_kernel.tellg());
  f_kernel.seekg(0, std::ios::beg);
  kernel_str.assign((std::istreambuf_iterator<char>(f_kernel)),
             std::istreambuf_iterator<char>());

  //std::cout << "CL Kernel: " <<  kernel_str.c_str() << std::endl;

  size_t source_size = kernel_str.length();

  //
  // Get platform and device information
  //

  cl_platform_id platform_id = NULL;
  cl_device_id device_id = NULL;
  cl_uint ret_num_devices;
  cl_uint ret_num_platforms;

  cl_int ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);

  //ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, 1, &device_id, &ret_num_devices);
  ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &ret_num_devices);
  if (ret < 0) std::cout << "ret: " << ret << " num devices: " << ret_num_devices << std::endl;

  // Create an OpenCL context
  cl_context context = clCreateContext( NULL, 1, &device_id, NULL, NULL, &ret);
  if (ret < 0) std::cout << "cl context ret: " << ret << std::endl;

  // Create a command queue
  cl_command_queue command_queue = clCreateCommandQueue(context, device_id, 0, &ret);
  if (ret < 0) std::cout << "cl command queue ret: " << ret << std::endl;

  //
  // Create memory buffers on the device for arguement to the  kernel
  //

  // Read

  cl_mem Ro_mo = clCreateBuffer(context, CL_MEM_READ_ONLY, No*3*sizeof(double), NULL, &ret);
  if (ret < 0) std::cout << "cl Ro mem ret: " << ret << std::endl;
  cl_int No_mo = No;

  cl_double a_mo = a;
  cl_double b_mo = b;
  cl_double dt_mo = dt;
  cl_int nt_mo = nt;
  cl_double delay_mo = delay;
  cl_double v_mo = v;
  cl_double cp_mo = cp;

  // Write

  cl_mem H_mo = clCreateBuffer(context, CL_MEM_WRITE_ONLY, nt * No * sizeof(double), NULL, &ret);
  if (ret < 0) std::cout << "cl H mem ret: " << ret << std::endl;

  //
  // Copy input arguments to their respective memory buffers
  //

  ret = clEnqueueWriteBuffer(command_queue, Ro_mo, CL_TRUE, 0, No*3 * sizeof(double), Ro,  0, NULL, NULL);

  //
  // Create a program from the kernel source
  //

  const char *k_c_str = kernel_str.c_str();
  cl_program program = clCreateProgramWithSource(context, 1, (const char **) &k_c_str, (const size_t *) &source_size, &ret);
  if (ret < 0) std::cout << "cl program ret: " << ret << std::endl;

  //
  // Build the program
  //

  ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
  if (ret < 0) std::cout << "cl build program ret: " << ret << std::endl;

  if (ret == CL_BUILD_PROGRAM_FAILURE) {
    // Determine the size of the log
    size_t log_size;
    clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);

    // Allocate memory for the log
    char *log = (char *) malloc(log_size);

    // Get the log
    clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);

    // Print the log
    printf("%s\n", log);

    free(log);
  }

  // Create the OpenCL kernel
  cl_kernel kernel = clCreateKernel(program, "rect_sir", &ret);
  if (ret < 0) std::cout << "cl kernel ret: " << ret << std::endl;

  // Set the arguments of the kernel

  ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &Ro_mo);
  if (ret < 0) std::cout << "cl kernel arg Ro ret: " << ret << std::endl;

  ret = clSetKernelArg(kernel, 1, sizeof(cl_int), (void *) &No_mo);
  if (ret < 0) std::cout << "cl kernel arg No ret: " << ret << std::endl;

  ret = clSetKernelArg(kernel, 2, sizeof(cl_double), (void *) &a_mo);
  if (ret < 0) std::cout << "cl kernel arg a ret: " << ret << std::endl;
  ret = clSetKernelArg(kernel, 3, sizeof(cl_double), (void *) &b_mo);
  if (ret < 0) std::cout << "cl kernel arg b ret: " << ret << std::endl;

  ret = clSetKernelArg(kernel, 4, sizeof(cl_double), (void *) &dt_mo);
  if (ret < 0) std::cout << "cl kernel arg dt ret: " << ret << std::endl;

  ret = clSetKernelArg(kernel, 5, sizeof(cl_int), (void *) &nt_mo);
  if (ret < 0) std::cout << "cl kernel arg nt ret: " << ret << std::endl;

  ret = clSetKernelArg(kernel, 6, sizeof(cl_double), (void *) &delay_mo);
  if (ret < 0) std::cout << "cl kernel arg delay ret: " << ret << std::endl;
  ret = clSetKernelArg(kernel, 7, sizeof(cl_double), (void *) &v_mo);
  if (ret < 0) std::cout << "cl kernel arg v ret: " << ret << std::endl;
  ret = clSetKernelArg(kernel, 8, sizeof(cl_double), (void *) &cp_mo);
  if (ret < 0) std::cout << "cl kernel arg cp ret: " << ret << std::endl;

  ret = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &H_mo);
  if (ret < 0) std::cout << "cl kernel arg H ret: " << ret << std::endl;

  //
  // Execute the OpenCL kernel
  //

  size_t global_item_size = No; // Process the entire number of observation points
  size_t local_item_size = 64; // NB No must be dividable by the local item  size (work group size).
  ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL);
  if (ret < 0) std::cout << "cl enqueue ret: " << ret << std::endl;

  //
  // Read the memory buffer C on the device to the local variable C
  //

  ret = clEnqueueReadBuffer(command_queue, H_mo, CL_TRUE, 0, nt*No*sizeof(double), h, 0, NULL, NULL);
  if (ret < 0) std::cout << "cl read ret: " << ret << std::endl;

  //
  // Clean up
  //

  ret = clFlush(command_queue);
  ret = clFinish(command_queue);
  ret = clReleaseKernel(kernel);
  ret = clReleaseProgram(program);

  ret = clReleaseMemObject(Ro_mo);
  ret = clReleaseMemObject(H_mo);

  ret = clReleaseCommandQueue(command_queue);
  ret = clReleaseContext(context);

  int err = ret;

  return err;
}
