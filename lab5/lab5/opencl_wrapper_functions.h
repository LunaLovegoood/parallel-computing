#ifndef LAB5_OPENCL_WRAPPER_FUNCTIONS_H_
#define LAB5_OPENCL_WRAPPER_FUNCTIONS_H_

#include "CL/cl.hpp"

#include <string>
#include <vector>

bool CheckOpenCLSupport(const cl_device_type device_type);
cl::Program CreateProgram(
    const std::string &kernel_path,
    const std::string &build_options
);
cl::Device GetFirstGPUDevice();

#endif  // LAB5_OPENCL_WRAPPER_FUNCTIONS_H_
