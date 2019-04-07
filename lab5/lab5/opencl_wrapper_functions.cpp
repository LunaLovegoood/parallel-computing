#define CL_USE_DEPRECATED_OPENCL_2_0_APIS

#include "opencl_wrapper_functions.h"

#include <fstream>
#include <iostream>

namespace {

std::vector<cl::Platform> GetPlatforms();
std::vector<cl::Device> GetDevices(
    const cl::Platform &platform,
    const cl_device_type device_type
);
std::string ReadKernelSource(const std::string &kernel_path);
void FatalBuildError(cl_int err, cl::Program program);

} // namespace

bool CheckOpenCLSupport(const cl_device_type device_type) {
  // Check if there are available platforms
  auto platforms = GetPlatforms();
  if (platforms.size() == 0) {
    return false;
  }

  // Check if there are available devices
  auto devices = GetDevices(platforms.front(), device_type);
  if (devices.size() == 0) {
    return false;
  }

  return true;
}

cl::Program CreateProgram(
    const std::string &kernel_path,
    const std::string &build_options
) {
  std::string kernel_source = ReadKernelSource(kernel_path);
  auto device = GetFirstGPUDevice();

  cl::Program::Sources sources(1,
    std::make_pair(kernel_source.c_str(), kernel_source.length() + 1));

  cl::Context context(device);
  cl::Program program(context, sources);

  cl_int err = program.build(build_options.c_str());
  if (err) {
    FatalBuildError(err, program);
  }

  return program;
}

cl::Device GetFirstGPUDevice() {
  return GetDevices(GetPlatforms().front(), CL_DEVICE_TYPE_GPU).front();
}

namespace {

std::vector<cl::Platform> GetPlatforms() {
  std::vector<cl::Platform> platforms;
  cl::Platform::get(&platforms);
  return platforms;
}

std::vector<cl::Device> GetDevices(
    const cl::Platform &platform,
    const cl_device_type device_type
) {
  std::vector<cl::Device> devices;
  platform.getDevices(device_type, &devices);
  return devices;
}

std::string ReadKernelSource(const std::string &kernel_path) {
  std::ifstream kernel_file(kernel_path);
  return std::string(std::istreambuf_iterator<char>(kernel_file),
                     (std::istreambuf_iterator<char>()));
}

void FatalBuildError(cl_int err, cl::Program program) {
  auto context = program.getInfo<CL_PROGRAM_CONTEXT>();
  auto device = context.getInfo<CL_CONTEXT_DEVICES>().front();

  std::string err_log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
  std::cerr << err_log << std::endl;

  system("pause");
  exit(1);
}

} // namespace
