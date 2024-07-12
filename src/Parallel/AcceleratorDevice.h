/**
 * @file
 * This file is part of SeisSol.
 *
 * The files which include this header should be added to
 * `general-sycl-offloading` CMake target
 */

#ifndef ACCELERATOR_DEVICE_H
#define ACCELERATOR_DEVICE_H

#include <CL/sycl.hpp>
#include <device.h>
#include <string>
#include <vector>

#ifndef __DPCPP_COMPILER
namespace sycl = cl::sycl;
#endif

namespace seissol {
class AcceleratorDevice {
  public:
  static AcceleratorDevice& getInstance() {
    static AcceleratorDevice instance;
    return instance;
  }

  void bindAcceleratorDevice(int deviceId) {
    bindSyclDevice(deviceId);
    bindNativeDevice(deviceId);
  }

  sycl::device& getSyclDevice() { return syclDevice; }

  sycl::queue& getSyclDefaultQueue() { return syclDefaultQueue; }

  void printInfo();

  private:
  void bindNativeDevice(int deviceId);
  void bindSyclDevice(int deviceId);

  sycl::device syclDevice;
  sycl::queue syclDefaultQueue;

  std::vector<std::string> infoMessages;
  std::vector<std::string> warnMessages;
};
} // namespace seissol

#endif
