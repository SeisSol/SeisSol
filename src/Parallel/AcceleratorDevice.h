/**
 * @file
 * This file is part of SeisSol.
 *
 * The files which include this header should be added to
 * `general-sycl-offloading` CMake target
 */

#ifndef ACCELERATOR_DEVICE_H
#define ACCELERATOR_DEVICE_H

#include <device.h>
#include <CL/sycl.hpp>

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

  private:
  void bindNativeDevice(int mpiRank);
  void bindSyclDevice(int mpiRank);

  sycl::device syclDevice;
  sycl::queue syclDefaultQueue;
};
} // namespace seissol

#endif
