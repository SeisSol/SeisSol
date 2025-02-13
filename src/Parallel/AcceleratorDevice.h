// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_ACCELERATORDEVICE_H_
#define SEISSOL_SRC_PARALLEL_ACCELERATORDEVICE_H_

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

  sycl::queue getInorderSyclQueue();

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

#endif // SEISSOL_SRC_PARALLEL_ACCELERATORDEVICE_H_
