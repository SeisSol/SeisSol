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

namespace seissol {
class AcceleratorDevice {
  public:
  static AcceleratorDevice& getInstance() {
    static AcceleratorDevice instance;
    return instance;
  }

  void bindAcceleratorDevice(int deviceId) {
    bindNativeDevice(deviceId);
  }

  private:
  void bindNativeDevice(int mpiRank);

  void printInfo();

  private:
  void bindNativeDevice(int deviceId);
  void bindSyclDevice(int deviceId);

  std::vector<std::string> infoMessages;
  std::vector<std::string> warnMessages;
};
} // namespace seissol

#endif
