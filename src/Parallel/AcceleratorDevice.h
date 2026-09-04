// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_ACCELERATORDEVICE_H_
#define SEISSOL_SRC_PARALLEL_ACCELERATORDEVICE_H_

#include <string>
#include <vector>

namespace seissol {
class AcceleratorDevice {
  public:
  static AcceleratorDevice& getInstance() {
    static AcceleratorDevice instance;
    return instance;
  }

  void bindAcceleratorDevice(int deviceId) { bindNativeDevice(deviceId); }

  /**
   * The device this rank is bound to. The current device is thread-local in both CUDA and HIP,
   * so every thread that issues device work has to re-select it.
   */
  [[nodiscard]] int getDeviceId() const { return deviceId_; }

  void printInfo();

  private:
  void bindNativeDevice(int deviceId);

  int deviceId_{0};
  std::vector<std::string> infoMessages_;
  std::vector<std::string> warnMessages_;
};
} // namespace seissol

#endif // SEISSOL_SRC_PARALLEL_ACCELERATORDEVICE_H_
