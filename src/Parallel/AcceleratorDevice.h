// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
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

  void printInfo();

  private:
  void bindNativeDevice(int deviceId);

  std::vector<std::string> infoMessages;
  std::vector<std::string> warnMessages;
};
} // namespace seissol

#endif // SEISSOL_SRC_PARALLEL_ACCELERATORDEVICE_H_
