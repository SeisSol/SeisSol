// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Parallel/AcceleratorDevice.h"
#include "utils/logger.h"
#include <sstream>
#include <string>

#include <device.h>

namespace seissol {

void AcceleratorDevice::bindNativeDevice(int deviceId) {
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  {
    std::ostringstream info;
    info << "Device API: " << device.api->getApiName();
    infoMessages.push_back(info.str());
  }
  {
    std::ostringstream info;
    info << "Device (rank=0): " << device.api->getDeviceName(deviceId);
    infoMessages.push_back(info.str());
  }
  device.api->setDevice(deviceId);
}

void AcceleratorDevice::printInfo() {
  for (const auto& warn : warnMessages) {
    logWarning() << warn.c_str();
  }
  for (const auto& info : infoMessages) {
    logInfo() << info.c_str();
  }
}
} // namespace seissol
