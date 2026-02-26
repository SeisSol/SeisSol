// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Parallel/AcceleratorDevice.h"

#include <string>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include <Device/device.h>
#include <sstream>
#endif

namespace seissol {

void AcceleratorDevice::bindNativeDevice(int deviceId) {
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  {
    std::ostringstream info;
    info << "Device API: " << device.api->getApiName();
    infoMessages_.push_back(info.str());
  }
  {
    std::ostringstream info;
    info << "Device (rank=0): " << device.api->getDeviceName(deviceId);
    infoMessages_.push_back(info.str());
  }
  device.api->setDevice(deviceId);
#endif
}

void AcceleratorDevice::printInfo() {
  for (const auto& warn : warnMessages_) {
    logWarning() << warn.c_str();
  }
  for (const auto& info : infoMessages_) {
    logInfo() << info.c_str();
  }
}
} // namespace seissol
