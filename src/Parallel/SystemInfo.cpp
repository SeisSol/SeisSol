// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "SystemInfo.h"

#include "Parallel/MPI.h"

#include <string>
#include <unistd.h>
#include <utils/logger.h>
#include <utils/stringutils.h>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol {

void SystemInfo::init() {
  loadHostInfo();
  loadCPUInfo();
  loadGPUInfo();
}

void SystemInfo::loadHostInfo() {
  std::string hostName(256, ' ');
  if (gethostname(hostName.data(), 256) != 0) {
    hostName = "unknown-host";
  } else {
    utils::StringUtils::rtrim(hostName);
    hostName.pop_back();
  }
  hostNames = Mpi::mpi.collectContainer(hostName);

  logInfo() << "Running on (rank=0):" << hostNames.front();
}

void SystemInfo::loadCPUInfo() {
  // do nothing. Yet.
}

void SystemInfo::loadGPUInfo() {
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  const auto pci = device.api->getPciAddress(0);
  const auto pcisNode = Mpi::mpi.collectContainer(pci, Mpi::mpi.sharedMemComm());
  pcis = Mpi::mpi.collectContainer(pci);

  logInfo() << "Device API:" << device.api->getApiName();
  logInfo() << "Device name (rank=0):" << device.api->getDeviceName(0);
  logInfo() << "Device PCI address (rank=0): " << pci;
  logInfo() << "Device PCI addresses (node of rank=0):" << pcisNode;
#endif
}

} // namespace seissol
