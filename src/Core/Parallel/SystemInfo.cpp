// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "SystemInfo.h"

#include "Core/Parallel/MPI.h"

#include <mpi.h>
#include <string>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol {

void SystemInfo::init() {
  loadHostInfo();
  // no dedicated CPU info printing yet (only GPU)
  loadGPUInfo();
}

void SystemInfo::loadHostInfo() {
  // (gethostname is POSIX; on Windows, it is part of Winsock and needs extra initialization)
  std::string hostName(MPI_MAX_PROCESSOR_NAME, '\0');
  int hostNameLength = 0;
  if (MPI_Get_processor_name(hostName.data(), &hostNameLength) != MPI_SUCCESS) {
    hostName = "unknown-host";
  } else {
    hostName.resize(hostNameLength);
  }
  hostNames_ = Mpi::mpi.collectContainer(hostName);

  // the names are gathered on rank 0 only
  if (Mpi::mpi.rank() == 0) {
    logInfo() << "Running on (rank=0):" << hostNames_.front();
  }
}

void SystemInfo::loadGPUInfo() {
#ifdef ACL_DEVICE
  const device::DeviceInstance& device = device::DeviceInstance::instance();
  const auto deviceId = device.api().getDeviceId();
  const auto pci = device.api().getPciAddress(deviceId);
  const auto pcisNode = Mpi::mpi.collectContainer(pci, Mpi::mpi.sharedMemComm());
  pcis_ = Mpi::mpi.collectContainer(pci);

  logInfo() << "Device API:" << device.api().getApiName();
  logInfo() << "Device name (rank=0):" << device.api().getDeviceName(deviceId);
  logInfo() << "Device PCI address (rank=0): " << pci;
  logInfo() << "Device PCI addresses (node of rank=0):" << pcisNode;
#endif
}

} // namespace seissol
