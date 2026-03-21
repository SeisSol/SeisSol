// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#include "MPI.h"

#include <algorithm>
#include <cctype>
#include <mpi.h>
#include <string>
#include <unistd.h>
#include <utils/env.h>
#include <utils/logger.h>
#include <utils/stringutils.h>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

void seissol::Mpi::init(int& argc, char**& argv) {
  // Note: Strictly speaking, we only require MPI_THREAD_MULTIPLE if using
  // a communication thread and/or async I/O.
  // The safer (and more sane) option is to enable it by default.
  const int required = MPI_THREAD_MULTIPLE;
  int provided = 0;
  MPI_Init_thread(&argc, &argv, required, &provided);

  setComm(MPI_COMM_WORLD);

  std::string hostName(256, ' ');
  if (gethostname(hostName.data(), 256) != 0) {
    hostName = "unknown-host";
  } else {
    utils::StringUtils::rtrim(hostName);
    hostName.pop_back();
  }
  hostNames = collectContainer(hostName);

  // Test this after setComm() to get the correct rank_
  if (provided < required) {
    logError() << utils::nospace << "Provided MPI thread support (" << provided
               << ") is less than required thread support (" << required << ").";
  }
}

void seissol::Mpi::setComm(MPI_Comm comm) {
  comm_ = comm;

  MPI_Comm_rank(comm, &rank_);
  MPI_Comm_size(comm, &size_);

  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &sharedMemComm_);
  MPI_Comm_rank(sharedMemComm_, &sharedMemMpiRank_);
  MPI_Comm_size(sharedMemComm_, &sharedMemMpiSize_);
}

void seissol::Mpi::bindAcceleratorDevice() {
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  device.api->setDevice(0);
#endif
}

void seissol::Mpi::printAcceleratorDeviceInfo() {
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  const auto pci = device.api->getPciAddress(0);
  const auto pcisNode = collectContainer(pci, sharedMemComm_);
  pcis = collectContainer(pci);
  logInfo() << "Device API:" << device.api->getApiName();
  logInfo() << "Device name (rank=0):" << device.api->getDeviceName(0);
  logInfo() << "Device PCI address (rank=0): " << pci;
  logInfo() << "Device PCI addresses (node of rank=0):" << pcisNode;
#endif
}

seissol::Mpi seissol::Mpi::mpi;
