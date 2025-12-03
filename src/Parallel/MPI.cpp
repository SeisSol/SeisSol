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
#include "Parallel/AcceleratorDevice.h"

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
  if (gethostname(const_cast<char*>(hostName.c_str()), 256) != 0) {
    hostName = "unknown-host";
  } else {
    utils::StringUtils::rtrim(hostName);
    hostName.pop_back();
  }
  hostNames = collectContainer(hostName);

  // Test this after setComm() to get the correct m_rank
  if (provided < required) {
    logError() << utils::nospace << "Provided MPI thread support (" << provided
               << ") is smaller than required thread support (" << required << ").";
  }
}

void seissol::Mpi::setComm(MPI_Comm comm) {
  m_comm = comm;

  MPI_Comm_rank(comm, &m_rank);
  MPI_Comm_size(comm, &m_size);

  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &m_sharedMemComm);
  MPI_Comm_rank(m_sharedMemComm, &m_sharedMemMpiRank);
  MPI_Comm_size(m_sharedMemComm, &m_sharedMemMpiSize);
}

void seissol::Mpi::bindAcceleratorDevice() {
#ifdef ACL_DEVICE
  auto& instance = seissol::AcceleratorDevice::getInstance();
  instance.bindAcceleratorDevice(0);
#endif
}

void seissol::Mpi::printAcceleratorDeviceInfo() {
#ifdef ACL_DEVICE
  auto& instance = seissol::AcceleratorDevice::getInstance();
  instance.printInfo();

  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  const auto pci = device.api->getPciAddress(0);
  const auto pcisNode = collectContainer(pci, m_sharedMemComm);
  pcis = collectContainer(pci);
  logInfo() << "Device PCI address (rank=0): " << pci;
  logInfo() << "Device PCI addresses (node of rank=0):" << pcisNode;
#endif
}

void seissol::Mpi::setDataTransferModeFromEnv() {
  const auto envVariable =
      utils::Env("SEISSOL_").getOptional<std::string>("PREFERRED_MPI_DATA_TRANSFER_MODE");
  if (envVariable.has_value()) {
    std::string option{envVariable.value()};
    std::transform(option.begin(), option.end(), option.begin(), [](unsigned char c) {
      return std::tolower(c);
    });

    if (option == "direct") {
      preferredDataTransferMode = DataTransferMode::Direct;
    } else if (option == "host") {
      preferredDataTransferMode = DataTransferMode::CopyInCopyOutHost;
    } else {
      logWarning() << "Ignoring `SEISSOL_PREFERRED_MPI_DATA_TRANSFER_MODE`."
                   << "Expected values: direct, host.";
      option = "direct";
    }
#ifndef ACL_DEVICE
    if (preferredDataTransferMode != DataTransferMode::Direct) {
      logWarning() << "The CPU version of SeisSol supports"
                   << "only the `direct` MPI transfer mode.";
      option = "direct";
      preferredDataTransferMode = DataTransferMode::Direct;
    }
#endif
    logInfo() << "Selected" << option << "MPI data transfer mode as the preferred one";
  }
}

seissol::Mpi seissol::Mpi::mpi;
