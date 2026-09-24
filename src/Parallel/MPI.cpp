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
  if (gethostname(hostName.data(), 256) != 0) {
    hostName = "unknown-host";
  } else {
    utils::StringUtils::rtrim(hostName);
    hostName.pop_back();
  }
  hostNames_ = collectContainer(hostName);

  // Test this after setComm() to get the correct rank_
  if (provided < required) {
    logError() << utils::nospace << "Provided MPI thread support (" << provided
               << ") is smaller than required thread support (" << required << ").";
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
  const auto pcisNode = collectContainer(pci, sharedMemComm_);
  pcis_ = collectContainer(pci);
  logInfo() << "Device PCI address (rank=0): " << pci;
  logInfo() << "Device PCI addresses (node of rank=0):" << pcisNode;
#endif
}

void seissol::Mpi::setDataTransferModeFromEnv() {
  auto env = utils::Env("SEISSOL_");
  const auto envVariable = env.getOptional<std::string>("TRANSFER_MODE");
  const auto envVariableMpi = env.getOptional<std::string>("PREFERRED_MPI_DATA_TRANSFER_MODE");
  if (envVariable.has_value() || envVariableMpi.has_value()) {
    std::string option{envVariable.value_or(envVariableMpi.value_or("direct"))};
    std::transform(option.begin(), option.end(), option.begin(), [](unsigned char c) {
      return std::tolower(c);
    });

    if (option == "direct") {
      preferredDataTransferMode_ = DataTransferMode::Direct;
    } else if (option == "host") {
      preferredDataTransferMode_ = DataTransferMode::CopyInCopyOutHost;
    } else if (option == "ccl") {
      preferredDataTransferMode_ = DataTransferMode::DirectCcl;
    } else {
      logWarning() << "Ignoring `SEISSOL_TRANSFER_MODE`."
                   << "Expected values: direct, host, ccl.";
      option = "direct";
    }
#ifndef ACL_DEVICE
    if (preferredDataTransferMode_ != DataTransferMode::Direct) {
      logWarning() << "The CPU version of SeisSol supports"
                   << "only the `direct` transfer mode.";
      option = "direct";
      preferredDataTransferMode_ = DataTransferMode::Direct;
    }
#endif
#ifndef USE_CCL
    if (preferredDataTransferMode_ == DataTransferMode::DirectCcl) {
      logWarning() << "This build of SeisSol does not support the `ccl` transfer mode.";
      option = "direct";
      preferredDataTransferMode_ = DataTransferMode::Direct;
    }
#endif
    logInfo() << "Selected" << option << "as data transfer mode between processes.";
  }
}

seissol::Mpi seissol::Mpi::mpi;
