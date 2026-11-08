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
#include <optional>
#include <string>
#include <string_view>
#include <utils/env.h>
#include <utils/logger.h>

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

  std::string hostName(MPI_MAX_PROCESSOR_NAME, '\0');
  int hostNameLength = 0;
  if (MPI_Get_processor_name(hostName.data(), &hostNameLength) != MPI_SUCCESS) {
    hostName = "unknown-host";
  } else {
    hostName.resize(hostNameLength);
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
  if (outputComm_ != MPI_COMM_NULL) {
    MPI_Comm_free(&outputComm_);
  }
  MPI_Comm_dup(comm, &outputComm_);

  MPI_Comm_rank(comm, &rank_);
  MPI_Comm_size(comm, &size_);

  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &sharedMemComm_);
  MPI_Comm_rank(sharedMemComm_, &sharedMemMpiRank_);
  MPI_Comm_size(sharedMemComm_, &sharedMemMpiSize_);
}

void seissol::Mpi::bindAcceleratorDevice() {
#ifdef ACL_DEVICE
  auto& instance = seissol::AcceleratorDevice::instance();
  instance.bindAcceleratorDevice(0);
#endif
}

void seissol::Mpi::printAcceleratorDeviceInfo() {
#ifdef ACL_DEVICE
  auto& instance = seissol::AcceleratorDevice::instance();
  instance.printInfo();

  const device::DeviceInstance& device = device::DeviceInstance::instance();
  const auto pci = device.api().getPciAddress(0);
  const auto pcisNode = collectContainer(pci, sharedMemComm_);
  pcis_ = collectContainer(pci);
  logInfo() << "Device PCI address (rank=0): " << pci;
  logInfo() << "Device PCI addresses (node of rank=0):" << pcisNode;
#endif
}

namespace {
/// the value of SEISSOL_TRANSFER_MODE (or of its older name), in lower case
std::optional<std::string> transferModeOption() {
  auto env = utils::Env("SEISSOL_");
  const auto envVariable = env.getOptional<std::string>("TRANSFER_MODE");
  const auto envVariableMpi = env.getOptional<std::string>("PREFERRED_MPI_DATA_TRANSFER_MODE");
  if (!envVariable.has_value() && !envVariableMpi.has_value()) {
    return std::nullopt;
  }
  std::string option{envVariable.value_or(envVariableMpi.value_or("direct"))};
  std::ranges::transform(option, option.begin(), [](unsigned char c) { return std::tolower(c); });
  return option;
}

std::optional<seissol::Mpi::DataTransferMode> transferModeOf(const std::string& option) {
  using DataTransferMode = seissol::Mpi::DataTransferMode;
  for (const auto mode : {DataTransferMode::Direct, DataTransferMode::CopyInCopyOutHost}) {
    if (option == seissol::Mpi::dataTransferModeName(mode)) {
      return mode;
    }
  }
  return std::nullopt;
}
} // namespace

std::optional<seissol::Mpi::DataTransferMode> seissol::Mpi::requestedDataTransferMode() {
  const auto option = transferModeOption();
  return option.has_value() ? transferModeOf(option.value()) : std::nullopt;
}

void seissol::Mpi::setDataTransferModeFromEnv() {
  if (auto chosen = transferModeOption(); chosen.has_value()) {
    std::string option = chosen.value();
    const auto mode = transferModeOf(option);
    if (mode.has_value()) {
      preferredDataTransferMode_ = mode.value();
    } else {
      logWarning() << "Ignoring `SEISSOL_TRANSFER_MODE`."
                   << "Expected values: direct, host.";
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
    logInfo() << "Selected" << option << "as data transfer mode between processes.";
  }
}

std::string_view seissol::Mpi::dataTransferModeName(DataTransferMode mode) {
  switch (mode) {
  case DataTransferMode::Direct:
    return "direct";
  case DataTransferMode::CopyInCopyOutHost:
    return "host";
  }
  return "unknown";
}

seissol::Mpi seissol::Mpi::mpi;
