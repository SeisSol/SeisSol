// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Helper.h"

#include "Parallel/MPI.h"

#include <utils/env.h>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol {

void printCommThreadInfo(const Mpi& mpiBasic, utils::Env& env) {
  const bool useThread = env.get<bool>("COMMTHREAD", true);
  if (mpiBasic.isSingleProcess()) {
    logInfo() << "Using polling for advancing MPI communication, due to having only "
                 "one MPI rank running.";
  } else if (useThread) {
    logInfo() << "Using a communication thread for advancing MPI communication.";
  } else {
    logInfo() << "Using polling for advancing MPI communication.";
  }
}

bool useCommThread(const Mpi& mpiBasic, utils::Env& env) {
  const bool useThread = env.get<bool>("COMMTHREAD", true);
  return useThread && !mpiBasic.isSingleProcess();
}

bool usePersistentMpi(utils::Env& env) { return env.get<bool>("MPI_PERSISTENT", true); }

void printPersistentMpiInfo(utils::Env& env) {
  if (usePersistentMpi(env)) {
    logInfo() << "Using persistent MPI routines.";
  } else {
    logInfo() << "Using asynchronous MPI routines.";
  }
}

bool useUSM(SEISSOL_GPU_PARAM utils::Env& env) {
#ifdef ACL_DEVICE
  return env.get<bool>("USM", device::DeviceInstance::getInstance().api->isUnifiedMemoryDefault());
#else
  return true;
#endif
}

bool useUSM() {
  utils::Env env("SEISSOL_");
  return useUSM(env);
}

void printUSMInfo(utils::Env& env) {
  if (useUSM(env)) {
    logInfo() << "Using unified buffers for CPU-GPU data.";
  } else {
    logInfo() << "Using separate buffers for CPU-GPU data.";
  }
}

bool useMPIUSM(SEISSOL_GPU_PARAM utils::Env& env) {
#ifdef ACL_DEVICE
  return env.get<bool>("USM_MPI",
                       device::DeviceInstance::getInstance().api->isUnifiedMemoryDefault());
#else
  return true;
#endif
}

bool useMPIUSM() {
  utils::Env env("SEISSOL_");
  return useMPIUSM(env);
}

void printMPIUSMInfo(utils::Env& env) {
  if (useMPIUSM(env)) {
    logInfo() << "Using unified buffers for CPU-GPU MPI data.";
  } else {
    logInfo() << "Using separate buffers for CPU-GPU MPI data.";
  }
}

bool useDeviceL2Compress(utils::Env& env) { return env.get<bool>("L2_COMPRESS", false); }

bool useDeviceL2Compress() {
  utils::Env env("SEISSOL_");
  return useDeviceL2Compress(env);
}

void printDeviceL2Compress(utils::Env& env) {
  if (useDeviceL2Compress(env)) {
    logInfo() << "Using L2 compression (if available).";
  }
}

DataTransferMode getDataTransferMode(utils::Env& env) {
  DataTransferMode preferredDataTransferMode = DataTransferMode::Direct;

  const auto envVariable = env.getOptional<std::string>("TRANSFER_MODE");
  const auto envVariable2 = env.getOptional<std::string>("PREFERRED_MPI_DATA_TRANSFER_MODE");

  if (envVariable.has_value() || envVariable2.has_value()) {
    std::string option{envVariable.value_or(envVariable2.value_or("direct"))};
    std::transform(option.begin(), option.end(), option.begin(), [](unsigned char c) {
      return std::tolower(c);
    });

    if (option == "direct") {
      preferredDataTransferMode = DataTransferMode::Direct;
    } else if (option == "host") {
      preferredDataTransferMode = DataTransferMode::CopyInCopyOutHost;
    } else {
      logWarning() << "Ignoring `TRANSFER_MODE`."
                   << "Expected values: direct, host.";
      option = "direct";
    }
#ifndef ACL_DEVICE
    if (preferredDataTransferMode != DataTransferMode::Direct) {
      logWarning() << "The CPU version of SeisSol supports"
                   << "only the `direct` inter-process transfer mode.";
      option = "direct";
      preferredDataTransferMode = DataTransferMode::Direct;
    }
#endif
    logInfo() << "Selected" << option << "as inter-process data transfer mode.";
  }

  return preferredDataTransferMode;
}

} // namespace seissol
