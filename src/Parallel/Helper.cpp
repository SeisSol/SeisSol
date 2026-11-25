// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Helper.h"

#include "Parallel/MPI.h"

#include <algorithm>
#include <cctype>
#include <optional>
#include <string>
#include <string_view>
#include <utils/env.h>
#include <utils/logger.h>

namespace seissol {

void printCommThreadInfo(const Mpi& mpi, utils::Env& env) {
  const bool useThread = commThreadRequested(env);
  if (mpi.isSingleProcess()) {
    logInfo() << "Using polling for advancing MPI communication, due to having only "
                 "one MPI rank running.";
  } else if (useThread) {
    logInfo() << "Using a communication thread for advancing MPI communication.";
  } else {
    logInfo() << "Using polling for advancing MPI communication.";
  }
}

bool useCommThread(const Mpi& mpi, utils::Env& env) {
  const bool useThread = commThreadRequested(env);
  return useThread && !mpi.isSingleProcess();
}

std::string_view dataTransferModeName(DataTransferMode mode) {
  switch (mode) {
  case DataTransferMode::Direct:
    return "direct";
  case DataTransferMode::CopyInCopyOutHost:
    return "host";
  }
  return "unknown";
}

DataTransferMode getDataTransferMode(utils::Env& env) {
  const auto envVariable = env.getOptional<std::string>("TRANSFER_MODE");
  const auto envVariableOld = env.getOptional<std::string>("PREFERRED_MPI_DATA_TRANSFER_MODE");
  if (!envVariable.has_value() && !envVariableOld.has_value()) {
    return DataTransferMode::Direct;
  }

  std::string option{envVariable.value_or(envVariableOld.value_or("direct"))};
  std::ranges::transform(option, option.begin(), [](unsigned char c) { return std::tolower(c); });

  auto mode = DataTransferMode::Direct;
  if (option == dataTransferModeName(DataTransferMode::CopyInCopyOutHost)) {
    mode = DataTransferMode::CopyInCopyOutHost;
  } else if (option != dataTransferModeName(DataTransferMode::Direct)) {
    logWarning() << "Ignoring `SEISSOL_TRANSFER_MODE`."
                 << "Expected values: direct, host.";
  }
  if (!DeviceBuild && mode != DataTransferMode::Direct) {
    logWarning() << "The CPU version of SeisSol supports"
                 << "only the `direct` transfer mode.";
    mode = DataTransferMode::Direct;
  }
  logInfo() << "Selected" << std::string(dataTransferModeName(mode))
            << "as data transfer mode between processes.";
  return mode;
}

} // namespace seissol
