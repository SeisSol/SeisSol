// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_HELPER_H_
#define SEISSOL_SRC_PARALLEL_HELPER_H_

#include "utils/env.h"
#include "utils/logger.h"

#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol {
template <typename T>
void printCommThreadInfo(const T& mpiBasic) {
  bool useThread = utils::Env::get<bool>("SEISSOL_COMMTHREAD", true);
  if (mpiBasic.isSingleProcess()) {
    logInfo() << "Using polling for advancing MPI communication, due to having only "
                 "one MPI rank running.";
  } else if (useThread) {
    logInfo() << "Using a communication thread for advancing MPI communication.";
  } else {
    logInfo() << "Using polling for advancing MPI communication.";
  }
}

template <typename T>
bool useCommThread(const T& mpiBasic) {
  bool useThread = utils::Env::get<bool>("SEISSOL_COMMTHREAD", true);
  return useThread && !mpiBasic.isSingleProcess();
}

inline bool usePersistentMpi() { return utils::Env::get<bool>("SEISSOL_MPI_PERSISTENT", true); }

inline void printPersistentMpiInfo() {
  if (usePersistentMpi()) {
    logInfo() << "Using persistent MPI routines.";
  } else {
    logInfo() << "Using asynchronous MPI routines.";
  }
}

#ifdef ACL_DEVICE
inline bool useUSM() {
  return utils::Env::get<bool>("SEISSOL_USM",
                               device::DeviceInstance::getInstance().api->isUnifiedMemoryDefault());
}

inline void printUSMInfo() {
  if (useUSM()) {
    logInfo() << "Using unified buffers for CPU-GPU data.";
  } else {
    logInfo() << "Using separate buffers for CPU-GPU data.";
  }
}

inline bool useMPIUSM() {
  return utils::Env::get<bool>("SEISSOL_USM_MPI",
                               device::DeviceInstance::getInstance().api->isUnifiedMemoryDefault());
}

inline void printMPIUSMInfo() {
  if (useMPIUSM()) {
    logInfo() << "Using unified buffers for CPU-GPU MPI data.";
  } else {
    logInfo() << "Using separate buffers for CPU-GPU MPI data.";
  }
}
#endif

} // namespace seissol

#endif // SEISSOL_SRC_PARALLEL_HELPER_H_
