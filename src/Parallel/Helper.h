// SPDX-FileCopyrightText: 2023 SeisSol Group
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
void printCommThreadInfo(const T& mpiBasic, utils::Env& env) {
  bool useThread = env.get<bool>("COMMTHREAD", true);
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
bool useCommThread(const T& mpiBasic, utils::Env& env) {
  bool useThread = env.get<bool>("COMMTHREAD", true);
  return useThread && !mpiBasic.isSingleProcess();
}

inline bool usePersistentMpi(utils::Env& env) { return env.get<bool>("MPI_PERSISTENT", true); }

inline void printPersistentMpiInfo(utils::Env& env) {
  if (usePersistentMpi(env)) {
    logInfo() << "Using persistent MPI routines.";
  } else {
    logInfo() << "Using asynchronous MPI routines.";
  }
}

#ifdef ACL_DEVICE
inline bool useUSM(utils::Env& env) {
  return env.get<bool>("USM", device::DeviceInstance::getInstance().api->isUnifiedMemoryDefault());
}

inline bool useUSM() {
  utils::Env env("SEISSOL_");
  return useUSM(env);
}

inline void printUSMInfo(utils::Env& env) {
  if (useUSM(env)) {
    logInfo() << "Using unified buffers for CPU-GPU data.";
  } else {
    logInfo() << "Using separate buffers for CPU-GPU data.";
  }
}

inline bool useMPIUSM(utils::Env& env) {
  return env.get<bool>("USM_MPI",
                       device::DeviceInstance::getInstance().api->isUnifiedMemoryDefault());
}

inline bool useMPIUSM() {
  utils::Env env("SEISSOL_");
  return useMPIUSM(env);
}

inline void printMPIUSMInfo(utils::Env& env) {
  if (useMPIUSM(env)) {
    logInfo() << "Using unified buffers for CPU-GPU MPI data.";
  } else {
    logInfo() << "Using separate buffers for CPU-GPU MPI data.";
  }
}
#endif

} // namespace seissol

#endif // SEISSOL_SRC_PARALLEL_HELPER_H_
