// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_HELPER_H_
#define SEISSOL_SRC_PARALLEL_HELPER_H_

#include "Common/Marker.h"

#include <utils/env.h>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol {
template <typename T>
void printCommThreadInfo(const T& mpiBasic, utils::Env& env) {
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

template <typename T>
bool useCommThread(const T& mpiBasic, utils::Env& env) {
  const bool useThread = env.get<bool>("COMMTHREAD", true);
  return useThread && !mpiBasic.isSingleProcess();
}

inline bool usePersistentMpi(utils::Env& env) { return env.get<bool>("MPI_PERSISTENT", true); }

/**
 * Whether the transfer modes that run on device streams (ccl, stream-mpi) use a stream per
 * direction between two time clusters, instead of one for all exchanges in a global order.
 */
inline bool useExchangePerDirection(utils::Env& env) {
  return env.get<bool>("EXCHANGE_PER_DIRECTION", env.get<bool>("CCL_PER_DIRECTION", false));
}

/**
 * Whether the clusters on the device only enqueue their work and wait for each other on the device,
 * instead of waiting for the completion of each of their steps. Implies scratchpads per layer.
 */
inline bool useConcurrentClusters(utils::Env& env) {
  return env.get<bool>("CONCURRENT_CLUSTERS", false);
}

/**
 * Whether the device work of regular super-timesteps gets recorded into graphs and replayed.
 * Requires the time stepping plan and concurrent clusters.
 */
inline bool useSuperStepGraphs(utils::Env& env) { return env.get<bool>("SUPERSTEP_GRAPHS", false); }

/**
 * Whether the clusters take their steps strictly along the time stepping plan, instead of whenever
 * they are ready.
 */
inline bool useTimeSteppingPlan(utils::Env& env) {
  return env.get<bool>("TIMESTEPPING_PLAN", false);
}

/**
 * Whether each layer gets scratchpads of its own, instead of all layers of a storage sharing one
 * set. Needed as soon as several layers are updated concurrently; costs the sum instead of the
 * maximum of the scratchpad demands.
 */
inline bool useScratchpadPerLayer(utils::Env& env) {
  return env.get<bool>("SCRATCHPAD_PER_LAYER", false);
}

inline void printPersistentMpiInfo(utils::Env& env) {
  if (usePersistentMpi(env)) {
    logInfo() << "Using persistent MPI routines.";
  } else {
    logInfo() << "Using asynchronous MPI routines.";
  }
}

inline bool useUSM(SEISSOL_GPU_PARAM utils::Env& env) {
#ifdef ACL_DEVICE
  return env.get<bool>("USM", device::DeviceInstance::getInstance().api->isUnifiedMemoryDefault());
#else
  return true;
#endif
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

inline bool useMPIUSM(SEISSOL_GPU_PARAM utils::Env& env) {
#ifdef ACL_DEVICE
  return env.get<bool>("USM_MPI",
                       device::DeviceInstance::getInstance().api->isUnifiedMemoryDefault());
#else
  return true;
#endif
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

inline bool useDeviceL2Compress(utils::Env& env) { return env.get<bool>("L2_COMPRESS", false); }

inline bool useDeviceL2Compress() {
  utils::Env env("SEISSOL_");
  return useDeviceL2Compress(env);
}

inline void printDeviceL2Compress(utils::Env& env) {
  if (useDeviceL2Compress(env)) {
    logInfo() << "Using L2 compression (if available).";
  }
}

} // namespace seissol

#endif // SEISSOL_SRC_PARALLEL_HELPER_H_
