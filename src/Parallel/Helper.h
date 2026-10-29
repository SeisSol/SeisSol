// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_HELPER_H_
#define SEISSOL_SRC_PARALLEL_HELPER_H_

#include "Common/Marker.h"

#include <algorithm>
#include <cstddef>
#include <string>
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
 * Whether the transfer modes that run on device streams (ccl, stream-mpi, shmem) exchange only
 * within a node on them, and with the processes on other nodes with GPU-aware MPI from the host, as
 * the direct mode does (experimental).
 */
inline bool useExchangeSplitByNode(utils::Env& env) {
  return env.get<bool>("EXCHANGE_SPLIT_BY_NODE", false);
}

/**
 * Into how many lanes the transfer modes that run on device streams (ccl, stream-mpi) split their
 * global order (experimental): streams of the highest priority, which run at the same time, with
 * the groups of each in the global order.
 */
inline std::size_t exchangeLanes(utils::Env& env) {
  return std::max<std::size_t>(env.get<std::size_t>("EXCHANGE_LANES", 1), 1);
}

/**
 * On how many queues of their own the clusters run their work on the device, 0 for none (see
 * device::AbstractAPI::enableOwnQueues(); HIP only): the other stream runtimes share one more.
 */
inline int ownQueueCount(utils::Env& env) { return std::max(env.get<int>("OWN_QUEUES", 0), 0); }

/**
 * Whether the clusters on the device only enqueue their work and wait for each other on the device,
 * instead of waiting for the completion of each of their steps. Implies scratchpads per layer.
 */
inline bool useConcurrentClusters(utils::Env& env) {
  return env.get<bool>("CONCURRENT_CLUSTERS", false);
}

/**
 * Whether concurrent clusters wait for each other by what their work does with the shared data:
 * before reading, for the latest writes of the neighbors, before writing, for their latest reads;
 * instead of for their latest actions of any kind.
 */
inline bool useAccessEvents(utils::Env& env) { return env.get<bool>("ACCESS_EVENTS", true); }

/**
 * Whether the device work of regular super-timesteps gets recorded into graphs and replayed.
 * Requires the time stepping plan and concurrent clusters.
 */
inline bool useSuperStepGraphs(utils::Env& env) { return env.get<bool>("SUPERSTEP_GRAPHS", false); }

/**
 * Whether the work on the host that reads what the device has written (the samples of the
 * receivers and fault receivers) runs on the thread that enqueues the work of the device once the
 * device is done with it (see parallel::runtime::HostActor), instead of in host functions on
 * streams; also for recorded super-timesteps, after each of their replays.
 */
inline bool useHostActor(utils::Env& env) { return env.get<bool>("HOST_ACTOR", false); }

/**
 * Whether the work that depends on the halo exchange on device streams waits only for the groups
 * of its own transport, instead of all groups launched so far (experimental; one stream in the
 * global order only).
 */
inline bool useExchangeTransportEvents(utils::Env& env) {
  return env.get<bool>("EXCHANGE_TRANSPORT_EVENTS", false);
}

/**
 * Whether the independent branches of the work of a cluster (e.g. the neighbor fluxes of the face
 * relations of a face) become nodes of their own in the recordings of super-timesteps, so that they
 * can run concurrently. Otherwise, they run one after another, as one node: the graphs get fewer
 * nodes, which makes them quicker to record and to launch (with HIP, each launch costs the host
 * time per node).
 */
inline bool useSuperStepBranches(utils::Env& env) {
  return env.get<bool>("SUPERSTEP_BRANCHES", true);
}

/**
 * Which super-timesteps get recorded into graphs with SEISSOL_SUPERSTEP_GRAPHS=1 (experimental).
 */
enum class SuperStepGraphMode {
  //! the regular ones, the second time they come up, to replay them from then on
  Replay,
  //! besides, all others that can be, each launched once right after its recording
  Always,
  //! every one that can be, anew, and launched once; none is replayed (to measure the recordings)
  Rebuild
};

inline SuperStepGraphMode superStepGraphMode(utils::Env& env) {
  const auto mode = env.get<std::string>("SUPERSTEP_GRAPH_MODE", "replay");
  if (mode == "always") {
    return SuperStepGraphMode::Always;
  }
  if (mode == "rebuild") {
    return SuperStepGraphMode::Rebuild;
  }
  if (mode != "replay") {
    logWarning() << "Unknown SEISSOL_SUPERSTEP_GRAPH_MODE" << mode
                 << "; known are replay, always and rebuild. Taking replay.";
  }
  return SuperStepGraphMode::Replay;
}

/**
 * Whether the outputs process their synchronization points on a thread of their own, while the
 * device goes on with the next epoch (experimental).
 */
inline bool useOutputThread(utils::Env& env) { return env.get<bool>("OUTPUT_THREAD", false); }

/**
 * Whether the thread of the outputs copies their data to the host from snapshots on the device,
 * while the device goes on with the next epoch (see useOutputThread()).
 */
inline bool useOutputSnapshot(utils::Env& env) { return env.get<bool>("OUTPUT_SNAPSHOT", true); }

/**
 * Whether the storages keep snapshots of the data that changes (see
 * initializer::DualMemoryContainer::snapshot): for the thread of the outputs, unless
 * SEISSOL_OUTPUT_SNAPSHOT=0, with all clusters on the device.
 */
inline bool keepOutputSnapshots([[maybe_unused]] utils::Env& env) {
#ifdef ACL_DEVICE
  return useOutputThread(env) && useOutputSnapshot(env) &&
         env.get<std::string>("DEVICE_HOST_SWITCH", "none") == "none";
#else
  return false;
#endif
}

/**
 * Whether the clusters run their work on the host in a task graph (see
 * parallel::runtime::HostScheduling), instead of each loop in a parallel region of its own.
 */
inline bool useHostTaskGraph(utils::Env& env) { return env.get<bool>("HOST_TASK_GRAPH", false); }

/**
 * The cells or faces that each task of a loop takes in a task graph on the host.
 */
inline std::size_t hostTaskSize(utils::Env& env) {
  return env.get<std::size_t>("HOST_TASK_SIZE", 32);
}

/**
 * How many actions may wait for their work in a task graph on the host at most, before the main
 * thread waits for them to catch up (0 for no limit).
 */
inline std::size_t hostActionLimit(utils::Env& env) {
  return env.get<std::size_t>("HOST_ACTION_LIMIT", 256);
}

/**
 * Whether the halo exchange gets tested before the time stepping (see solver::testHaloExchange()).
 */
inline bool useHaloCheck(utils::Env& env) { return env.get<bool>("HALO_CHECK", true); }

/**
 * Whether the loops that the clusters enqueue on the device get timed for the loop statistics.
 */
inline bool useDeviceLoopStatistics(utils::Env& env) {
  return env.get<bool>("DEVICE_LOOP_STATISTICS", true);
}

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

/**
 * Kill switch for explicit graph node construction. When off, graphs that would be built node
 * by node fall back to whole-stream capture, which is what every graph used before.
 */
inline bool useGraphNodes(SEISSOL_GPU_PARAM utils::Env& env) {
#ifdef ACL_DEVICE
  return env.get<bool>("DEVICE_GRAPH_NODES", true) &&
         device::DeviceInstance::instance().api().isCapableOfGraphNodes();
#else
  return false;
#endif
}

inline bool useGraphNodes() {
  utils::Env env("SEISSOL_");
  return useGraphNodes(env);
}

/**
 * Reports how the compute graphs get built on a device; and the options that only concern devices
 * when there is none.
 */
inline void printDeviceOptionInfo(utils::Env& env) {
#ifdef ACL_DEVICE
  if (useGraphNodes(env)) {
    logInfo() << "Building compute graphs node by node where possible.";
  } else {
    logInfo() << "Capturing all compute graphs from streams.";
  }
#else
  if (useConcurrentClusters(env)) {
    logWarning() << "SEISSOL_CONCURRENT_CLUSTERS has no effect without a device.";
  }
  if (useScratchpadPerLayer(env)) {
    logWarning() << "SEISSOL_SCRATCHPAD_PER_LAYER has no effect without a device.";
  }
#endif
}

inline void printPersistentMpiInfo(utils::Env& env) {
  if (usePersistentMpi(env)) {
    logInfo() << "Using persistent MPI routines.";
  } else {
    logInfo() << "Using asynchronous MPI routines.";
  }
}

/**
 * Whether the host and the device share the buffers of the data (unified memory), instead of each
 * having its own. Off by default, also on devices that share their memory with the host (e.g. Grace
 * Hopper, MI300A, APUs): separate buffers keep what the host reads (e.g. for the outputs) apart
 * from what the device works on, and leave the device its own memory.
 */
inline bool useUSM(SEISSOL_GPU_PARAM utils::Env& env) {
#ifdef ACL_DEVICE
  return env.get<bool>("USM", false);
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
#ifdef ACL_DEVICE
    if (device::DeviceInstance::instance().api().isUnifiedMemoryDefault()) {
      logInfo() << "The device shares its memory with the host; SEISSOL_USM=1 would save the copies"
                << "the host keeps of the data.";
    }
#endif
  }
}

/**
 * Whether the buffers that are also used for MPI transfers are unified memory; off by default, as
 * useUSM().
 */
inline bool useMPIUSM(SEISSOL_GPU_PARAM utils::Env& env) {
#ifdef ACL_DEVICE
  return env.get<bool>("USM_MPI", false);
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
