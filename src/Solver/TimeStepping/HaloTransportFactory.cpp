// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "HaloTransportFactory.h"

#include "Kernels/Common.h"
#include "Parallel/MPI.h"
#include "Parallel/Runtime/HostRuntime.h"
#include "Solver/TimeStepping/Halo/HaloCommunication.h"
#include "Solver/TimeStepping/Halo/HaloTransport.h"
#include "Solver/TimeStepping/Halo/Mpi/MpiHaloTransport.h"
#include "Solver/TimeStepping/Halo/SplitHaloTransport.h"
#include "Solver/TimeStepping/Halo/Stream/CclExchangeScheduler.h"
#include "Solver/TimeStepping/Halo/Stream/ExchangeScheduler.h"
#include "Solver/TimeStepping/Halo/Stream/ShmemExchangeScheduler.h"
#include "Solver/TimeStepping/Halo/Stream/StreamMpiExchangeScheduler.h"
#include "Solver/TimeStepping/Halo/TaskHaloTransport.h"

#include <cstddef>
#include <memory>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include "Solver/TimeStepping/Halo/Mpi/StagedMpiHaloTransport.h"
#endif

namespace seissol::solver {

HaloTransportFactory::HaloTransportFactory(Mpi::DataTransferMode mode,
                                           bool persistent,
                                           bool perDirection,
                                           std::size_t clusterCount,
                                           std::size_t configCount,
                                           std::size_t lanes)
    : mode_(mode), persistent_(persistent) {
  const auto order = perDirection ? LaunchOrder::PerDirection : LaunchOrder::Global;
  if (lanes > 1 && (perDirection || mode_ == Mpi::DataTransferMode::DirectShmem)) {
    logWarning() << "Lanes of the halo exchange (SEISSOL_EXCHANGE_LANES) are only for ccl and"
                 << "stream-mpi in the global order; the exchange takes none.";
    lanes = 1;
  }
  if (mode_ == Mpi::DataTransferMode::DirectCcl) {
    scheduler_ = std::make_unique<CclExchangeScheduler>(clusterCount, order, configCount, lanes);
  } else if (mode_ == Mpi::DataTransferMode::DirectStreamMpi) {
    scheduler_ =
        std::make_unique<StreamMpiExchangeScheduler>(clusterCount, order, configCount, lanes);
  } else if (mode_ == Mpi::DataTransferMode::DirectShmem) {
    if (perDirection) {
      logWarning() << "The shmem transfer mode exchanges in the global order only.";
    }
    scheduler_ = std::make_unique<ShmemExchangeScheduler>(clusterCount, configCount);
  }
#ifdef SEISSOL_HOST_TASK_GRAPH
  // (with a device, the task graph only takes the work on the host for the device, see
  // parallel::runtime::HostActor; the clusters there do not wait for each other in it)
  if (parallel::runtime::HostRuntime::scheduling() ==
          parallel::runtime::HostScheduling::TaskGraph &&
      !isDeviceOn()) {
    taskEngine_ = std::make_unique<TaskHaloEngine>();
  }
#endif
}

HaloTransportFactory::~HaloTransportFactory() = default;

void HaloTransportFactory::splitByNode(std::vector<bool> onNode) { onNode_ = std::move(onNode); }

std::unique_ptr<HaloTransport> HaloTransportFactory::create(const RemoteClusterPair& regions,
                                                            std::size_t cluster,
                                                            std::size_t otherCluster) {
  switch (mode_) {
  case Mpi::DataTransferMode::Direct:
#ifdef SEISSOL_HOST_TASK_GRAPH
    if (taskEngine_ != nullptr) {
      // the exchanges become tasks of the task graph on the host
      return std::make_unique<TaskHaloTransport>(
          std::make_unique<MpiHaloTransport>(regions, persistent_), *taskEngine_);
    }
#endif
    return std::make_unique<MpiHaloTransport>(regions, persistent_);
#ifdef ACL_DEVICE
  case Mpi::DataTransferMode::CopyInCopyOutHost:
    return std::make_unique<StagedMpiHaloTransport>(regions, persistent_);
#endif
  case Mpi::DataTransferMode::DirectCcl:
  case Mpi::DataTransferMode::DirectStreamMpi:
  case Mpi::DataTransferMode::DirectShmem:
    if (split()) {
      // within the node on the device streams, to the other nodes from the host
      auto [onNode, offNode] = SplitHaloTransport::split(regions, onNode_);
      return std::make_unique<SplitHaloTransport>(
          std::make_unique<ScheduledTransport>(*scheduler_, onNode, cluster, otherCluster),
          std::make_unique<MpiHaloTransport>(offNode, persistent_));
    }
    return std::make_unique<ScheduledTransport>(*scheduler_, regions, cluster, otherCluster);
  default:
    logError() << "The requested data transfer mode is not available in this build.";
    return nullptr;
  }
}

} // namespace seissol::solver
