// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "HaloTransportFactory.h"

#include "Parallel/MPI.h"
#include "Solver/TimeStepping/Halo/HaloCommunication.h"
#include "Solver/TimeStepping/Halo/HaloTransport.h"
#include "Solver/TimeStepping/Halo/Mpi/MpiHaloTransport.h"
#include "Solver/TimeStepping/Halo/Stream/CclExchangeScheduler.h"
#include "Solver/TimeStepping/Halo/Stream/ExchangeScheduler.h"
#include "Solver/TimeStepping/Halo/Stream/ShmemExchangeScheduler.h"
#include "Solver/TimeStepping/Halo/Stream/StreamMpiExchangeScheduler.h"

#include <cstddef>
#include <memory>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include "Solver/TimeStepping/Halo/Mpi/StagedMpiHaloTransport.h"
#endif

namespace seissol::time_stepping {

HaloTransportFactory::HaloTransportFactory(Mpi::DataTransferMode mode,
                                           bool persistent,
                                           bool perDirection,
                                           std::size_t clusterCount)
    : mode_(mode), persistent_(persistent) {
  const auto order = perDirection ? LaunchOrder::PerDirection : LaunchOrder::Global;
  if (mode_ == Mpi::DataTransferMode::DirectCcl) {
    scheduler_ = std::make_unique<CclExchangeScheduler>(clusterCount, order);
  } else if (mode_ == Mpi::DataTransferMode::DirectStreamMpi) {
    scheduler_ = std::make_unique<StreamMpiExchangeScheduler>(clusterCount, order);
  } else if (mode_ == Mpi::DataTransferMode::DirectShmem) {
    if (perDirection) {
      logWarning() << "The shmem transfer mode exchanges in the global order only.";
    }
    scheduler_ = std::make_unique<ShmemExchangeScheduler>(clusterCount);
  }
}

HaloTransportFactory::~HaloTransportFactory() = default;

std::unique_ptr<HaloTransport> HaloTransportFactory::create(
    const solver::RemoteClusterPair& regions, std::size_t cluster, std::size_t otherCluster) {
  switch (mode_) {
  case Mpi::DataTransferMode::Direct:
    return std::make_unique<MpiHaloTransport>(regions, persistent_);
#ifdef ACL_DEVICE
  case Mpi::DataTransferMode::CopyInCopyOutHost:
    return std::make_unique<StagedMpiHaloTransport>(regions, persistent_);
#endif
  case Mpi::DataTransferMode::DirectCcl:
  case Mpi::DataTransferMode::DirectStreamMpi:
  case Mpi::DataTransferMode::DirectShmem:
    return std::make_unique<ScheduledTransport>(*scheduler_, regions, cluster, otherCluster);
  default:
    logError() << "The requested data transfer mode is not available in this build.";
    return nullptr;
  }
}

} // namespace seissol::time_stepping
