// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "HaloTransport.h"

#include "Parallel/MPI.h"
#include "Solver/TimeStepping/CclExchangeScheduler.h"
#include "Solver/TimeStepping/ExchangeScheduler.h"
#include "Solver/TimeStepping/HaloCommunication.h"
#include "Solver/TimeStepping/MpiHaloTransport.h"

#include <cstddef>
#include <memory>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include "Solver/TimeStepping/StagedMpiHaloTransport.h"
#endif

namespace seissol::time_stepping {

HaloTransportFactory::HaloTransportFactory(Mpi::DataTransferMode mode,
                                           bool persistent,
                                           std::size_t clusterCount)
    : mode_(mode), persistent_(persistent) {
  if (mode_ == Mpi::DataTransferMode::DirectCcl) {
    scheduler_ = std::make_unique<CclExchangeScheduler>(clusterCount);
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
    return std::make_unique<ScheduledTransport>(*scheduler_, regions, cluster, otherCluster);
  default:
    logError() << "The requested data transfer mode is not available in this build.";
    return nullptr;
  }
}

} // namespace seissol::time_stepping
