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
#include "Solver/TimeStepping/Halo/TaskHaloTransport.h"

#include <memory>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include "Solver/TimeStepping/Halo/Mpi/StagedMpiHaloTransport.h"
#endif

namespace seissol::solver {

HaloTransportFactory::HaloTransportFactory(Mpi::DataTransferMode mode, bool persistent)
    : mode_(mode), persistent_(persistent) {
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

std::unique_ptr<HaloTransport> HaloTransportFactory::create(const RemoteClusterPair& regions) {
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
  default:
    logError() << "The requested data transfer mode is not available in this build.";
    return nullptr;
  }
}

} // namespace seissol::solver
