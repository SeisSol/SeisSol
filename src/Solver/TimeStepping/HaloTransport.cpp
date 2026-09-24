// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "HaloTransport.h"

#include "Parallel/MPI.h"
#include "Solver/TimeStepping/HaloCommunication.h"
#include "Solver/TimeStepping/MpiHaloTransport.h"

#include <memory>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include "Solver/TimeStepping/StagedMpiHaloTransport.h"
#endif

namespace seissol::time_stepping {

std::unique_ptr<HaloTransport> createHaloTransport(const solver::RemoteClusterPair& regions,
                                                   Mpi::DataTransferMode mode,
                                                   bool persistent) {
  switch (mode) {
  case Mpi::DataTransferMode::Direct:
    return std::make_unique<MpiHaloTransport>(regions, persistent);
#ifdef ACL_DEVICE
  case Mpi::DataTransferMode::CopyInCopyOutHost:
    return std::make_unique<StagedMpiHaloTransport>(regions, persistent);
#endif
  default:
    logError() << "The requested MPI data transfer mode is not available in this build.";
    return nullptr;
  }
}

} // namespace seissol::time_stepping
