// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_GHOSTTIMECLUSTERFACTORY_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_GHOSTTIMECLUSTERFACTORY_H_

#include "Solver/TimeStepping/DirectGhostTimeCluster.h"
#include "Solver/TimeStepping/HaloCommunication.h"
#ifdef ACL_DEVICE
#include "Solver/TimeStepping/GhostTimeClusterWithCopy.h"
#ifdef USE_CCL
#include "Solver/TimeStepping/CCLNeighborCluster.h"
#endif
#ifdef USE_SHMEM
#include "Solver/TimeStepping/ShmemCluster.h"
#endif
#endif // ACL_DEVICE
#include "Parallel/MPI.h"
#include "memory"

namespace seissol::time_stepping {
struct GhostTimeClusterFactory {
  public:
  static std::vector<void*> setup(Mpi::DataTransferMode mode, std::size_t clusters) {
    if (mode == Mpi::DataTransferMode::DirectCcl) {
      const auto commCount = (clusters * (clusters + 1)) / 2;

#if defined(ACL_DEVICE) && defined(USE_CCL)
      return createComms(commCount);
#endif
    }

    return {};
  }

  static std::unique_ptr<AbstractTimeCluster> get(double maxTimeStepSize,
                                                  int timeStepRate,
                                                  int globalTimeClusterId,
                                                  int otherGlobalTimeClusterId,
                                                  const solver::HaloCommunication& meshStructure,
                                                  Mpi::DataTransferMode mode,
                                                  const std::vector<void*>& comms,
                                                  bool persistent) {
    switch (mode) {
#ifdef ACL_DEVICE
    case Mpi::DataTransferMode::CopyInCopyOutHost: {
      using GhostClusterT = GhostTimeClusterWithCopy<Mpi::DataTransferMode::CopyInCopyOutHost>;
      return std::make_unique<GhostClusterT>(maxTimeStepSize,
                                             timeStepRate,
                                             globalTimeClusterId,
                                             otherGlobalTimeClusterId,
                                             meshStructure,
                                             persistent);
    }
#endif // ACL_DEVICE
#if defined(ACL_DEVICE) && defined(USE_CCL)
    case Mpi::DataTransferMode::DirectCcl: {
      return std::make_unique<CCLNeighborCluster>(maxTimeStepSize,
                                                  timeStepRate,
                                                  globalTimeClusterId,
                                                  otherGlobalTimeClusterId,
                                                  meshStructure,
                                                  persistent,
                                                  comms);
    }
#endif // defined(ACL_DEVICE) && defined(USE_CCL)
#if defined(ACL_DEVICE) && defined(USE_SHMEM)
    case Mpi::DataTransferMode::DirectShmem: {
      return std::make_unique<ShmemCluster>(maxTimeStepSize,
                                            timeStepRate,
                                            globalTimeClusterId,
                                            otherGlobalTimeClusterId,
                                            meshStructure,
                                            persistent);
    }
#endif // defined(ACL_DEVICE) && defined(USE_SHMEM)
    case Mpi::DataTransferMode::Direct: {
      return std::make_unique<DirectGhostTimeCluster>(maxTimeStepSize,
                                                      timeStepRate,
                                                      globalTimeClusterId,
                                                      otherGlobalTimeClusterId,
                                                      meshStructure,
                                                      persistent);
    }
    default: {
      return nullptr;
    }
    }
  }
};
} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_GHOSTTIMECLUSTERFACTORY_H_
