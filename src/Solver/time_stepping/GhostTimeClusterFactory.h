// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIME_STEPPING_GHOSTTIMECLUSTERFACTORY_H_
#define SEISSOL_SRC_SOLVER_TIME_STEPPING_GHOSTTIMECLUSTERFACTORY_H_

#include "Solver/time_stepping/DirectGhostTimeCluster.h"
#ifdef ACL_DEVICE
#include "Solver/time_stepping/GhostTimeClusterWithCopy.h"
#endif // ACL_DEVICE
#include "Parallel/MPI.h"
#include "memory"

namespace seissol::time_stepping {
struct GhostTimeClusterFactory {
  public:
  static std::unique_ptr<AbstractGhostTimeCluster> get(double maxTimeStepSize,
                                                       int timeStepRate,
                                                       int globalTimeClusterId,
                                                       int otherGlobalTimeClusterId,
                                                       const MeshStructure* meshStructure,
                                                       MPI::DataTransferMode mode,
                                                       bool persistent) {
    switch (mode) {
#ifdef ACL_DEVICE
    case MPI::DataTransferMode::CopyInCopyOutHost: {
      using ghostCluster_t = GhostTimeClusterWithCopy<MPI::DataTransferMode::CopyInCopyOutHost>;
      return std::make_unique<ghostCluster_t>(maxTimeStepSize,
                                              timeStepRate,
                                              globalTimeClusterId,
                                              otherGlobalTimeClusterId,
                                              meshStructure,
                                              persistent);
    }
#endif // ACL_DEVICE
    case MPI::DataTransferMode::Direct: {
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

#endif // SEISSOL_SRC_SOLVER_TIME_STEPPING_GHOSTTIMECLUSTERFACTORY_H_

