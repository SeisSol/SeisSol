#pragma once

#include "Solver/Clustering/Communication/DirectGhostTimeCluster.h"
#include <Solver/Clustering/AbstractTimeCluster.h>
#ifdef ACL_DEVICE
#include "Solver/Clustering/Communication/GhostTimeClusterWithCopy.h"
#endif // ACL_DEVICE
#include "Parallel/MPI.h"
#include "memory"

#if defined(ACL_DEVICE) && defined(USE_CCL)
#include "CCLCluster.hpp"
#include "CCLSetup.hpp"
#endif

namespace seissol::time_stepping {
struct GhostTimeClusterFactory {
  private:
  static inline std::vector<void*> comms;
  static inline std::size_t globalClusterCount;

  public:
  static void prepare(MPI::DataTransferMode mode, std::size_t clusterCount) {
    globalClusterCount = clusterCount;
#ifdef ACL_DEVICE
    if (mode == MPI::DataTransferMode::DirectCCL) {
      const auto commCount = (clusterCount * (clusterCount + 1)) / 2;
      comms = seissol::solver::clustering::communication::createComms(commCount);
    }
#endif
  }

  static std::unique_ptr<AbstractTimeCluster> get(double maxTimeStepSize,
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
#ifdef USE_CCL
    case MPI::DataTransferMode::DirectCCL: {
      using ghostCluster_t = CCLCluster;
      auto baseCluster = std::min(globalTimeClusterId, otherGlobalTimeClusterId);
      auto otherCluster = std::max(globalTimeClusterId, otherGlobalTimeClusterId);
      auto offset = (baseCluster * (baseCluster + 1)) / 2;
      return std::make_unique<ghostCluster_t>(maxTimeStepSize,
                                              timeStepRate,
                                              globalTimeClusterId,
                                              otherGlobalTimeClusterId,
                                              meshStructure,
                                              comms[offset + otherCluster],
                                              comms[offset + otherCluster]);
    }
#endif
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
