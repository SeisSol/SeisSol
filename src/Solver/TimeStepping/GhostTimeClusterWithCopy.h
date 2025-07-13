// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_GHOSTTIMECLUSTERWITHCOPY_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_GHOSTTIMECLUSTERWITHCOPY_H_

#include "Parallel/MPI.h"
#include "Solver/TimeStepping/AbstractGhostTimeCluster.h"
#include <device.h>

namespace seissol::time_stepping {
template <MPI::DataTransferMode CommType>
class GhostTimeClusterWithCopy : public AbstractGhostTimeCluster {
  public:
  GhostTimeClusterWithCopy(double maxTimeStepSize,
                           int timeStepRate,
                           int globalTimeClusterId,
                           int otherGlobalTimeClusterId,
                           const MeshStructure* meshStructure,
                           bool persistent);
  ~GhostTimeClusterWithCopy() override;

  GhostTimeClusterWithCopy(const GhostTimeClusterWithCopy<CommType>&) = delete;
  GhostTimeClusterWithCopy(const GhostTimeClusterWithCopy<CommType>&&) = delete;
  GhostTimeClusterWithCopy& operator=(const GhostTimeClusterWithCopy<CommType>&) = delete;
  GhostTimeClusterWithCopy& operator=(const GhostTimeClusterWithCopy<CommType>&&) = delete;

  void sendCopyLayer() override;
  void receiveGhostLayer() override;

  bool testForGhostLayerReceives() override;
  bool testReceiveQueue();

  void finalize() override;

  std::list<int> prefetchCopyLayer();
  void prefetchGhostRegion(int region);

  private:
  unsigned int numberOfRegions{};
  std::vector<real*> duplicatedCopyRegions;
  std::vector<real*> duplicatedGhostRegions;

  std::vector<void*> prefetchCopyRegionsStreams;
  std::vector<void*> prefetchGhostRegionsStreams;

  enum class ReceiveState { RequiresMpiTesting, RequiresPrefetchTesting, Ready };
  std::vector<ReceiveState> receiveRegionsStates{};

  device::DeviceInstance& device = device::DeviceInstance::getInstance();

  bool persistent;
};
} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_GHOSTTIMECLUSTERWITHCOPY_H_
