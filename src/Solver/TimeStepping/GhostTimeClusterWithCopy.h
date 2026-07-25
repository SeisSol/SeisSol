// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_GHOSTTIMECLUSTERWITHCOPY_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_GHOSTTIMECLUSTERWITHCOPY_H_

#include "Parallel/Helper.h"
#include "Parallel/MPI.h"
#include "Solver/TimeStepping/AbstractGhostTimeCluster.h"
#include "Solver/TimeStepping/HaloCommunication.h"

#include <Device/device.h>

namespace seissol::time_stepping {
template <DataTransferMode CommType>
class GhostTimeClusterWithCopy : public AbstractGhostTimeCluster {
  public:
  GhostTimeClusterWithCopy(double maxTimeStepSize,
                           std::uint64_t timeStepRate,
                           std::size_t globalTimeClusterId,
                           std::size_t otherGlobalTimeClusterId,
                           const std::string& displayName,
                           const std::string& otherDisplayName,
                           const seissol::solver::HaloCommunication& meshStructure,
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

  std::list<std::size_t> prefetchCopyLayer();
  void prefetchGhostRegion(std::size_t region);

  private:
  std::vector<void*> duplicatedCopyRegions_;
  std::vector<void*> duplicatedGhostRegions_;

  std::vector<void*> prefetchCopyRegionsStreams_;
  std::vector<void*> prefetchGhostRegionsStreams_;

  enum class ReceiveState { RequiresMpiTesting, RequiresPrefetchTesting, Ready };
  std::vector<ReceiveState> receiveRegionsStates_;

  device::DeviceInstance& device_ = device::DeviceInstance::getInstance();

  bool persistent_;
};
} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_GHOSTTIMECLUSTERWITHCOPY_H_
