// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_CCLNEIGHBORCLUSTER_H_
#define SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_CCLNEIGHBORCLUSTER_H_

#include "NeighborCluster.h"
#include <DataTypes.h>
#include <Parallel/Runtime/Stream.h>
#include <atomic>
#include <mpi.h>
#include <vector>

namespace seissol::solver::clustering::communication {

class CCLNeighborCluster : public NeighborCluster {
  public:
  CCLNeighborCluster(double maxTimeStepSize,
                     long timeStepRate,
                     const std::vector<RemoteCluster>& sends,
                     const std::vector<RemoteCluster>& receives,
                     void* comm,
                     const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
                     double priority);

  ~CCLNeighborCluster() override;

  bool emptyStep(ComputeStep step) const override { return step != ComputeStep::Communicate; }

  void runCompute(ComputeStep step) override;

  private:
  std::vector<RemoteCluster> remote;
  std::vector<bool> isSend;
  std::vector<void*> memoryHandles;
  void* comm;
  device::DeviceGraphHandle handle;
};

} // namespace seissol::solver::clustering::communication
#endif // SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_CCLNEIGHBORCLUSTER_H_
