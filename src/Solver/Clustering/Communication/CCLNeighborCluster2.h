// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_CCLNEIGHBORCLUSTER2_H_
#define SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_CCLNEIGHBORCLUSTER2_H_

#include "NeighborCluster.h"
#include <DataTypes.h>
#include <Parallel/Runtime/Stream.h>
#include <atomic>
#include <mpi.h>
#include <vector>

namespace seissol::solver::clustering::communication {

class CCLNeighborCluster : public SendNeighborCluster {
  public:
  bool poll() override;
  void start(parallel::runtime::StreamRuntime& runtime) override;
  void stop(parallel::runtime::StreamRuntime& runtime) override;
  bool blocking() override { return false; }

  CCLNeighborCluster(const std::vector<RemoteCluster>& remoteSend,
                     const std::vector<RemoteCluster>& remoteRecv,
                     void* comm,
                     const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
                     double priority);
  ~CCLNeighborCluster() override;

  private:
  std::vector<RemoteCluster> remote;
  std::vector<bool> isSend;
  std::vector<void*> memoryHandles;
  void* comm;
  device::DeviceGraphHandle handle;
};

class CCLNeighborDummyCluster : public RecvNeighborCluster {
  public:
  bool poll() override { return true; }
  void start(parallel::runtime::StreamRuntime& runtime) override {}
  void stop(parallel::runtime::StreamRuntime& runtime) override {}
  bool blocking() override { return false; }

  CCLNeighborDummyCluster(const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
                          double priority)
      : RecvNeighborCluster(cpuExecutor, priority) {}
  ~CCLNeighborDummyCluster() override = default;
};

} // namespace seissol::solver::clustering::communication
#endif // SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_CCLNEIGHBORCLUSTER2_H_
