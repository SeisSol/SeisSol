// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_CCLNEIGHBORCLUSTER_H_
#define SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_CCLNEIGHBORCLUSTER_H_

#include "NeighborCluster.h"
#include <Parallel/Runtime/Stream.h>
#include <atomic>
#include <mpi.h>
#include <vector>

namespace seissol::solver::clustering::communication {

class CCLSendNeighborCluster : public SendNeighborCluster {
  public:
  bool poll() override;
  void start(parallel::runtime::StreamRuntime& runtime) override;
  void stop(parallel::runtime::StreamRuntime& runtime) override;

  CCLSendNeighborCluster(const std::vector<RemoteCluster>& remote,
                         void* comm,
                         const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
                         double priority);
  ~CCLSendNeighborCluster() override;

  private:
  std::vector<RemoteCluster> remote;
  std::vector<void*> memoryHandles;
  void* comm;
};

class CCLRecvNeighborCluster : public RecvNeighborCluster {
  public:
  bool poll() override;
  void start(parallel::runtime::StreamRuntime& runtime) override;
  void stop(parallel::runtime::StreamRuntime& runtime) override;

  CCLRecvNeighborCluster(const std::vector<RemoteCluster>& remote,
                         void* comm,
                         const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
                         double priority);
  ~CCLRecvNeighborCluster() override;

  private:
  std::vector<RemoteCluster> remote;
  std::vector<void*> memoryHandles;
  void* comm;
};

} // namespace seissol::solver::clustering::communication
#endif // SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_CCLNEIGHBORCLUSTER_H_
