// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_DIRECTMPINEIGHBORCLUSTERBLOCKING_H_
#define SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_DIRECTMPINEIGHBORCLUSTERBLOCKING_H_

#include "NeighborCluster.h"
#include <Parallel/Host/CpuExecutor.h>
#include <Parallel/Runtime/Stream.h>
#include <mpi.h>
#include <vector>

namespace seissol::solver::clustering::communication {

class DirectMPINeighborClusterBlocking {
  public:
  bool poll();
  void start(parallel::runtime::StreamRuntime& runtime);

  DirectMPINeighborClusterBlocking() = default;
  ~DirectMPINeighborClusterBlocking();

  protected:
  std::size_t enqueueCounter{0};
  std::size_t enqueueCounterTrue{0};
  std::size_t enqueueCounterPoll{0};
  std::size_t dequeueCounterPoll{0};
  std::vector<MPI_Request> requests;
  std::mutex requestMutex;
};

class DirectMPISendNeighborClusterBlocking : public SendNeighborCluster,
                                             DirectMPINeighborClusterBlocking {
  public:
  DirectMPISendNeighborClusterBlocking(
      const std::vector<RemoteCluster>& remote,
      const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
      double priority);
  ~DirectMPISendNeighborClusterBlocking() override = default;
  bool blocking() override { return true; }

  bool poll() override;
  void start(parallel::runtime::StreamRuntime& runtime) override;
  void stop(parallel::runtime::StreamRuntime& runtime) override;
};

class DirectMPIRecvNeighborClusterBlocking : public RecvNeighborCluster,
                                             DirectMPINeighborClusterBlocking {
  public:
  DirectMPIRecvNeighborClusterBlocking(
      const std::vector<RemoteCluster>& remote,
      const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
      double priority);
  ~DirectMPIRecvNeighborClusterBlocking() override = default;
  bool blocking() override { return true; }

  bool poll() override;
  void start(parallel::runtime::StreamRuntime& runtime) override;
  void stop(parallel::runtime::StreamRuntime& runtime) override;
};

/*
class CopyMPISendNeighborCluster : public SendNeighborCluster {
public:
    CopyMPISendNeighborCluster(const std::vector<RemoteCluster>& remote);
    ~CopyMPISendNeighborCluster() override;

    bool poll() override;
    void start(parallel::runtime::StreamRuntime& runtime) override;
private:
    std::vector<MPI_Request> requests;
    std::vector<void*> copyStreams;
    std::mutex requestMutex;
};

class CopyMPIRecvNeighborCluster : public RecvNeighborCluster {
public:
    CopyMPIRecvNeighborCluster(const std::vector<RemoteCluster>& remote);
    ~CopyMPIRecvNeighborCluster() override;

    bool poll() override;
    void start(parallel::runtime::StreamRuntime& runtime) override;
private:
    std::vector<MPI_Request> requests;
    std::vector<void*> copyStreams;
    std::mutex requestMutex;
};
*/

} // namespace seissol::solver::clustering::communication
#endif // SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_DIRECTMPINEIGHBORCLUSTERBLOCKING_H_
