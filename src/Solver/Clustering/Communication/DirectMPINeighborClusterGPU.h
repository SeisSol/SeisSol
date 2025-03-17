// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_DIRECTMPINEIGHBORCLUSTER_H_
#define SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_DIRECTMPINEIGHBORCLUSTER_H_

#include "NeighborCluster.h"
#include <Parallel/Runtime/Stream.h>
#include <mpi.h>
#include <vector>

namespace seissol::solver::clustering::communication {

class DirectMPISendNeighborClusterGPU : public SendNeighborCluster {
  public:
  bool poll() override;
  void start(parallel::runtime::StreamRuntime& runtime) override;
  void stop(parallel::runtime::StreamRuntime& runtime) override;

  DirectMPISendNeighborClusterGPU(const std::vector<RemoteCluster>& remote,
                                  const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
                                  double priority);
  ~DirectMPISendNeighborClusterGPU() override;

  private:
  std::vector<MPI_Request> requests;
  std::vector<int> status;
  std::mutex requestMutex;
  uint32_t progressRestart{0};
  uint32_t* progressEnd;
  uint32_t progressStart{0};
};

class DirectMPIRecvNeighborClusterGPU : public RecvNeighborCluster {
  public:
  bool poll() override;
  void start(parallel::runtime::StreamRuntime& runtime) override;
  void stop(parallel::runtime::StreamRuntime& runtime) override;

  DirectMPIRecvNeighborClusterGPU(const std::vector<RemoteCluster>& remote,
                                  const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
                                  double priority);
  ~DirectMPIRecvNeighborClusterGPU() override;

  private:
  std::vector<MPI_Request> requests;
  std::vector<int> status;
  std::mutex requestMutex;
  uint32_t progressRestart{0};
  uint32_t* progressEnd;
  uint32_t progressStart{0};
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
#endif // SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_DIRECTMPINEIGHBORCLUSTER_H_
