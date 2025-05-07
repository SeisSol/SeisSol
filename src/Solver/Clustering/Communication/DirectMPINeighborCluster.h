// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_DIRECTMPINEIGHBORCLUSTER_H_
#define SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_DIRECTMPINEIGHBORCLUSTER_H_

#include "NeighborCluster.h"
#include <Parallel/Host/CpuExecutor.h>
#include <Parallel/Runtime/Stream.h>
#include <Solver/Clustering/ActorState.h>
#include <mpi.h>
#include <vector>

#include "utils/logger.h"

namespace seissol::solver::clustering::communication {

class DirectMPINeighborClusterUnidirectional {
  public:
  bool poll();
  void start(parallel::runtime::StreamRuntime& runtime);

  DirectMPINeighborClusterUnidirectional(const std::vector<MPI_Request>& requests)
      : requests(requests) {}
  ~DirectMPINeighborClusterUnidirectional();

  protected:
  std::size_t enqueueCounter{0};
  std::size_t enqueueCounterTrue{0};
  std::size_t enqueueCounterPoll{0};
  std::size_t dequeueCounterPoll{0};
  std::vector<MPI_Request> requests;
  std::mutex requestMutex;
};

class DirectMPINeighborCluster : public NeighborCluster {
  private:
  static std::vector<MPI_Request> initSends(const std::vector<RemoteCluster>& sends) {
    std::vector<MPI_Request> sendRequests(sends.size());
    for (std::size_t i = 0; i < sends.size(); ++i) {
      MPI_Send_init(sends[i].data,
                    sends[i].size,
                    MPI::precisionToMpiType(sends[i].datatype),
                    sends[i].rank,
                    sends[i].tag,
                    MPI::mpi.comm(),
                    &sendRequests[i]);
    }
    return sendRequests;
  }

  static std::vector<MPI_Request> initReceives(const std::vector<RemoteCluster>& receives) {
    std::vector<MPI_Request> recvRequests(receives.size());
    for (std::size_t i = 0; i < receives.size(); ++i) {
      MPI_Recv_init(receives[i].data,
                    receives[i].size,
                    MPI::precisionToMpiType(receives[i].datatype),
                    receives[i].rank,
                    receives[i].tag,
                    MPI::mpi.comm(),
                    &recvRequests[i]);
    }
    return recvRequests;
  }

  public:
  DirectMPINeighborCluster(double maxTimeStepSize,
                           long timeStepRate,
                           const std::vector<RemoteCluster>& sends,
                           const std::vector<RemoteCluster>& receives,
                           const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
                           double priority)
      : NeighborCluster(maxTimeStepSize, timeStepRate, sends, receives, cpuExecutor, priority),
        send(initSends(sends)), recv(initReceives(receives)) {}

  ~DirectMPINeighborCluster() override = default;

  bool emptyStep(ComputeStep step) const override { return step == ComputeStep::Interact; }

  void runCompute(ComputeStep step) override {
    if (step == ComputeStep::Predict) {
      recv.start(streamRuntime);
    }
    if (step == ComputeStep::Communicate) {
      send.start(streamRuntime);
    }
  }

  bool pollCompute(ComputeStep step) override {
    const auto sendPollResult = send.poll();
    const auto recvPollResult = recv.poll();
    if (step == ComputeStep::Communicate) {
      return recvPollResult;
    }
    if (step == ComputeStep::Correct) {
      return sendPollResult;
    }
    return true;
  }

  bool poll() override {
    const auto sendPollResult = send.poll();
    const auto recvPollResult = recv.poll();
    return recvPollResult && sendPollResult;
  }

  DirectMPINeighborClusterUnidirectional send;
  DirectMPINeighborClusterUnidirectional recv;
};

} // namespace seissol::solver::clustering::communication
#endif // SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_DIRECTMPINEIGHBORCLUSTER_H_
