// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_NEIGHBORCLUSTER_H_
#define SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_NEIGHBORCLUSTER_H_

#include <Initializer/TimeStepping/ClusterLayout.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Common.h>
#include <Parallel/Host/CpuExecutor.h>
#include <Parallel/Runtime/Stream.h>
#include <Solver/Clustering/AbstractTimeCluster.h>
#include <Solver/Clustering/ActorState.h>

#include "Parallel/MPI.h"

namespace seissol::solver::clustering::communication {

struct RemoteCluster {
  void* data;
  std::size_t size;
  MPI_Datatype datatype;
  int rank;
  int tag;
};

struct HaloCommunication {
  std::vector<std::vector<RemoteCluster>> ghost;
  std::vector<std::vector<RemoteCluster>> copy;
};

HaloCommunication getHaloCommunication(const initializer::ClusterLayout& layout,
                                       const MeshStructure* structure);

class NeighborCluster : public AbstractTimeCluster {
  public:
  NeighborCluster(double maxTimeStepSize,
                  long timeStepRate,
                  const std::vector<RemoteCluster>& sends,
                  const std::vector<RemoteCluster>& receives,
                  const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
                  double priority)
      : AbstractTimeCluster(maxTimeStepSize,
                            timeStepRate,
                            isDeviceOn() ? Executor::Device : Executor::Host,
                            cpuExecutor,
                            priority) {}

  void start() override {}
  bool emptyStep(ComputeStep step) const override { return step != ComputeStep::Communicate; }

  virtual bool poll() { return true; }

  ~NeighborCluster() override = default;

  void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) override {
    logWarning(true) << "No update since " << timeSinceLastUpdate.count()
                     << "[s] for global cluster at rate " << timeStepRate << " at state "
                     << actorStateToString(state)
                     << " predTime = " << ct.time.at(ComputeStep::Predict)
                     << " predictionsSinceSync = "
                     << ct.computeSinceLastSync.at(ComputeStep::Predict)
                     << " corrTime = " << ct.time.at(ComputeStep::Correct)
                     << " correctionsSinceSync = "
                     << ct.computeSinceLastSync.at(ComputeStep::Correct)
                     << " stepsTillSync = " << ct.stepsUntilSync
                     << " maySync = " << maySynchronize();
    for (auto& neighbor : neighbors) {
      logWarning(true) << "Neighbor with rate = " << neighbor.ct.timeStepRate
                       << "PredTime = " << neighbor.ct.time.at(ComputeStep::Predict)
                       << "CorrTime = " << neighbor.ct.time.at(ComputeStep::Correct)
                       << "predictionsSinceSync = "
                       << neighbor.ct.computeSinceLastSync.at(ComputeStep::Predict)
                       << "correctionsSinceSync = "
                       << neighbor.ct.computeSinceLastSync.at(ComputeStep::Correct);
    }
  }

  std::string description() const override { return "communication"; }

  protected:
  std::vector<RemoteCluster> sends;
  std::vector<RemoteCluster> receives;
};

} // namespace seissol::solver::clustering::communication
#endif // SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_NEIGHBORCLUSTER_H_
