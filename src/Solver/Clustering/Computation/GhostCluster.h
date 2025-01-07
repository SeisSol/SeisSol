// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Common/Executor.h>
#include <Kernels/Common.h>
#include <Solver/Clustering/AbstractTimeCluster.h>
#include <Solver/Clustering/ActorState.h>
#include <Solver/Clustering/Communication/NeighborCluster.h>
#include <memory>
namespace seissol::solver::clustering::computation {

class GhostCluster : public AbstractTimeCluster {
  public:
  GhostCluster(double maxTimeStepSize,
               long timeStepRate,
               std::shared_ptr<communication::RecvNeighborCluster> neighbor,
               const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
               double priority)
      : AbstractTimeCluster(maxTimeStepSize, timeStepRate, isDeviceOn() ? Executor::Device : Executor::Host, cpuExecutor, priority),
        neighbor(neighbor) {}

  LayerType getLayerType() const override { return Ghost; }

  protected:
  bool emptyStep(ComputeStep step) const override { return step != ComputeStep::Predict; }

  void start() override {}

  void runCompute(ComputeStep step) override {
    neighbor->startFrom(streamRuntime);
    neighbor->stopTo(streamRuntime);
  }

  void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) override {
    const auto rank = seissol::MPI::mpi.rank();
    logWarning(rank) << "No update since " << timeSinceLastUpdate.count()
                     << "[s] for global cluster " << 0 << " with local cluster id " << 0
                     << " at state " << actorStateToString(state)
                     << " predTime = " << ct.time.at(ComputeStep::Predict)
                     << " predictionsSinceSync = "
                     << ct.computeSinceLastSync.at(ComputeStep::Predict)
                     << " corrTime = " << ct.time.at(ComputeStep::Correct)
                     << " correctionsSinceSync = "
                     << ct.computeSinceLastSync.at(ComputeStep::Correct)
                     << " stepsTillSync = " << ct.stepsUntilSync
                     << " maySync = " << maySynchronize();
    for (auto& neighbor : neighbors) {
      logWarning(rank) << "Neighbor with rate = " << neighbor.ct.timeStepRate
                       << "PredTime = " << neighbor.ct.time.at(ComputeStep::Predict)
                       << "CorrTime = " << neighbor.ct.time.at(ComputeStep::Correct)
                       << "predictionsSinceSync = "
                       << neighbor.ct.computeSinceLastSync.at(ComputeStep::Predict)
                       << "correctionsSinceSync = "
                       << neighbor.ct.computeSinceLastSync.at(ComputeStep::Correct);
    }
  }

  std::string description() const override { return "ghost-cell"; }

  private:
  std::shared_ptr<communication::RecvNeighborCluster> neighbor;
};

} // namespace seissol::solver::clustering::computation
