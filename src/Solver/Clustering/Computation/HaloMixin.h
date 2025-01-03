// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Solver/Clustering/AbstractTimeCluster.h>
#include <Solver/Clustering/ActorState.h>
#include <Solver/Clustering/Communication/NeighborCluster.h>
#include <type_traits>
namespace seissol::solver::clustering::computation {

// TODO: make mixin?
template <typename Base, ComputeStep SendStep, ComputeStep ArriveStep>
class CopyCluster : public Base {
  public:
  // TODO: constructors
  LayerType getLayerType() const override { return Copy; }

  protected:
  static_assert(std::is_base_of_v<AbstractTimeCluster, Base>,
                "Base needs to inherit from AbstractTimeCluster");

  void start() override { dataSent = false; }

  void runCompute(ComputeStep step) override {
    if (step == ArriveStep && dataSent) {
      neighbor->stopTo(this->streamRuntime);
    }

    Base::runCompute(step);

    if (step == SendStep) {
      neighbor->startFrom(this->streamRuntime);
      dataSent = true;
    }
  }

  void reset() override {
    if (dataSent) {
      neighbor->stopTo(this->streamRuntime);
      dataSent = false;
    }
  }

  private:
  std::shared_ptr<communication::SendNeighborCluster> neighbor;
  bool dataSent{false};
};

template <ComputeStep SendStep, ComputeStep ArriveStep>
class GhostCluster : public AbstractTimeCluster {
  public:
  GhostCluster(double maxTimeStepSize,
               long timeStepRate,
               std::shared_ptr<communication::RecvNeighborCluster> neighbor,
               const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor)
      : AbstractTimeCluster(maxTimeStepSize, timeStepRate, Executor::Host, cpuExecutor),
        neighbor(neighbor) {}

  LayerType getLayerType() const override { return Ghost; }

  protected:
  bool emptyStep(ComputeStep step) const override { return step != SendStep && step != ArriveStep; }

  void start() override {}

  void runCompute(ComputeStep step) override {
    if (step == SendStep) {
      neighbor->startFrom(streamRuntime);
    }
    if (step == ArriveStep) {
      neighbor->stopTo(streamRuntime);
    }
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

  private:
  std::shared_ptr<communication::RecvNeighborCluster> neighbor;
};

} // namespace seissol::solver::clustering::computation
